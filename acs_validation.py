"""ACS validation (needs CENSUS_API_KEY in the environment; the key is never printed or written).
Compares modeled transit share by origin zone (real/proposed network, observed LODES flows) with ACS 2022 5-yr transit commute share (B08301_010E/B08301_001E)
aggregated from tracts to zones via block-centroid crosswalk. Run: python3 -I acs_validation.py"""
import os, sys, json, math, urllib.request, urllib.parse
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np, pandas as pd
import cities, model

KEY = os.environ.get("CENSUS_API_KEY")
if not KEY and os.path.exists(os.path.join(ROOT, ".census_key")):      # git-ignored local file, never committed
    KEY = open(os.path.join(ROOT, ".census_key")).read().strip()
if not KEY:
    sys.exit("No Census key (CENSUS_API_KEY or .census_key); skipping ACS validation.")
STATES = {"baltimore_lodes": ("24", "md", ["510", "005", "003", "027", "025"]), "triangle_lodes": ("37", "nc", ["183", "063", "135", "037", "101"])}
REAL = {"baltimore_lodes": ["Metro SubwayLink", "Light RailLink", "Red Line (planned)"], "triangle_lodes": ["Proposed commuter rail (Durham-Garner-Clayton)"]}
out = {}
for name, (fips, st, counties) in STATES.items():
    rows = []
    for cty in counties:
        q = urllib.parse.urlencode({"get": "B08301_001E,B08301_010E", "for": "tract:*", "in": f"state:{fips} county:{cty}", "key": KEY})
        d = json.load(urllib.request.urlopen("https://api.census.gov/data/2022/acs/acs5?" + q, timeout=60))
        rows += [dict(tract=r[2] + r[3] + r[4], tot=float(r[0]), transit=float(r[1])) for r in d[1:]]
    acs = pd.DataFrame(rows).set_index("tract")
    c = cities.get_city(name)
    xw = pd.read_csv(os.path.join(ROOT, "data/raw", st, "xwalk.csv.gz"), usecols=["trct", "blklatdd", "blklondd"], dtype={"trct": str})
    x = c_ax = None
    z = np.load(os.path.join(ROOT, "data", name.split("_")[0] + "_lodes.npz"))
    X = float(z["ax"]) + (xw.blklondd.values - float(z["lon0"])) * 111.32 * math.cos(math.radians(float(z["lat0"])))
    Y = float(z["ay"]) - (xw.blklatdd.values - float(z["lat0"])) * 110.57
    xw["r"], xw["c"] = np.floor(Y / c.cell).astype(int), np.floor(X / c.cell).astype(int)
    xw = xw[(xw.r >= 0) & (xw.r < c.R) & (xw.c >= 0) & (xw.c < c.C)]
    xw["zone"] = (xw.r // c.block) * (-(-c.C // c.block)) + xw.c // c.block
    tz = xw.groupby(["trct", "zone"]).size().reset_index(name="n")
    tz = tz.join(acs, on="trct").dropna()
    w = tz.n / tz.groupby("trct").n.transform("sum")
    tz["tot"] *= w; tz["transit"] *= w
    zs = tz.groupby("zone")[["tot", "transit"]].sum()
    zs["acs_share"] = zs.transit / zs.tot
    # modeled origin-zone transit share (same computation as model.evaluate, keeping the origin dimension)
    rn = c.real_networks(); L = [l for k in REAL[name] for l in rn[k]]
    fleet = model.fleet_for(c.mid_budget)
    idx, hw, T, pred, nid = model.station_times(c, L, 5., fleet)
    S, W = model._access(c, L, idx)
    Bm = (W[:, :, None] + T[S]).min(1); Cm = (Bm[:, S] + W[None]).min(-1)
    fin = np.isfinite(Cm)
    p = np.where(fin, 1 / (1 + np.exp(-(c.beta * np.where(fin, c.drive_t - Cm, 0) + c.asc_z[:, None]))), 0.0)
    zshare = (c.flow * p).sum(1) / c.flow.sum(1).clip(1e-9)
    Cz = -(-c.C // c.block)
    zone_id = (np.floor(c.zy / c.cell / c.block).astype(int)) * Cz + np.floor(c.zx / c.cell / c.block).astype(int)
    mod = pd.DataFrame({"zone": zone_id, "modeled_share": zshare, "workers": c.zH}).set_index("zone")
    j = zs.join(mod, how="inner")
    j = j[(j.tot > 200) & (j.workers > 200)]
    wc = lambda a, b, w: float(np.cov(a, b, aweights=w)[0, 1] / math.sqrt(np.cov(a, aweights=w) * np.cov(b, aweights=w)))
    out[name] = dict(n_zones=int(len(j)), acs_overall_transit_share=float(zs.transit.sum() / zs.tot.sum()), modeled_mode_share_all_trips=float((c.flow * p).sum() / c.total_trips),
                     modeled_mean_origin_share=float(np.average(j.modeled_share, weights=j.workers)), acs_mean_zone_share=float(np.average(j.acs_share, weights=j.tot)),
                     pearson=float(np.corrcoef(j.acs_share, j.modeled_share)[0, 1]), weighted_pearson=wc(j.acs_share.values, j.modeled_share.values, j.workers.values),
                     spearman=float(j.acs_share.corr(j.modeled_share, method="spearman")))
    j.to_csv(os.path.join(ROOT, "results", f"acs_zones_{name}.csv"))
json.dump(out, open(os.path.join(ROOT, "results", "acs_validation.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
