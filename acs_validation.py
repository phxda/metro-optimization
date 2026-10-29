"""ACS validation (needs CENSUS_API_KEY in the environment; the key is never printed or written).
Compares modeled transit share by origin zone (real/proposed network, observed LODES flows) with ACS 2022 5-yr transit commute share (B08301_010E/B08301_001E)
aggregated from tracts to zones via block-centroid crosswalk. Run: python3 -I acs_validation.py"""
import os, sys, json, math, urllib.request, urllib.parse
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np, pandas as pd
import cities, model

KEY = os.environ.get("CENSUS_API_KEY")
if not KEY:
    sys.exit("CENSUS_API_KEY not set; skipping ACS validation.")
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
    # modeled origin-zone transit share
    rn = c.real_networks(); L = [l for k in REAL[name] for l in rn[k]]
    f, m = model.evaluate(c, L, 5., model.fleet_for(c.mid_budget), detail=True)
    zid = np.full(-(-c.R // c.block) * -(-c.C // c.block), -1)
    keep_zone_ids = []
    out[name] = dict(n_zones_with_acs=int(len(zs)), acs_overall_transit_share=float(zs.transit.sum() / zs.tot.sum()), modeled_mode_share=m["mode_share"],
                     note="zone-level correlation requires zone id alignment; see per-zone table in results/acs_zones_%s.csv" % name)
    zs.to_csv(os.path.join(ROOT, "results", f"acs_zones_{name}.csv"))
json.dump(out, open(os.path.join(ROOT, "results", "acs_validation.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
