"""v2 robustness analysis: 10 seeds, sensitivity variants, convergence check, post-hoc structure metrics."""
import sys, os, json, glob, math
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np, pandas as pd, networkx as nx
import cities, model
from analysis import boot_ci

BUD, PEN = ["low", "mid", "high"], ["lowpen", "highpen"]


def load_v1():
    rows = []
    for p in glob.glob(os.path.join(ROOT, "results", "runs", "*.json")):
        r = json.load(open(p)); rows.append(dict({k: r[k] for k in ("city", "budget", "pen", "mode", "seed")}, variant="base", gens=500,
                                                  lines=r["lines"], phase_fit=r["phase_fit"], **r["metrics"]))
    return rows


def load_v2(suite):
    rows = []
    for p in glob.glob(os.path.join(ROOT, "results", "v2", suite, "*.json")):
        r = json.load(open(p)); rows.append(dict({k: r[k] for k in ("city", "budget", "pen", "mode", "seed", "variant", "gens")},
                                                  lines=r["lines"], phase_fit=r["phase_fit"], **r["metrics"]))
    return rows


def paired(df, factor, a, b, metric, keys, rel=False):
    A = df[df[factor] == a].set_index(keys)[metric]; B = df[df[factor] == b].set_index(keys)[metric]
    j = pd.concat([A.rename("a"), B.rename("b")], axis=1).dropna()
    d = (j.b - j.a) / j.a.abs().clip(lower=1e-9) if rel else (j.b - j.a)
    lo, hi = boot_ci(d.values)
    return dict(n=len(d), mean=float(d.mean()), ci=(lo, hi), dz=float(d.mean() / d.std(ddof=1)) if d.std() > 0 else float("nan"))


def headline(df):
    """Headline effects for any subset of runs (all variants identical keys)."""
    out = {}
    k = ["city", "budget", "mode", "seed"]
    for m in ["n_lines", "mean_stations_per_line", "n_transfer_stations"]:
        pm = df.groupby("pen")[m].mean()
        out[f"pen_{m}_ratio"] = float(pm["highpen"] / pm["lowpen"])
    out["pen_fitness_rel"] = paired(df, "pen", "lowpen", "highpen", "fitness", k, rel=True)["mean"]
    k3 = ["city", "budget", "pen", "seed"]
    A = df[df["mode"] == "planned"].set_index(k3).fitness; B = df[df["mode"] == "phased"].set_index(k3).fitness
    g = (1 - B / A).dropna().reset_index(name="gap")
    out["phased_gap_mean"] = float(g.gap.mean()); out["phased_gap_median"] = float(g.gap.median())
    out["phased_gap_frac_worse"] = float((g.gap > 0).mean())
    for b in BUD:
        out[f"phased_gap_{b}"] = float(g[g.budget == b].gap.mean())
    out["gap_poly_minus_mono"] = float(g[g.city == "poly"].gap.mean() - g[g.city == "mono"].gap.mean())
    out["gap_river_minus_mono"] = float(g[g.city == "river"].gap.mean() - g[g.city == "mono"].gap.mean())
    pl = df[df["mode"] == "planned"]
    for c in ["mono", "baltimore", "triangle"]:
        for b in BUD:
            out[f"modeshare_{c}_{b}"] = float(pl[(pl.city == c) & (pl.budget == b)].mode_share.mean())
    return out


def real_ratio(variant):
    """Evolved (planned, mid, mean over seeds in df) vs real/proposed network under the variant's parameters."""
    res = {}
    for cn, rk in (("baltimore", ["Metro SubwayLink", "Light RailLink", "Red Line (planned)"]), ("triangle", ["Proposed commuter rail (Durham-Garner-Clayton)"])):
        c = cities.get_variant_city(cn, variant)
        rn = cities.get_city(cn).real_networks()
        L = [l for k in rk for l in rn[k]]
        fl = model.fleet_for(cities.get_city(cn).mid_budget, c.fleet_per_unit)
        res[cn] = {p: model.evaluate(c, L, pv, fl) for p, pv in (("lowpen", 5.0), ("highpen", 20.0))}
    return res


def structure(lines, city):
    G = nx.Graph()
    for ln in lines:
        for a, b in zip(ln[:-1], ln[1:]):
            G.add_edge(tuple(a), tuple(b))
    xy = {n: city.xy(*n) for n in G}
    rings = 0
    for cyc in nx.cycle_basis(G):
        pts = np.array([xy[n] for n in cyc])
        area = 0.5 * abs(np.dot(pts[:, 0], np.roll(pts[:, 1], 1)) - np.dot(pts[:, 1], np.roll(pts[:, 0], 1)))
        if area >= 6.0:
            rings += 1
    cnt = {}
    for ln in lines:
        for c in set(map(tuple, ln)):
            cnt[c] = cnt.get(c, 0) + 1
    hubs = [n for n, v in cnt.items() if v >= 2 or G.degree(n) >= 3]
    out = dict(rings_area6=rings, core_frac=np.nan, branch_mean_len=np.nan, tip_over_core_dist=np.nan, n_branches=np.nan)
    comp = max(nx.connected_components(G), key=len)
    H = G.subgraph(comp)
    hubs = [h for h in hubs if h in H]
    if len(hubs) >= 2:
        core = set(hubs)
        for a in hubs:
            for b in hubs:
                if a < b:
                    core.update(nx.shortest_path(H, a, b))
    elif len(hubs) == 1:
        core = {hubs[0]}
    else:
        return out
    R = H.copy(); R.remove_nodes_from(core)
    branches = list(nx.connected_components(R))
    cen = np.array([np.average([xy[n][i] for n in core]) for i in (0, 1)])
    out["core_frac"] = len(core) / H.number_of_nodes()
    out["n_branches"] = len(branches)
    out["branch_mean_len"] = float(np.mean([len(b) for b in branches])) if branches else 0.0
    if branches:
        tips = [n for b in branches for n in b if H.degree(n) == 1]
        if tips:
            out["tip_over_core_dist"] = float(np.mean([np.hypot(*(np.array(xy[n]) - cen)) for n in tips]) /
                                               max(np.mean([np.hypot(*(np.array(xy[n]) - cen)) for n in core]), 0.5))
    return out


def main():
    S = {}
    v1 = pd.DataFrame(load_v1()); s2 = pd.DataFrame(load_v2("seeds"))
    base10 = pd.concat([v1, s2], ignore_index=True)
    S["n_base10"] = len(base10)
    S["base3_headline"] = headline(v1)
    S["base10_headline"] = headline(base10)
    # seed noise with 10 seeds
    pl = base10[base10["mode"] == "planned"]
    cv = pl.groupby(["city", "budget", "pen"]).fitness.agg(lambda x: x.std() / x.mean())
    S["base10_planned_cv"] = dict(median=float(cv.median()), max=float(cv.max()))
    # Baltimore / Triangle ratio vs real, 10 seeds
    rr = real_ratio("base")
    S["real_vs_evolved_base10"] = {}
    for cn in rr:
        for p in PEN:
            f = pl[(pl.city == cn) & (pl.budget == "mid") & (pl.pen == p)].fitness
            S["real_vs_evolved_base10"][f"{cn}|{p}"] = dict(real=rr[cn][p], evolved_mean=float(f.mean()), ratio=float(f.mean() / rr[cn][p]),
                                                           ratio_min=float(f.min() / rr[cn][p]), n=int(len(f)))
    # convergence
    lg = pd.DataFrame(load_v2("long"))
    if len(lg):
        k = ["city", "budget", "pen", "seed"]
        a = v1[(v1["mode"] == "planned") & (v1.pen == "lowpen")].set_index(k).fitness
        b = lg.set_index(k).fitness
        r = (b / a).dropna()
        S["long_vs_500"] = dict(n=int(len(r)), mean_ratio=float(r.mean()), median_ratio=float(r.median()), min=float(r.min()), max=float(r.max()),
                                frac_improved=float((r > 1.0005).mean()),
                                by_city=r.groupby(level="city").mean().round(4).to_dict(),
                                cv_1500=float(lg.groupby(["city", "budget"]).fitness.agg(lambda x: x.std() / x.mean()).median()),
                                cv_500_same=float(v1[(v1["mode"] == "planned") & (v1.pen == "lowpen")].groupby(["city", "budget"]).fitness.agg(lambda x: x.std() / x.mean()).median()))
    # sensitivity
    sn = pd.DataFrame(load_v2("sens"))
    sens = {"base(seeds0-1)": headline(v1[v1.seed < 2])}
    if len(sn):
        for v in sorted(sn.variant.unique()):
            sub = sn[sn.variant == v]
            h = headline(sub)
            rr = real_ratio(v)
            for cn in rr:
                for p in PEN:
                    f = sub[(sub.city == cn) & (sub.budget == "mid") & (sub["mode"] == "planned") & (sub.pen == p)].fitness
                    h[f"ratio_vs_real_{cn}_{p}"] = float(f.mean() / rr[cn][p]) if rr[cn][p] > 0 else float("nan")
                    h[f"real_fit_{cn}_{p}"] = float(rr[cn][p])
            sens[v] = h
    rb = real_ratio("base")
    for cn in rb:
        for p in PEN:
            f = v1[(v1.seed < 2) & (v1.city == cn) & (v1.budget == "mid") & (v1["mode"] == "planned") & (v1.pen == p)].fitness
            sens["base(seeds0-1)"][f"ratio_vs_real_{cn}_{p}"] = float(f.mean() / rb[cn][p])
            sens["base(seeds0-1)"][f"real_fit_{cn}_{p}"] = float(rb[cn][p])
    S["sensitivity"] = sens
    # post-hoc structure on planned base10 networks
    rows = []
    for _, r in pl.iterrows():
        c = cities.get_city(r.city)
        st = structure([[tuple(x) for x in l] for l in r.lines], c)
        rows.append(dict(city=r.city, budget=r.budget, pen=r.pen, seed=r.seed, **st))
    st = pd.DataFrame(rows)
    st.to_csv(os.path.join(ROOT, "results", "structure_posthoc.csv"), index=False)
    S["structure_by_city_budget"] = st.groupby(["city", "budget"])[["rings_area6", "core_frac", "branch_mean_len", "tip_over_core_dist", "n_branches"]].mean().round(3).reset_index().to_dict("records")
    S["structure_frac_tips_beyond_core"] = float((st.tip_over_core_dist.dropna() > 1).mean())
    S["structure_median_tip_over_core"] = float(st.tip_over_core_dist.median())
    json.dump(S, open(os.path.join(ROOT, "results", "stats_v2.json"), "w"), indent=1, default=float)
    return S


if __name__ == "__main__":
    S = main()
    print(json.dumps({k: S[k] for k in S if k in ("n_base10", "base3_headline", "base10_headline", "base10_planned_cv", "real_vs_evolved_base10", "long_vs_500")}, indent=1, default=float)[:6000])
