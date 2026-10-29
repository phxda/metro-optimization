"""Build results.csv, figures/, stats.json from results/runs/*.json."""
import sys, os, json, glob, math
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
import pandas as pd
import cities, model, render
import matplotlib.pyplot as plt

BUD = ["low", "mid", "high"]
PEN = ["lowpen", "highpen"]
SYN = ["mono", "coastal", "river", "poly"]
NAMES = {"mono": "Monocentric", "coastal": "Linear coastal", "river": "River-split", "poly": "Polycentric",
         "baltimore": "Baltimore*", "triangle": "Triangle*"}
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "figure.dpi": 130})


def load():
    rows = []
    for p in sorted(glob.glob(os.path.join(ROOT, "results", "runs", "*.json"))):
        r = json.load(open(p))
        m = r["metrics"]
        row = {k: r[k] for k in ("id", "city", "budget", "pen", "mode", "seed", "budget_units", "fleet", "penalty",
                                 "evals", "secs", "fit_cross_penalty")}
        row.update(m)
        row["budget_used_frac"] = m["cost"] / r["budget_units"]
        for i, v in enumerate(r["phase_fit"]):
            row[f"phase{i+1}_fitness"] = v
        row["lines_json"] = json.dumps(r["lines"])
        rows.append(row)
    return pd.DataFrame(rows)


def boot_ci(x, n=4000, seed=0):
    x = np.asarray(x, float)
    if len(x) < 2:
        return (float("nan"),) * 2
    rng = np.random.default_rng(seed)
    b = rng.choice(x, (n, len(x))).mean(1)
    return float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def paired(df, factor, a, b, metric, keys, rel=False):
    """Paired differences (b - a) of `metric` across runs that match on `keys`. Returns summary dict."""
    A = df[df[factor] == a].set_index(keys)[metric]
    B = df[df[factor] == b].set_index(keys)[metric]
    j = pd.concat([A.rename("a"), B.rename("b")], axis=1).dropna()
    d = (j.b - j.a) / j.a.abs().clip(lower=1e-9) if rel else (j.b - j.a)
    sd = d.std(ddof=1)
    lo, hi = boot_ci(d.values)
    return dict(n=int(len(d)), mean=float(d.mean()), median=float(d.median()), ci_lo=lo, ci_hi=hi,
                dz=float(d.mean() / sd) if sd > 0 else float("nan"), frac_pos=float((d > 0).mean()),
                frac_neg=float((d < 0).mean()))


def seed_consistency(df, factor, a, b, metric, cond_keys):
    """Fraction of conditions where all 3 seeds agree on sign of (b-a) (paired by seed)."""
    keys = cond_keys + ["seed"]
    A = df[df[factor] == a].set_index(keys)[metric]
    B = df[df[factor] == b].set_index(keys)[metric]
    d = (B - A).dropna().reset_index(name="d")
    g = d.groupby(cond_keys)["d"]
    allpos = g.apply(lambda s: (s > 0).all()).mean()
    allneg = g.apply(lambda s: (s < 0).all()).mean()
    return float(allpos), float(allneg)


def main():
    df = load()
    df.drop(columns=["lines_json"]).to_csv(os.path.join(ROOT, "results.csv"), index=False)
    S = {}
    cond = ["city", "budget", "pen", "seed"]

    # ---------- Q1: topology vs budget and shape (planned, both penalties) ----------
    pl = df[df["mode"] == "planned"]
    g = pl.groupby(["city", "budget", "pen"])[["n_lines", "mean_stations_per_line", "n_transfer_stations", "transfers_per_line",
                                              "rings", "mode_share", "fitness", "water_hops", "longest_share", "track_km",
                                              "n_stations", "mean_transfers_per_trip"]].mean()
    S["q1_table"] = g.reset_index().round(3).to_dict("records")
    # P1a line-count ratio high/low budget
    r = {}
    for c in cities.CITY_ORDER:
        for p in PEN:
            lo_ = pl[(pl.city == c) & (pl.budget == "low") & (pl.pen == p)].n_lines.mean()
            hi_ = pl[(pl.city == c) & (pl.budget == "high") & (pl.pen == p)].n_lines.mean()
            r[f"{c}|{p}"] = hi_ / lo_
    S["P1a_line_ratio_high_over_low"] = r
    # P1b rings per city/budget
    S["P1b_rings_mono"] = {b: pl[(pl.city == "mono") & (pl.budget == b)].groupby("pen").rings.apply(list).to_dict() for b in BUD}
    S["P1c_coastal"] = {b: dict(angle=pl[(pl.city == "coastal") & (pl.budget == b)].dominant_angle.round(1).tolist(),
                                rings=pl[(pl.city == "coastal") & (pl.budget == b)].rings.tolist()) for b in BUD}
    S["P1d_river_water_hops"] = {b: pl[(pl.city == "river") & (pl.budget == b)].groupby("pen").water_hops.apply(list).to_dict() for b in BUD}
    S["P1e_transfers_per_line_syn"] = {b: pl[(pl.city.isin(SYN)) & (pl.budget == b)].groupby("city").transfers_per_line.mean().round(3).to_dict() for b in BUD}
    ms = pl.groupby(["city", "budget"]).mode_share.mean().unstack()[BUD]
    ms["gain_low_mid"] = ms["mid"] - ms["low"]
    ms["gain_mid_high"] = ms["high"] - ms["mid"]
    S["P1f_mode_share"] = ms.round(4).reset_index().to_dict("records")

    # ---------- Q2: transfer penalty ----------
    keys = ["city", "budget", "mode", "seed"]
    q2 = {}
    for met in ["n_lines", "mean_stations_per_line", "n_transfer_stations", "transfers_per_line", "longest_line_km",
                "mean_transfers_per_trip", "rings", "fitness", "track_km", "n_stations"]:
        q2[met] = dict(abs=paired(df, "pen", "lowpen", "highpen", met, keys),
                       rel=paired(df, "pen", "lowpen", "highpen", met, keys, rel=True))
    df["km_per_line"] = df.track_km / df.n_lines
    q2["km_per_line"] = dict(abs=paired(df, "pen", "lowpen", "highpen", "km_per_line", keys),
                             rel=paired(df, "pen", "lowpen", "highpen", "km_per_line", keys, rel=True))
    S["q2_pooled"] = q2
    S["q2_seed_consistency"] = {m: seed_consistency(df, "pen", "lowpen", "highpen", m, ["city", "budget", "mode"])
                                for m in ["n_lines", "mean_stations_per_line", "n_transfer_stations", "km_per_line"]}
    mid = df[(df.budget == "mid") & (df["mode"] == "planned")]
    S["P2b_mid_lines_by_city"] = {c: dict(low=mid[(mid.city == c) & (mid.pen == "lowpen")].n_lines.mean(),
                                         high=mid[(mid.city == c) & (mid.pen == "highpen")].n_lines.mean()) for c in cities.CITY_ORDER}
    # cross-evaluation: low-penalty-optimized network scored under the high penalty vs high-penalty-optimized
    lo = df[df.pen == "lowpen"].set_index(keys)
    hi = df[df.pen == "highpen"].set_index(keys)
    j = pd.DataFrame({"hi_opt_under_hi": hi.fitness, "lo_opt_under_hi": lo.fit_cross_penalty,
                      "lo_opt_under_lo": lo.fitness, "hi_opt_under_lo": hi.fit_cross_penalty}).dropna()
    S["q2_cross"] = dict(
        ratio_hiopt_over_loopt_under_hi=float((j.hi_opt_under_hi / j.lo_opt_under_hi.clip(lower=1e-9)).median()),
        ratio_hi_fit_over_lo_fit=float((j.hi_opt_under_hi / j.lo_opt_under_lo.clip(lower=1e-9)).median()),
        frac_hi_lower=float((j.hi_opt_under_hi < j.lo_opt_under_lo).mean()),
        regret_lo_net_under_hi_rel=paired(pd.DataFrame({"f": ["a"] * len(j) + ["b"] * len(j), "v": list(j.hi_opt_under_hi) + list(j.lo_opt_under_hi),
                                                          "k": list(range(len(j))) * 2}), "f", "a", "b", "v", ["k"], rel=True))

    # ---------- Q3: phased vs planned ----------
    keys3 = ["city", "budget", "pen", "seed"]
    gap = {}
    A = df[df["mode"] == "planned"].set_index(keys3)
    B = df[df["mode"] == "phased"].set_index(keys3)
    j3 = pd.DataFrame({"pl": A.fitness, "ph": B.fitness}).dropna().reset_index()
    j3["gap"] = 1 - j3.ph / j3.pl
    cg = j3.groupby(["city", "budget", "pen"]).gap.agg(["mean", "min", "max"]).reset_index()
    S["q3_gap_by_condition"] = cg.round(4).to_dict("records")
    S["q3_gap_summary"] = dict(median_cond_gap=float(cg["mean"].median()), mean_gap=float(j3.gap.mean()),
                              worst_cond_gap=float(cg["mean"].max()), worst_run_gap=float(j3.gap.max()),
                              frac_runs_phased_worse=float((j3.gap > 0).mean()), ci=boot_ci(j3.gap.values),
                              frac_conditions_phased_worse_all_seeds=float((j3.groupby(["city", "budget", "pen"]).gap.apply(lambda s: (s > 0).all())).mean()))
    S["q3_gap_by_budget"] = j3.groupby("budget").gap.agg(["mean", "median"]).round(4).to_dict()
    S["q3_gap_by_city"] = j3.groupby("city").gap.agg(["mean", "median"]).round(4).to_dict()
    S["q3_gap_by_city_budget"] = j3.groupby(["city", "budget"]).gap.mean().round(4).unstack().to_dict()
    q3s = {}
    for met in ["n_lines", "mean_stations_per_line", "longest_share", "rings", "n_transfer_stations", "transfers_per_line", "mode_share", "track_km", "budget_used_frac"]:
        q3s[met] = dict(abs=paired(df, "mode", "planned", "phased", met, keys3), rel=paired(df, "mode", "planned", "phased", met, keys3, rel=True))
    S["q3_structure"] = q3s
    S["q3_structure_seed_consistency"] = {m: seed_consistency(df, "mode", "planned", "phased", m, ["city", "budget", "pen"])
                                          for m in ["n_lines", "mean_stations_per_line", "longest_share", "rings"]}

    # ---------- Q4/Q5 are computed in analysis_real.py-like block below ----------
    real = json.load(open(os.path.join(ROOT, "results", "real_eval.json")))
    S["real_eval"] = {c: {k: dict(fitness=v["metrics"]["fitness"], cost=v["metrics"]["cost"], mode_share=v["metrics"]["mode_share"],
                                  transfer_faster=v["metrics"]["share_flow_transit_faster"], n_lines=v["metrics"]["n_lines"],
                                  n_stations=v["metrics"]["n_stations"], track_km=v["metrics"]["track_km"]) for k, v in d.items()}
                      for c, d in real.items()}
    LAND = {
        "baltimore": {"Downtown": (17.3, 24.3), "Homewood/JHU": (17.5, 18.5), "JHH/East Balt": (20.6, 23.6), "West Balt": (10, 23),
                      "Owings Mills": (5.5, 12.8), "Hunt Valley": (20.5, 2.5), "BWI": (9.5, 36), "Towson": (20, 10)},
        "triangle": {"Raleigh DT": (38, 29), "RTP": (24.5, 22), "Durham": (14, 12), "Chapel Hill": (6, 22),
                     "Wake Forest": (43, 6), "Cary": (29, 33), "Garner": (42, 37)},
    }
    S["landmarks"] = {}
    for cn in ["baltimore", "triangle"]:
        c = cities.get_city(cn)
        sub = df[(df.city == cn) & (df["mode"] == "planned")]
        out = {}
        for _, r in sub.iterrows():
            pass
        S["landmarks"][cn] = {}
        for b in BUD:
            for p in PEN:
                rows = sub[(sub.budget == b) & (sub.pen == p)]
                served = {k: [] for k in LAND[cn]}
                for rid in rows.id:
                    L = json.load(open(os.path.join(ROOT, "results", "runs", rid + ".json")))["lines"]
                    pts = np.array([c.xy(*s) for l in L for s in l])
                    for k, (x, y) in LAND[cn].items():
                        served[k].append(float(np.hypot(pts[:, 0] - x, pts[:, 1] - y).min()))
                S["landmarks"][cn][f"{b}|{p}"] = {k: [round(v, 2) for v in vs] for k, vs in served.items()}
    # evolved vs real
    ev = {}
    for cn, rk in (("baltimore", "Metro + Light Rail + Red Line"), ("triangle", "Proposed commuter rail")):
        sub = df[(df.city == cn) & (df["mode"] == "planned")]
        for p in PEN:
            realf = real[cn][f"{rk}|{p}"]["metrics"]
            for b in BUD:
                f = sub[(sub.budget == b) & (sub.pen == p)].fitness
                ev[f"{cn}|{p}|{b}"] = dict(evolved_mean=float(f.mean()), evolved_min=float(f.min()), evolved_max=float(f.max()),
                                          real=realf["fitness"], ratio_mean=float(f.mean() / realf["fitness"]),
                                          ratio_min=float(f.min() / realf["fitness"]),
                                          mode_share=float(sub[(sub.budget == b) & (sub.pen == p)].mode_share.mean()),
                                          transit_faster=float(sub[(sub.budget == b) & (sub.pen == p)].share_flow_transit_faster.mean()),
                                          real_mode_share=realf["mode_share"])
    S["q45_evolved_vs_real"] = ev
    # per-worker savings
    S["per_worker_hours_mid"] = {cn: float(df[(df.city == cn) & (df.budget == "mid") & (df["mode"] == "planned") & (df.pen == "lowpen")].fitness.mean() / cities.get_city(cn).H.sum()) for cn in ["baltimore", "triangle"]}
    S["per_budget_unit_mid"] = {cn: float(df[(df.city == cn) & (df.budget == "mid") & (df["mode"] == "planned") & (df.pen == "lowpen")].fitness.mean() / cities.get_city(cn).mid_budget) for cn in ["baltimore", "triangle"]}
    # Triangle sensitivity of the Q5 conclusion to the driving-speed assumption (fixed evolved/proposed networks re-scored)
    sens = {}
    t = cities.make_triangle()
    rn = t.real_networks()
    prop = [l for l in rn["Proposed commuter rail (Durham-Garner-Clayton)"]]
    tsub = df[(df.city == "triangle") & (df.budget == "mid") & (df["mode"] == "planned") & (df.pen == "lowpen")]
    evolved = [json.load(open(os.path.join(ROOT, "results", "runs", rid + ".json")))["lines"] for rid in tsub.id]
    for dk in (45, 50, 55, 60):
        t.drive_kmh = dk
        t._prep_zones()
        fl = model.fleet_for(cities.get_city("triangle").mid_budget)
        fp, mp = model.evaluate(t, prop, 5., fl, detail=True)
        evs = [model.evaluate(t, [[tuple(s) for s in l] for l in L], 5., fl, detail=True) for L in evolved]
        sens[str(dk)] = dict(proposed=fp, proposed_share=mp["mode_share"], evolved=[e[0] for e in evs], evolved_share=[e[1]["mode_share"] for e in evs])
    S["triangle_driving_speed_sensitivity"] = sens

    json.dump(S, open(os.path.join(ROOT, "results", "stats.json"), "w"), indent=1, default=float)
    figures(df, real)
    return df, S


def get_lines(df, **kw):
    q = df
    for k, v in kw.items():
        q = q[q[k] == v]
    r = q.iloc[0]
    return json.loads(r.lines_json), r


def figures(df, real):
    df = load()
    # --- Fig 1: shapes x budgets
    fig, axs = plt.subplots(4, 3, figsize=(7.2, 9))
    for i, c in enumerate(SYN):
        for j, b in enumerate(BUD):
            L, r = get_lines(df, city=c, budget=b, pen="lowpen", mode="planned", seed=0)
            L = [[tuple(s) for s in l] for l in L]
            render.draw(axs[i, j], cities.get_city(c), L, title=f"{NAMES[c]} | {b} budget ({r.budget_units:.0f}) | {int(r.n_lines)} lines, "
                        f"{int(r.n_transfer_stations)} xfer, ring {int(r.rings)}", lw=1.8, ms=3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig1_shapes_by_budget.png"), dpi=130)
    plt.close(fig)

    # --- Fig 2: planned vs phased
    keys3 = ["city", "budget", "pen", "seed"]
    A = df[df["mode"] == "planned"].set_index(keys3)
    B = df[df["mode"] == "phased"].set_index(keys3)
    j3 = pd.DataFrame({"pl": A.fitness, "ph": B.fitness}).dropna().reset_index()
    j3["gap"] = 1 - j3.ph / j3.pl
    fig = plt.figure(figsize=(7.2, 6.2))
    ax = fig.add_subplot(2, 1, 1)
    order = cities.CITY_ORDER
    cols = {"low": "#9ecae1", "mid": "#3182bd", "high": "#08306b"}
    for ci, c in enumerate(order):
        for bi, b in enumerate(BUD):
            sub = j3[(j3.city == c) & (j3.budget == b)]
            x = ci + (bi - 1) * 0.25
            ax.scatter(np.full(len(sub), x) + np.random.default_rng(0).uniform(-.04, .04, len(sub)), sub.gap * 100, s=8, color=cols[b], alpha=.7,
                       label=b if ci == 0 else None)
            ax.plot([x - .09, x + .09], [sub.gap.mean() * 100] * 2, color="k", lw=1.2)
    ax.axhline(0, color="grey", lw=.6)
    ax.set_xticks(range(6)); ax.set_xticklabels([NAMES[c] for c in order])
    ax.set_ylabel("Phased shortfall vs planned (% of fitness)")
    ax.set_title("Each dot: one (penalty, seed) pair; bar: mean over those pairs. Positive = phased worse")
    ax.legend(title="budget", ncol=3, fontsize=7)
    for k, (c, b) in enumerate([("poly", "mid"), ("river", "high")]):
        for m, mode in enumerate(["planned", "phased"]):
            ax2 = fig.add_subplot(2, 4, 5 + k * 2 + m)
            L, r = get_lines(df, city=c, budget=b, pen="lowpen", mode=mode, seed=0)
            render.draw(ax2, cities.get_city(c), [[tuple(s) for s in l] for l in L],
                        title=f"{NAMES[c]} {b}, {mode}\nfit {r.fitness:.0f}, {int(r.n_lines)} lines", lw=1.6, ms=3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_planned_vs_phased.png"), dpi=130)
    plt.close(fig)

    # --- Fig 3/4: real cities
    for cn, rkey, tag in (("baltimore", "Metro + Light Rail + Red Line", "fig3_baltimore"), ("triangle", "Proposed commuter rail", "fig4_triangle")):
        c = cities.get_city(cn)
        fig, axs = plt.subplots(1, 3, figsize=(7.4, 3.2), gridspec_kw={"width_ratios": [1, 1, 1.1]})
        rl = real[cn][f"{rkey}|lowpen"]
        render.draw(axs[0], c, [[tuple(s) for s in l] for l in rl["lines"]],
                    title=f"{'Real+planned' if cn=='baltimore' else 'Proposed'} (stylized)\nfit {rl['metrics']['fitness']:.0f}, cost {rl['metrics']['cost']:.0f}", lw=1.8)
        L, r = get_lines(df, city=cn, budget="mid", pen="lowpen", mode="planned", seed=0)
        render.draw(axs[1], c, [[tuple(s) for s in l] for l in L], title=f"Evolved, mid budget, seed 0\nfit {r.fitness:.0f}, cost {r.cost:.0f}", lw=1.8)
        ax = axs[2]
        labels, vals = [], []
        for k, v in real[cn].items():
            if v["pen"] == "lowpen":
                labels.append(v["label"].replace(" + ", "+").replace("Light Rail", "LR").replace("Line", "Ln")); vals.append(v["metrics"]["fitness"])
        ax.barh(range(len(vals)), vals, color="#bbbbbb")
        base = len(vals)
        for k, b in enumerate(BUD):
            f = df[(df.city == cn) & (df.budget == b) & (df.pen == "lowpen") & (df["mode"] == "planned")].fitness
            ax.barh(base + k, f.mean(), color="#d62728" if b == "mid" else "#f4a6a6")
            ax.scatter(f.values, np.full(len(f), base + k), color="k", s=7, zorder=3)
            labels.append(f"Evolved {b}")
        ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=6.5)
        ax.set_xlabel("Fitness (person-h/day saved)"); ax.invert_yaxis()
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, tag + ".png"), dpi=130)
        plt.close(fig)

    # --- Fig 5: topology metrics vs budget (planned)
    pl = df[df["mode"] == "planned"]
    mets = [("n_lines", "Lines"), ("mean_stations_per_line", "Stations per line"), ("n_transfer_stations", "Transfer stations"),
            ("rings", "Independent cycles (rings)"), ("mode_share", "Transit mode share"), ("longest_share", "Longest line's share of track")]
    fig, axs = plt.subplots(2, 3, figsize=(7.4, 5))
    cc = dict(zip(cities.CITY_ORDER, ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"]))
    for ax, (m, t) in zip(axs.ravel(), mets):
        for c in cities.CITY_ORDER:
            for p, ls in (("lowpen", "-"), ("highpen", "--")):
                y = [pl[(pl.city == c) & (pl.budget == b) & (pl.pen == p)][m].mean() for b in BUD]
                ax.plot(range(3), y, ls, color=cc[c], marker="o", ms=3, lw=1.2, label=f"{NAMES[c]}" if p == "lowpen" else None)
        ax.set_xticks(range(3)); ax.set_xticklabels(["low", "mid", "high"]); ax.set_title(t)
    axs[0, 0].legend(fontsize=6)
    fig.suptitle("Planned networks, mean of 3 seeds. Solid: low transfer penalty; dashed: high", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig5_topology_vs_budget.png"), dpi=130)
    plt.close(fig)

    # --- Fig 6: transfer penalty effect (paired differences)
    keys = ["city", "budget", "mode", "seed"]
    fig, axs = plt.subplots(1, 4, figsize=(7.4, 2.6))
    for ax, (m, t) in zip(axs, [("n_lines", "Lines"), ("mean_stations_per_line", "Stations/line"),
                                ("n_transfer_stations", "Transfer stations"), ("mean_transfers_per_trip", "Transfers per rider trip")]):
        A = df[df.pen == "lowpen"].set_index(keys)[m]; B = df[df.pen == "highpen"].set_index(keys)[m]
        d = (B - A).dropna().reset_index()
        for ci, c in enumerate(cities.CITY_ORDER):
            v = d[d.city == c][m]
            ax.scatter(np.full(len(v), ci) + np.random.default_rng(1).uniform(-.15, .15, len(v)), v, s=5, color=cc[c], alpha=.6)
            ax.plot([ci - .25, ci + .25], [v.mean()] * 2, color="k", lw=1)
        ax.axhline(0, color="grey", lw=.6)
        ax.set_xticks(range(6)); ax.set_xticklabels(["Mono", "Coast", "River", "Poly", "Balt", "Tri"], fontsize=6)
        ax.set_title(t + " (high - low)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig6_transfer_penalty.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
