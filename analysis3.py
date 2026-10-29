"""Real-data (LODES + TIGER) analysis: demand-model validation, evolved vs real networks, structure, figures."""
import sys, os, json, glob, copy
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cities, model, render
from analysis import boot_ci

BUD, PEN = ["low", "mid", "high"], ["lowpen", "highpen"]
PV = {"lowpen": 5.0, "highpen": 20.0}
FIG = os.path.join(ROOT, "figures")
REAL = {"baltimore_lodes": ["Metro SubwayLink", "Light RailLink", "Red Line (planned)"],
        "triangle_lodes": ["Proposed commuter rail (Durham-Garner-Clayton)"]}


def load():
    rows = []
    for p in glob.glob(os.path.join(ROOT, "results", "v2", "real", "*.json")):
        r = json.load(open(p))
        rows.append(dict({k: r[k] for k in ("city", "budget", "pen", "mode", "seed", "budget_units")}, lines=r["lines"], phase_fit=r["phase_fit"], **r["metrics"]))
    return pd.DataFrame(rows)


def gravity_twin(c):
    """Same grid, homes and jobs, but flows from the gravity model with gamma fitted to the observed mean trip length."""
    d = np.hypot(c.zx[:, None] - c.zx[None], c.zy[:, None] - c.zy[None])
    obs = c.flow
    w = obs / obs.sum()
    mean_obs = float((w * d).sum())
    best = None
    for g in np.arange(2, 40, 1.0):
        t = copy.deepcopy(c)
        t.od = None; t.gamma = g
        dd = d.copy(); np.fill_diagonal(dd, 0.5 * c.block * c.cell)
        f = np.exp(-dd / g) * c.zJ[None]; np.fill_diagonal(f, 0)
        fl = c.zH[:, None] * f / f.sum(1, keepdims=True).clip(1e-12)
        fl = fl * (obs.sum() / fl.sum())
        m = float((fl / fl.sum() * d).sum())
        if best is None or abs(m - mean_obs) < best[0]:
            best = (abs(m - mean_obs), g, fl, m)
    _, g, fl, m = best
    t = copy.deepcopy(c); t.flow = fl; t.total_trips = fl.sum()
    return t, g, m, mean_obs


def main():
    S = {}
    df = load()
    S["n"] = len(df)
    pl = df[df["mode"] == "planned"]
    # --- validation of demand model
    S["validation"] = {}
    twins = {}
    for n in REAL:
        c = cities.get_city(n)
        t, g, m, mo = gravity_twin(c)
        twins[n] = t
        a, b = c.flow.ravel(), t.flow.ravel()
        mask = (a > 0) | (b > 0)
        d = np.hypot(c.zx[:, None] - c.zx[None], c.zy[:, None] - c.zy[None]).ravel()
        S["validation"][n] = dict(gamma_fit=float(g), mean_trip_km_obs=mo, mean_trip_km_gravity=m,
                                  corr_flows=float(np.corrcoef(a, b)[0, 1]),
                                  corr_log=float(np.corrcoef(np.log1p(a), np.log1p(b))[0, 1]),
                                  share_le5km_obs=float(a[d <= 5].sum() / a.sum()), share_le5km_gravity=float(b[d <= 5].sum() / b.sum()),
                                  share_ge20km_obs=float(a[d >= 20].sum() / a.sum()), share_ge20km_gravity=float(b[d >= 20].sum() / b.sum()),
                                  zones=int(c.Z), trips_per_day=float(c.total_trips))
        # real network under observed vs gravity demand
        rn = c.real_networks(); L = [l for k in REAL[n] for l in rn[k]]
        fl = model.fleet_for(c.mid_budget)
        S["validation"][n]["real_fit_observed"] = model.evaluate(c, L, 5., fl)
        S["validation"][n]["real_fit_gravity"] = model.evaluate(t, L, 5., fl)
        fo = model.evaluate(c, L, 5., fl, detail=True)[1]
        S["validation"][n]["real_mode_share_observed"] = fo["mode_share"]
    # --- evolved vs real
    S["vs_real"] = {}
    for n in REAL:
        c = cities.get_city(n); rn = c.real_networks()
        for p in PEN:
            for part, ks in (("all", REAL[n]),):
                L = [l for k in ks for l in rn[k]]
                rf = model.evaluate(c, L, PV[p], model.fleet_for(c.mid_budget))
                for b in BUD:
                    f = pl[(pl.city == n) & (pl.budget == b) & (pl.pen == p)].fitness
                    S["vs_real"][f"{n}|{p}|{b}"] = dict(real=rf, evolved_mean=float(f.mean()), evolved_min=float(f.min()), evolved_max=float(f.max()),
                                                       ratio=float(f.mean() / rf), ratio_min=float(f.min() / rf), n=int(len(f)),
                                                       mode_share=float(pl[(pl.city == n) & (pl.budget == b) & (pl.pen == p)].mode_share.mean()))
        # components
        S.setdefault("real_components", {})[n] = {k: model.evaluate(c, rn[k], 5., model.fleet_for(c.mid_budget)) for k in rn}
        if n == "baltimore_lodes":
            ml = [l for k in ("Metro SubwayLink", "Light RailLink") for l in rn[k]]
            S["real_components"][n]["Metro+LRT"] = model.evaluate(c, ml, 5., model.fleet_for(c.mid_budget))
    # --- evolved networks rescored under gravity demand (does demand realism change what is good?)
    S["rescored_under_gravity"] = {}
    for n in REAL:
        sub = pl[(pl.city == n) & (pl.budget == "mid") & (pl.pen == "lowpen")]
        fl = model.fleet_for(cities.get_city(n).mid_budget)
        S["rescored_under_gravity"][n] = [dict(observed=float(r.fitness), gravity=model.evaluate(twins[n], [[tuple(s) for s in l] for l in r.lines], 5., fl)) for _, r in sub.iterrows()]
    # --- landmarks served
    S["landmarks"] = {}
    for n in REAL:
        c = cities.get_city(n)
        for b in BUD:
            sub = pl[(pl.city == n) & (pl.budget == b) & (pl.pen == "lowpen")]
            res = {k: [] for k in c.landmarks}
            for _, r in sub.iterrows():
                pts = np.array([c.xy(*s) for l in r.lines for s in l])
                for k, (x, y) in c.landmarks.items():
                    res[k].append(round(float(np.hypot(pts[:, 0] - x, pts[:, 1] - y).min()), 1))
            S["landmarks"][f"{n}|{b}"] = res
    # --- effects in real-data cities (5 seeds)
    k = ["city", "budget", "mode", "seed"]
    def pr(a, b, factor, metric, keys, rel=False):
        A = df[df[factor] == a].set_index(keys)[metric]; B = df[df[factor] == b].set_index(keys)[metric]
        j = pd.concat([A.rename("a"), B.rename("b")], axis=1).dropna()
        d = (j.b - j.a) / j.a.abs().clip(lower=1e-9) if rel else (j.b - j.a)
        lo, hi = boot_ci(d.values)
        return dict(n=len(d), mean=float(d.mean()), ci=(lo, hi))
    S["pen_effect"] = {m: pr("lowpen", "highpen", "pen", m, k) for m in ["n_lines", "mean_stations_per_line", "n_transfer_stations", "rings"]}
    S["pen_effect"]["fitness_rel"] = pr("lowpen", "highpen", "pen", "fitness", k, rel=True)
    k3 = ["city", "budget", "pen", "seed"]
    A = df[df["mode"] == "planned"].set_index(k3).fitness; B = df[df["mode"] == "phased"].set_index(k3).fitness
    g = (1 - B / A).dropna().reset_index(name="gap")
    S["phased_gap"] = dict(mean=float(g.gap.mean()), median=float(g.gap.median()), frac_worse=float((g.gap > 0).mean()), ci=boot_ci(g.gap.values),
                           by_city=g.groupby("city").gap.mean().round(4).to_dict(), by_budget=g.groupby("budget").gap.mean().round(4).to_dict())
    S["table"] = df[df["mode"] == "planned"].groupby(["city", "budget"])[["n_lines", "n_transfer_stations", "rings", "mode_share", "fitness", "budget_used_frac" if "budget_used_frac" in df else "cost"]].mean().round(3).reset_index().to_dict("records")
    json.dump(S, open(os.path.join(ROOT, "results", "stats_real.json"), "w"), indent=1, default=float)
    # --- figure
    fig, axs = plt.subplots(2, 4, figsize=(9.5, 5.2))
    for i, n in enumerate(REAL):
        c = cities.get_city(n); rn = c.real_networks()
        render.draw(axs[i, 0], c, [l for k2 in REAL[n] for l in rn[k2]], title=f"{n.split('_')[0].title()} real/proposed (approx.)\nfit {S['vs_real'][n+'|lowpen|mid']['real']:.0f}", lw=1.6)
        for j, b in enumerate(["low", "mid", "high"]):
            r = pl[(pl.city == n) & (pl.budget == b) & (pl.pen == "lowpen") & (pl.seed == 0)].iloc[0]
            render.draw(axs[i, 1 + j], c, [[tuple(s) for s in l] for l in r.lines], title=f"Evolved {b}, seed 0\nfit {r.fitness:.0f}, {int(r.n_lines)} lines", lw=1.6)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig7_real_data.png"), dpi=130); plt.close(fig)
    return S


if __name__ == "__main__":
    S = main()
    print(json.dumps({k: S[k] for k in ("n", "validation", "vs_real")}, indent=1, default=float)[:5000])
