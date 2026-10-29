"""v2 robustness suites (resumable). Suites: seeds (seeds 3-9, base), sens (8 perturbed variants, seeds 0-1), long (1500 gens, planned, lowpen, seeds 0-2).
Usage: python3 run_sweep2.py --suite sens [--procs 4]"""
import sys, os, json, time, itertools, argparse
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "src"))
from multiprocessing import Pool
import cities, model, evo

OUT = os.path.join(ROOT, "results", "v2")
BUDGETS = {"low": 0.5, "mid": 1.0, "high": 2.0}
PENS = {"lowpen": 5.0, "highpen": 20.0}


def tasks(suite):
    C, B, P = cities.CITY_ORDER, list(BUDGETS), list(PENS)
    if suite == "seeds":
        it = itertools.product(["base"], C, B, P, ["planned", "phased"], range(3, 10), [500])
    elif suite == "sens":
        it = itertools.product([v for v in cities.VARIANTS if v != "base"], C, B, P, ["planned", "phased"], range(2), [500])
    elif suite == "long":
        it = itertools.product(["base"], C, B, ["lowpen"], ["planned"], range(3), [1500])
    return [dict(variant=v, city=c, budget=b, pen=p, mode=m, seed=s, gens=g) for v, c, b, p, m, s, g in it]


def tid(t):
    return f"{t['variant']}_{t['city']}_{t['budget']}_{t['pen']}_{t['mode']}_s{t['seed']}_g{t['gens']}"


def run_one(t):
    d = os.path.join(OUT, t["suite"])
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, tid(t) + ".json")
    if os.path.exists(path):
        return "skip"
    t0 = time.time()
    city = cities.get_variant_city(t["city"], t["variant"])
    pen = PENS[t["pen"]]
    net, phase_fit, evals = evo.optimize(t["city"], BUDGETS[t["budget"]], pen, t["mode"], t["seed"], t["gens"], t["variant"])
    lines = [l.s for l in net]
    budget = cities.get_city(t["city"]).mid_budget * BUDGETS[t["budget"]]
    fleet = model.fleet_for(budget, city.fleet_per_unit)
    f, m = model.evaluate(city, lines, pen, fleet, detail=True)
    rec = dict(t, id=tid(t), budget_units=budget, fleet=fleet, penalty=pen, lines=lines, metrics=m,
               phase_fit=phase_fit, evals=evals, secs=round(time.time() - t0, 1))
    json.dump(rec, open(path + ".tmp", "w"))
    os.replace(path + ".tmp", path)
    return "done"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", required=True)
    ap.add_argument("--procs", type=int, default=4)
    a = ap.parse_args()
    ts = [dict(t, suite=a.suite) for t in tasks(a.suite)]
    t0 = time.time()
    with Pool(a.procs) as p:
        res = list(p.imap_unordered(run_one, ts, chunksize=1))
    print(f"suite {a.suite}: {len(ts)} tasks, {res.count('done')} done, {res.count('skip')} skipped, {time.time()-t0:.0f}s", flush=True)
