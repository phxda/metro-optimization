"""Sweep driver. Resumable: one JSON per run in results/runs/. Usage: python3 run_sweep.py [--chunk i --of n] [--gens G]"""
import sys, os, json, time, itertools, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
from multiprocessing import Pool
import cities, model, evo

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "runs")
BUDGETS = {"low": 0.5, "mid": 1.0, "high": 2.0}
PENS = {"lowpen": 5.0, "highpen": 20.0}
MODES = ["planned", "phased"]
SEEDS = [0, 1, 2]


def all_tasks():
    return [dict(city=c, budget=b, pen=p, mode=m, seed=s)
            for c, b, p, m, s in itertools.product(cities.CITY_ORDER, BUDGETS, PENS, MODES, SEEDS)]


def tid(t):
    return f"{t['city']}_{t['budget']}_{t['pen']}_{t['mode']}_s{t['seed']}"


def run_one(args):
    t, gens = args
    path = os.path.join(OUT, tid(t) + ".json")
    if os.path.exists(path):
        return tid(t), "skip"
    t0 = time.time()
    city = cities.get_city(t["city"])
    pen, other = PENS[t["pen"]], PENS["highpen" if t["pen"] == "lowpen" else "lowpen"]
    net, phase_fit, evals = evo.optimize(t["city"], BUDGETS[t["budget"]], pen, t["mode"], t["seed"], gens)
    lines = [l.s for l in net]
    budget = city.mid_budget * BUDGETS[t["budget"]]
    fleet = model.fleet_for(budget)
    f, m = model.evaluate(city, lines, pen, fleet, detail=True)
    rec = dict(t, id=tid(t), budget_units=budget, fleet=fleet, penalty=pen, lines=lines, metrics=m,
               phase_fit=phase_fit, evals=evals, fit_cross_penalty=model.evaluate(city, lines, other, fleet),
               secs=round(time.time() - t0, 1), gens=gens)
    with open(path + ".tmp", "w") as fh:
        json.dump(rec, fh)
    os.replace(path + ".tmp", path)
    return tid(t), "done"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--gens", type=int, default=500)
    ap.add_argument("--procs", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    tasks = all_tasks()[a.chunk::a.of]
    t0 = time.time()
    with Pool(a.procs) as p:
        res = list(p.imap_unordered(run_one, [(t, a.gens) for t in tasks]))
    print(f"chunk {a.chunk}/{a.of}: {sum(r[1]=='done' for r in res)} done, {sum(r[1]=='skip' for r in res)} skipped, {time.time()-t0:.0f}s")
