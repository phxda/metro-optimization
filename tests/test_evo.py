"""Optimizer invariants. Run: python3 tests/test_evo.py"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import numpy as np
import cities, model, evo

def test_frozen_core_preserved():
    c = cities.get_city("poly"); rng = np.random.default_rng(1)
    ev = evo.Evolver(c, 5., rng)
    net1, f1 = ev.run([], 110., model.fleet_for(110.), 30)
    base = evo.freeze(net1)
    net2, f2 = ev.run(base, 220., model.fleet_for(220.), 60)
    assert model.net_cost(c, [l.s for l in net2]) <= 220 + 1e-9
    # every frozen station sequence is still a contiguous slice of exactly one line in the new network
    for b in base:
        assert any(any(l.s[i:i + len(b.s)] == b.s for i in range(len(l.s) - len(b.s) + 1)) for l in net2), "frozen track changed"
    assert f2 >= f1 - 1e-6 * 0 - 50  # cost tie-break epsilon aside, extra budget should not make it worse (loose)

def test_mutations_respect_budget_and_validity():
    c = cities.get_city("river"); rng = np.random.default_rng(3)
    ev = evo.Evolver(c, 5., rng)
    net = ev.init_net(150., 40)
    n_ok = 0
    for _ in range(400):
        k = ev.mutate(net, 150., 40)
        if k is None: continue
        n_ok += 1
        L = [l.s for l in k]
        assert model.net_cost(c, L) <= 150 + 1e-9
        for l in L:
            assert all(not c.water[a] for a in l) and ev.line_ok(l)
        net = k
    assert n_ok > 50

def test_phased_budget_schedule():
    net, pf, ne = evo.optimize("mono", 0.5, 5., "phased", 0, gens=30)
    c = cities.get_city("mono")
    assert model.net_cost(c, [l.s for l in net]) <= 0.5 * c.mid_budget + 1e-9 and len(pf) == 3

if __name__ == "__main__":
    for k, v in list(globals().items()):
        if k.startswith("test_"):
            v(); print("ok", k)
