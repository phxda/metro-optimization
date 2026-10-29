"""Evolutionary optimizer: (mu+lambda) over networks (lists of ordered station lists)."""
import math
import numpy as np
import model
from cities import get_city

MU, LAM = 16, 32


class Line:
    __slots__ = ("s", "lo", "hi")

    def __init__(self, s, lo=0, hi=0):
        self.s, self.lo, self.hi = list(s), lo, hi

    def copy(self):
        return Line(self.s, self.lo, self.hi)

    def movable(self, i):
        return i < self.lo or i >= self.hi


def key(net):
    return tuple((tuple(l.s), l.lo, l.hi) for l in net)


class Evolver:
    def __init__(self, city, pen, rng):
        self.c, self.pen, self.rng = city, pen, rng
        self.cache = {}
        self.evals = 0
        w = city.demand_w.ravel()
        self.cells_w = w / w.sum()
        self.cells_idx = np.arange(city.R * city.C)

    # -------- helpers --------
    def hop_ok(self, a, b):
        km = self.c.seg(a, b)[0]
        return self.c.min_hop - 1e-9 <= km <= self.c.max_hop + 1e-9

    def line_ok(self, s):
        if len(s) < 2 or len(set(s)) != len(s):
            return False
        return all(self.hop_ok(a, b) for a, b in zip(s[:-1], s[1:]))

    def random_cell(self):
        i = self.rng.choice(self.cells_idx, p=self.cells_w)
        return divmod(int(i), self.c.C)

    def near_cell(self, p, rmin, rmax):
        c = self.c
        r0, c0 = p
        R = int(math.ceil(rmax))
        rs = np.arange(max(0, r0 - R), min(c.R, r0 + R + 1))
        cs = np.arange(max(0, c0 - R), min(c.C, c0 + R + 1))
        rr, cc = np.meshgrid(rs, cs, indexing="ij")
        d = np.hypot(rr - r0, cc - c0)
        ok = (d >= rmin) & (d <= rmax) & c.land[rr, cc]
        if not ok.any():
            return None
        w = (c.demand_w[rr, cc] + 1.0)[ok] ** 0.7
        k = self.rng.choice(len(w), p=w / w.sum())
        return int(rr[ok][k]), int(cc[ok][k])

    def maybe_snap(self, p, net, skip):
        """With some probability, snap p to a nearby station of another line (creates a transfer)."""
        if self.rng.random() > 0.3:
            return p
        best, bd = None, 3.01
        for j, l in enumerate(net):
            if j == skip:
                continue
            for q in l.s:
                d = math.hypot(q[0] - p[0], q[1] - p[1])
                if 0 < d < bd:
                    best, bd = q, d
        return best if best is not None else p

    def cap_cells(self):
        """hop length in cells to sample for new stations"""
        cell = self.c.cell
        return max(1.0, 0.9 / cell), self.c.max_hop / cell

    # -------- mutations --------
    def m_move(self, net):
        cand = [(j, i) for j, l in enumerate(net) for i in range(len(l.s)) if l.movable(i)]
        if not cand:
            return None
        j, i = cand[self.rng.integers(len(cand))]
        l = net[j]
        rmin, rmax = 1.0, min(3.0, self.cap_cells()[1])
        p = self.near_cell(l.s[i], rmin, rmax)
        if p is None:
            return None
        p = self.maybe_snap(p, net, j)
        s = l.s[:i] + [p] + l.s[i + 1:]
        if not self.line_ok(s):
            return None
        l.s = s
        return net

    def m_extend(self, net):
        movable = [j for j, l in enumerate(net)]
        if not movable:
            return None
        j = movable[self.rng.integers(len(movable))]
        l = net[j]
        mn, mx = self.cap_cells()
        r = self.rng.random()
        if r < 0.7 or len(l.s) < 2:       # extend at an end
            end = self.rng.integers(2)
            ref = l.s[-1] if end else l.s[0]
            prev = l.s[-2] if end else l.s[1]
            p = self.near_cell(ref, mn, min(mx, 5.0 / self.c.cell + 2))
            if p is None:
                return None
            if self.rng.random() < 0.5:   # prefer continuing direction
                v = np.array(ref) - np.array(prev)
                q = np.array(p) - np.array(ref)
                if v @ q < 0:
                    p2 = (int(ref[0] + (ref[0] - p[0])), int(ref[1] + (ref[1] - p[1])))
                    if 0 <= p2[0] < self.c.R and 0 <= p2[1] < self.c.C and self.c.land[p2]:
                        p = p2
            p = self.maybe_snap(p, net, j)
            s = (l.s + [p]) if end else ([p] + l.s)
            if not self.line_ok(s):
                return None
            if not end and l.hi > l.lo:
                l.lo += 1; l.hi += 1
            l.s = s
            return net
        # insert an interior station
        pos = [i for i in range(1, len(l.s)) if i <= l.lo or i >= l.hi]
        if not pos:
            return None
        i = pos[self.rng.integers(len(pos))]
        a, b = l.s[i - 1], l.s[i]
        mid = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
        p = self.near_cell((int(round(mid[0])), int(round(mid[1]))), 0.0, 1.5)
        if p is None:
            return None
        p = self.maybe_snap(p, net, j)
        s = l.s[:i] + [p] + l.s[i:]
        if not self.line_ok(s):
            return None
        if i <= l.lo and l.hi > l.lo:
            l.lo += 1; l.hi += 1
        l.s = s
        return net

    def m_trim(self, net):
        cand = [(j, i) for j, l in enumerate(net) for i in range(len(l.s)) if l.movable(i)]
        if not cand:
            return None
        j, i = cand[self.rng.integers(len(cand))]
        l = net[j]
        if l.hi > l.lo and i < l.lo and i != 0:
            return None        # only trim ahead of a frozen core from its outer end
        s = l.s[:i] + l.s[i + 1:]
        if len(s) < 2:
            if l.hi > l.lo:
                return None
            del net[j]
            return net
        if not self.line_ok(s):
            return None
        if i < l.lo:
            l.lo -= 1; l.hi -= 1
        l.s = s
        return net

    def new_line(self, net):
        if len(net) >= self.c.max_lines:
            return None
        st = []
        if net and self.rng.random() < 0.5:
            l = net[self.rng.integers(len(net))]
            p0 = l.s[self.rng.integers(len(l.s))]
        else:
            p0 = self.random_cell()
            if self.c.water[p0]:
                return None
        st = [p0]
        mn, mx = self.cap_cells()
        n = int(self.rng.integers(2, 5))
        ang = self.rng.random() * 2 * math.pi
        for _ in range(n - 1):
            step = self.rng.uniform(mn, min(mx, 3.5 / self.c.cell + 1))
            p = (int(round(st[-1][0] + step * math.sin(ang))), int(round(st[-1][1] + step * math.cos(ang))))
            ang += self.rng.normal(0, 0.4)
            if not (0 <= p[0] < self.c.R and 0 <= p[1] < self.c.C) or self.c.water[p]:
                break
            st.append(p)
        if not self.line_ok(st):
            return None
        net.append(Line(st))
        return net

    def m_delete(self, net):
        c = [j for j, l in enumerate(net) if l.hi == l.lo]
        if not c:
            return None
        del net[c[self.rng.integers(len(c))]]
        return net

    def m_split(self, net):
        c = [j for j, l in enumerate(net) if l.hi == l.lo and len(l.s) >= 4 and len(net) < self.c.max_lines]
        if not c:
            return None
        j = c[self.rng.integers(len(c))]
        l = net[j]
        i = int(self.rng.integers(1, len(l.s) - 1))
        share = self.rng.random() < 0.5
        a = l.s[:i + 1] if share else l.s[:i]
        b = l.s[i:] if share else l.s[i:]
        if len(a) < 2 or len(b) < 2:
            return None
        net[j] = Line(a)
        net.append(Line(b))
        return net

    def m_merge(self, net):
        c = [j for j, l in enumerate(net) if l.hi == l.lo]
        if len(c) < 2:
            return None
        j, k = self.rng.choice(c, 2, replace=False)
        A, B = net[j].s, net[k].s
        best = None
        for ra in (False, True):
            for rb in (False, True):
                a = A[::-1] if ra else A
                b = B[::-1] if rb else B
                d = math.hypot(a[-1][0] - b[0][0], a[-1][1] - b[0][1])
                if best is None or d < best[0]:
                    best = (d, a, b)
        d, a, b = best
        s = a + (b[1:] if a[-1] == b[0] else b)
        if not self.line_ok(s):
            return None
        net[j] = Line(s)
        del net[k]
        return net

    OPS = [("m_move", .30), ("m_extend", .24), ("m_trim", .10), ("new_line", .10), ("m_delete", .06),
           ("m_split", .10), ("m_merge", .10)]

    def mutate(self, parent, cap, fleet):
        names = [o[0] for o in self.OPS]
        pr = np.array([o[1] for o in self.OPS])
        pr = pr / pr.sum()
        for _ in range(12):
            net = [l.copy() for l in parent]
            k = 1 + int(self.rng.geometric(0.55) - 1)
            ok = True
            for _ in range(min(k, 3)):
                op = getattr(self, names[self.rng.choice(len(names), p=pr)])
                res = op(net)
                if res is None:
                    ok = False
                    break
                net = res
            if not ok or not net:
                continue
            if model.net_cost(self.c, [l.s for l in net]) > cap + 1e-9:
                continue
            if len(net) * model.V_MIN > fleet:
                continue
            return net
        return None

    def fit(self, net, fleet):
        k = (key(net), fleet)
        v = self.cache.get(k)
        if v is None:
            v = model.evaluate(self.c, [l.s for l in net], self.pen, fleet)
            v -= 1e-3 * model.net_cost(self.c, [l.s for l in net])   # tie-break toward cheaper networks
            self.evals += 1
            if len(self.cache) > 20000:
                self.cache.clear()
            self.cache[k] = v
        return v

    def init_net(self, cap, fleet):
        for _ in range(50):
            net = []
            for _ in range(int(self.rng.integers(1, 3))):
                r = self.new_line(net)
                if r is not None:
                    net = r
            if net and model.net_cost(self.c, [l.s for l in net]) <= cap and len(net) * model.V_MIN <= fleet:
                return net
        return []

    def run(self, base, cap, fleet, gens):
        """Evolve starting from `base` (frozen cores respected). Returns best net."""
        pop = []
        if base:
            pop.append([l.copy() for l in base])
        while len(pop) < MU:
            if base:
                n = self.mutate(base, cap, fleet)
                pop.append(n if n is not None else [l.copy() for l in base])
            else:
                pop.append(self.init_net(cap, fleet))
        sc = [self.fit(n, fleet) for n in pop]
        for _ in range(gens):
            kids = []
            for _ in range(LAM):
                par = pop[self.rng.integers(len(pop))]
                if self.rng.random() < 0.5:   # tournament bias toward good parents
                    par2 = pop[self.rng.integers(len(pop))]
                    if self.fit(par2, fleet) > self.fit(par, fleet):
                        par = par2
                k = self.mutate(par, cap, fleet)
                if k is not None:
                    kids.append(k)
            allp = pop + kids
            allsc = sc + [self.fit(n, fleet) for n in kids]
            order = np.argsort(allsc)[::-1]
            pop, sc, seen = [], [], set()
            for o in order:
                kk = key(allp[o])
                if kk in seen:
                    continue
                seen.add(kk)
                pop.append(allp[o]); sc.append(allsc[o])
                if len(pop) == MU:
                    break
        return pop[0], sc[0]


def freeze(net):
    return [Line(l.s, 0, len(l.s)) for l in net]


def optimize(city_name, budget_factor, pen, mode, seed, gens=200):
    city = get_city(city_name)
    budget = city.mid_budget * budget_factor
    rng = np.random.default_rng(seed * 1009 + 7)
    ev = Evolver(city, pen, rng)
    full_fleet = model.fleet_for(budget)
    if mode == "planned":
        net, f = ev.run([], budget, full_fleet, gens)
        phase_fit = [f]
    else:
        net, phase_fit = [], []
        for p in (1, 2, 3):
            cap = budget * p / 3
            fl = model.fleet_for(cap)
            base = freeze(net) if net else []
            net, f = ev.run(base, cap, fl, max(1, gens // 3))
            phase_fit.append(f)
    return net, phase_fit, ev.evals
