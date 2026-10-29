"""Network cost and fitness (rider time saved vs driving)."""
import math
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from cities import (COST_KM, COST_STATION, COST_WATER_KM, LAYOVER_MIN, V_MIN, HEADWAY_MIN,
                    HEADWAY_MAX, BETA, ACCESS_K)

INF = float("inf")


def net_cost(city, lines):
    """lines: list of station lists [(r,c),...]. Unique stations charged once; track charged per line-hop."""
    km = water = 0.0
    st = set()
    for ln in lines:
        st.update(ln)
        for a, b in zip(ln[:-1], ln[1:]):
            k, w = city.seg(a, b)
            km += k
            water += w
    return COST_KM * km + COST_STATION * len(st) + COST_WATER_KM * water


def fleet_for(budget):
    return int(round(0.2 * budget))


def line_times(city, ln):
    """per-hop in-vehicle minutes (incl. dwell at arrival station)."""
    return [city.seg(a, b)[0] / city.transit_kmh * 60 + city.dwell for a, b in zip(ln[:-1], ln[1:])]


def headways(city, lines, fleet):
    """Equal-headway allocation with a V_MIN floor. Returns headways (min) or None if infeasible."""
    n = len(lines)
    if n == 0:
        return []
    if fleet < V_MIN * n:
        return None
    cyc = np.array([2 * sum(line_times(city, ln)) + 2 * LAYOVER_MIN for ln in lines])
    veh = np.zeros(n)
    free = np.ones(n, bool)
    left = float(fleet)
    for _ in range(n + 1):
        share = left * cyc * free / (cyc * free).sum()
        low = free & (share < V_MIN)
        if not low.any():
            veh[free] = share[free]
            break
        veh[low] = V_MIN
        left -= V_MIN * low.sum()
        free &= ~low
        if not free.any():
            break
    hw = cyc / np.maximum(veh, 1e-9)
    return list(np.clip(hw, HEADWAY_MIN, HEADWAY_MAX))


def station_graph(lines):
    idx = {}
    for ln in lines:
        for c in ln:
            if c not in idx:
                idx[c] = len(idx)
    return idx


def _cyclomatic(lines, idx):
    edges = set()
    for ln in lines:
        for a, b in zip(ln[:-1], ln[1:]):
            edges.add((min(idx[a], idx[b]), max(idx[a], idx[b])))
    par = list(range(len(idx)))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for a, b in edges:
        par[f(a)] = f(b)
    comps = len({f(i) for i in range(len(idx))})
    return len(edges) - len(idx) + comps, comps


def station_times(city, lines, pen, fleet, want_pred=False):
    """P x P generalized minutes between physical stations (wait*w + ride + pen*transfers)."""
    idx = station_graph(lines)
    P = len(idx)
    hw = headways(city, lines, fleet)
    if hw is None:
        return None
    src, dst, wt = [], [], []
    nid = P
    for ln, h in zip(lines, hw):
        base = nid
        times = line_times(city, ln)
        for i, c in enumerate(ln):
            s = idx[c]
            src.append(s); dst.append(base + i); wt.append(city.wait_w * h / 2 + pen / 2 + 1e-6)
            src.append(base + i); dst.append(s); wt.append(pen / 2 + 1e-6)
        for i, t in enumerate(times):
            src += [base + i, base + i + 1]; dst += [base + i + 1, base + i]; wt += [t, t]
        nid += len(ln)
    G = coo_matrix((wt, (src, dst)), shape=(nid, nid)).tocsr()
    if want_pred:
        D, pred = dijkstra(G, indices=np.arange(P), return_predecessors=True)
    else:
        D, pred = dijkstra(G, indices=np.arange(P)), None
    T = D[:, :P] - pen
    np.fill_diagonal(T, INF)
    return idx, hw, T, pred, nid


def transfer_counts(P, pred, nid):
    """X[i,j] = number of transfers on the shortest path i->j (station nodes visited minus 2)."""
    X = np.zeros((P, P))
    for i in range(P):
        cnt = np.full(nid, -1, int)
        cnt[i] = 0
        p = pred[i]
        for j in range(nid):
            if cnt[j] >= 0 or p[j] < 0:
                continue
            stack = []
            k = j
            while cnt[k] < 0 and p[k] >= 0:
                stack.append(k)
                k = p[k]
            base = cnt[k]
            for k2 in reversed(stack):
                base += 1 if k2 < P else 0
                cnt[k2] = base
        X[i] = np.maximum(cnt[:P] - 1, 0)
    return X


def _access(city, lines, idx):
    pts = np.array([city.xy(*c) for c in idx])
    d = np.hypot(city.zx[:, None] - pts[None, :, 0], city.zy[:, None] - pts[None, :, 1])
    off = 0.3 * city.block * city.cell
    d = np.sqrt(d * d + off * off)
    P = d.shape[1]
    K = min(ACCESS_K, P)
    S = np.argpartition(d, K - 1, axis=1)[:, :K] if K < P else np.tile(np.arange(P), (d.shape[0], 1))
    dk = np.take_along_axis(d, S, 1)
    W = np.where(dk <= city.access_radius + off, city.access_w * dk / city.walk_kmh * 60, INF)
    return S, W


def evaluate(city, lines, pen, fleet, detail=False):
    """Returns fitness (person-hours/day of generalized time saved vs driving); -inf if infeasible."""
    lines = [ln for ln in lines if len(ln) >= 2]
    if not lines:
        return (0.0, {}) if detail else 0.0
    res = station_times(city, lines, pen, fleet, want_pred=detail)
    if res is None:
        return (-INF, {}) if detail else -INF
    idx, hw, T, pred, nid = res
    P = len(idx)
    S, W = _access(city, lines, idx)
    Bm = (W[:, :, None] + T[S]).min(1)                      # Z x P
    B_arg = (W[:, :, None] + T[S]).argmin(1) if detail else None
    Cm = (Bm[:, S] + W[None]).min(-1)                       # Z x Z
    fin = np.isfinite(Cm)
    dlt = np.where(fin, city.drive_t - Cm, 0.0)
    u = BETA * dlt + city.asc_z[:, None]
    p = np.where(fin, 1.0 / (1.0 + np.exp(-u)), 0.0)
    sav = city.flow * p * dlt
    fit = float(sav.sum() / 60.0)
    if not detail:
        return fit
    # --- detail metrics ---
    riders = city.flow * p
    tot_r = riders.sum()
    C_arg = (Bm[:, S] + W[None]).argmin(-1)
    X = transfer_counts(P, pred, nid)
    sp = np.take_along_axis(S, C_arg, 1)                    # dest station index for each (z, z')
    so = B_arg[np.arange(city.Z)[:, None], sp]
    so_st = S[np.arange(city.Z)[:, None], so]
    xf = X[so_st, sp]
    cyc, comps = _cyclomatic(lines, idx)
    cnt = {}
    for ln in lines:
        for c in set(ln):
            cnt[c] = cnt.get(c, 0) + 1
    kms = [sum(city.seg(a, b)[0] for a, b in zip(ln[:-1], ln[1:])) for ln in lines]
    water_hops = sum(1 for ln in lines for a, b in zip(ln[:-1], ln[1:]) if city.seg(a, b)[1] > 0.3)
    long_i = int(np.argmax(kms))
    pts = np.array([city.xy(*c) for c in lines[long_i]])
    ang = 0.0
    if len(pts) > 1:
        u_, s_, vt = np.linalg.svd(pts - pts.mean(0), full_matrices=False)
        ang = math.degrees(math.atan2(vt[0][1], vt[0][0])) % 180
    cov = (city.zH[(W.min(1) < INF)]).sum() / city.zH.sum()
    m = dict(
        fitness=fit, cost=net_cost(city, lines), n_lines=len(lines), n_stations=P,
        track_km=float(sum(kms)), mean_stations_per_line=float(np.mean([len(l) for l in lines])),
        max_stations_per_line=max(len(l) for l in lines), longest_line_km=float(max(kms)),
        longest_share=float(max(kms) / sum(kms)),
        n_transfer_stations=sum(1 for v in cnt.values() if v >= 2),
        transfers_per_line=sum(1 for v in cnt.values() if v >= 2) / len(lines),
        mean_line_multiplicity=float(np.mean(list(cnt.values()))),
        rings=int(cyc), components=int(comps),
        mode_share=float(tot_r / city.total_trips), riders=float(tot_r),
        mean_headway=float(np.mean(hw)), mean_transfers_per_trip=float((riders * xf).sum() / max(tot_r, 1e-9)),
        share_flow_transit_faster=float(city.flow[fin & (dlt > 0)].sum() / city.total_trips),
        water_hops=int(water_hops), dominant_angle=float(ang), coverage=float(cov),
        mean_saved_min_per_rider=float(sav.sum() / max(tot_r, 1e-9)),
    )
    return fit, m
