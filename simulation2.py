# Extended simulator with deadlines/time windows and priority-weighted lateness (P50/P95/P99)
# - Keeps H3-Route priority-aware; GreedyNearest ignores priority/deadlines; FIFO appends.
# - Plots P50, P95, P99 curves of priority-weighted lateness as requests stream in.
# - Also prints a concise summary table and saves it to CSV.

import math, random, time
from dataclasses import dataclass
from bisect import bisect_left
from typing import List, Tuple, Dict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from caas_jupyter_tools import display_dataframe_to_user

random.seed(1234)
np.random.seed(1234)

# ----------------------------
# Hilbert (2D) utilities (fast)
# ----------------------------
def _rot(n, x, y, rx, ry):
    if ry == 0:
        if rx == 1:
            x = n - 1 - x
            y = n - 1 - y
        x, y = y, x
    return x, y

def hilbert_xy2d(order_pow2:int, x:int, y:int) -> int:
    n = 1 << order_pow2
    d = 0
    s = n >> 1
    xx, yy = x, y
    while s > 0:
        rx = 1 if (xx & s) else 0
        ry = 1 if (yy & s) else 0
        d += s * s * ((3 * rx) ^ ry)
        xx, yy = _rot(s, xx, yy, rx, ry)
        s >>= 1
    return d

def hilbert_key_2d(xf: float, yf: float, order_pow2:int=16) -> int:
    eps = 1e-12
    xf = min(max(xf, 0.0), 1.0 - eps)
    yf = min(max(yf, 0.0), 1.0 - eps)
    n = 1 << order_pow2
    xi = int(xf * n)
    yi = int(yf * n)
    return hilbert_xy2d(order_pow2, xi, yi)

# ----------------------------
# Core structures & helpers
# ----------------------------
@dataclass
class Request:
    rid: int
    x: float
    y: float
    priority: int        # 1 (highest) .. 3 (lowest)
    arrival_time: float  # when it appears to the algorithm
    deadline: float      # absolute time by which service is desired

DEPOT = (0.5, 0.5)
SPEED = 1.0      # units / time
SERVICE = 0.05   # service time per stop

def euclid(a:Tuple[float,float], b:Tuple[float,float]) -> float:
    return math.hypot(a[0]-b[0], a[1]-b[1])

def travel_time(a:Tuple[float,float], b:Tuple[float,float]) -> float:
    return euclid(a,b) / SPEED

def priority_weight(p:int) -> int:
    return {1:4, 2:2, 3:1}.get(p,1)

def compute_etas(route: List[Request], start_time: float = 0.0) -> List[float]:
    """Compute ETA for each stop in current order, starting at start_time from the depot."""
    etas = []
    t = start_time + travel_time(DEPOT, (route[0].x, route[0].y)) if route else start_time
    for i, r in enumerate(route):
        if i == 0:
            t = start_time + travel_time(DEPOT, (r.x, r.y))
        else:
            prev = route[i-1]
            t += SERVICE + travel_time((prev.x, prev.y), (r.x, r.y))
        etas.append(t)
    return etas

def lateness_metrics(route: List[Request], etas: List[float]) -> Tuple[float,float,float]:
    """Return P50, P95, P99 of priority-weighted lateness across current route."""
    if not route:
        return 0.0, 0.0, 0.0
    vals = []
    for r, eta in zip(route, etas):
        late = max(0.0, eta - r.deadline)
        vals.append(priority_weight(r.priority) * late)
    arr = np.array(vals)
    return np.percentile(arr, 50), np.percentile(arr, 95), np.percentile(arr, 99)

# ----------------------------
# Generating online requests with deadlines
# ----------------------------
def generate_requests_with_deadlines(n:int=100, interarrival:float=0.5) -> List[Request]:
    reqs = []
    for i in range(n):
        # Geometry
        x, y = float(np.random.rand()), float(np.random.rand())
        # Priority mix
        p = int(np.random.choice([1,2,3], p=[0.25,0.35,0.40]))
        # Arrival time
        arrival = i * interarrival + float(np.random.rand())*0.05  # small jitter
        # Baseline difficulty: deadline = arrival + base * distance_from_depot + slack(priority)
        base_dist = euclid((x,y), DEPOT)
        # Priority-dependent slack (tight for P1)
        slack = {1: 0.35, 2: 0.9, 3: 1.6}[p]
        # Add small randomness so instances differ
        slack *= (0.9 + 0.2*np.random.rand())
        # Slightly superlinear travel scaling to make long trips harder
        deadline = arrival + 1.2*base_dist + slack
        reqs.append(Request(i, x, y, p, arrival, deadline))
    return reqs

# ----------------------------
# Algorithms
# ----------------------------
class H3Route:
    """Priority-aware: distance + disruption + *lateness penalty* in insertion scoring; Hilbert backbone with urgency warp."""
    def __init__(self, freeze_F:int=6, k:int=16, order_pow2:int=14, warp_alpha:int=2):
        self.route: List[Request] = []
        self.sorted_hilbert: List[Tuple[int, int]] = []
        self.index_in_route: Dict[int, int] = {}
        self.freeze_F = freeze_F
        self.k = k
        self.order_pow2 = order_pow2
        self.warp_alpha = warp_alpha
        self.prefix_times: List[float] = []  # ETAs for current route from t=0

    def _warped_key(self, r: Request) -> int:
        base = hilbert_key_2d(r.x, r.y, self.order_pow2)
        # urgency warp: push high-priority earlier
        warp = self.warp_alpha * 10**6 * priority_weight(r.priority)
        return max(0, base - warp)

    def _insert_sorted_hilbert(self, r: Request, rid_idx:int):
        key = self._warped_key(r)
        pos = bisect_left(self.sorted_hilbert, (key, -1))
        self.sorted_hilbert.insert(pos, (key, rid_idx))

    def _update_indices_after(self, start:int, delta:int):
        for rid, idx in list(self.index_in_route.items()):
            if idx >= start:
                self.index_in_route[rid] = idx + delta
        for i in range(len(self.sorted_hilbert)):
            key, idx = self.sorted_hilbert[i]
            if idx >= start:
                self.sorted_hilbert[i] = (key, idx + delta)

    def _candidate_edges_from_neighbors(self, new_r: Request) -> List[int]:
        key = self._warped_key(new_r)
        pos = bisect_left(self.sorted_hilbert, (key, -1))
        idxs = []
        lo = max(0, pos - self.k//2)
        hi = min(len(self.sorted_hilbert), pos + self.k//2)
        neighbors = self.sorted_hilbert[lo:hi]
        for _, route_idx in neighbors:
            idxs.append(route_idx)
        idxs.append(len(self.route)-1)
        idxs = sorted(set(i for i in idxs if i >= -1))
        return idxs

    def _disruption_penalty(self, insert_after_idx:int) -> float:
        return 3.0 if insert_after_idx < self.freeze_F - 1 else 0.0

    def _compute_prefix_etas(self):
        self.prefix_times = compute_etas(self.route, 0.0)

    def _eta_if_inserted_after(self, idx:int, r:Request) -> float:
        # ETA to r if inserted after position idx in current route (starting from t=0)
        if idx < 0:
            t_prev = 0.0
            prev_pt = DEPOT
        else:
            t_prev = self.prefix_times[idx]  # arrival time at u
            prev_pt = (self.route[idx].x, self.route[idx].y)
        t_prev += (SERVICE if idx >= 0 else 0.0)
        return t_prev + travel_time(prev_pt, (r.x, r.y))

    def insert(self, r: Request):
        if not self.route:
            self.route.append(r)
            self.index_in_route[r.rid] = 0
            self._insert_sorted_hilbert(r, 0)
            self._compute_prefix_etas()
            return

        # ensure prefix times ready
        self._compute_prefix_etas()

        candidates = self._candidate_edges_from_neighbors(r)
        best_cost = float('inf')
        best_pos = len(self.route)-1
        for idx in candidates:
            # distance delta
            if idx < 0:
                prev_pt = DEPOT
            else:
                prev_pt = (self.route[idx].x, self.route[idx].y)
            if idx + 1 < len(self.route):
                nxt_pt = (self.route[idx+1].x, self.route[idx+1].y)
            else:
                nxt_pt = None
            add = euclid(prev_pt, (r.x, r.y))
            if nxt_pt is not None:
                add += euclid((r.x, r.y), nxt_pt) - euclid(prev_pt, nxt_pt)
            dist_term = add

            # lateness penalty for r at this position
            eta_r = self._eta_if_inserted_after(idx, r)
            late_r = max(0.0, eta_r - r.deadline)
            lateness_term = 10.0 * priority_weight(r.priority) * late_r  # μ=10

            # disruption penalty
            disrupt = self._disruption_penalty(idx)

            cost = dist_term + lateness_term + disrupt
            if cost < best_cost:
                best_cost = cost
                best_pos = idx

        insert_at = best_pos + 1
        self.route.insert(insert_at, r)
        self._update_indices_after(insert_at, +1)
        self.index_in_route[r.rid] = insert_at
        self._insert_sorted_hilbert(r, insert_at)

        # tiny local improvement (2-opt) around insertion with time cap
        self._local_two_opt(center=insert_at, window=14, time_budget_ms=6)
        self._compute_prefix_etas()

    def _local_two_opt(self, center:int, window:int=14, time_budget_ms:int=6):
        if len(self.route) < 4:
            return
        n = len(self.route)
        a = max(0, center - window//2)
        b = min(n-1, center + window//2)
        start_time = time.perf_counter()
        improved = True
        while improved and (time.perf_counter() - start_time) * 1000.0 < time_budget_ms:
            improved = False
            for i in range(a, b-2):
                for j in range(i+2, b):
                    pA = DEPOT if i == -1 else (self.route[i].x, self.route[i].y)
                    pB = (self.route[i+1].x, self.route[i+1].y)
                    pC = (self.route[j].x, self.route[j].y)
                    pD = (self.route[j+1].x, self.route[j+1].y) if j+1 < n else None
                    before = euclid(pA, pB) + (euclid(pC, pD) if pD else 0.0)
                    after  = euclid(pA, pC) + (euclid(pB, pD) if pD else 0.0)
                    if after + 1e-12 < before:
                        self.route[i+1:j+1] = reversed(self.route[i+1:j+1])
                        improved = True

class GreedyNearestInsertion:
    """Distance-only; ignores priority and deadlines."""
    def __init__(self, freeze_F:int=6):
        self.route: List[Request] = []
        self.freeze_F = freeze_F

    def _delta_distance(self, insert_after_idx:int, r:Request) -> float:
        if insert_after_idx < 0:
            prev_pt = DEPOT
        else:
            prev = self.route[insert_after_idx]
            prev_pt = (prev.x, prev.y)
        if insert_after_idx + 1 < len(self.route):
            nxt = self.route[insert_after_idx+1]
            nxt_pt = (nxt.x, nxt.y)
        else:
            nxt_pt = None
        add = euclid(prev_pt, (r.x, r.y))
        if nxt_pt is not None:
            add += euclid((r.x, r.y), nxt_pt) - euclid(prev_pt, nxt_pt)
        return add

    def insert(self, r: Request):
        if not self.route:
            self.route.append(r)
            return
        best_cost = float('inf')
        best_pos = len(self.route)-1
        for idx in range(-1, len(self.route)):
            penalty = 3.0 if idx < self.freeze_F - 1 else 0.0
            delta = self._delta_distance(idx, r) + penalty
            if delta < best_cost:
                best_cost = delta
                best_pos = idx
        self.route.insert(best_pos+1, r)

class FifoAppend:
    """Arrival order only."""
    def __init__(self):
        self.route: List[Request] = []
    def insert(self, r: Request):
        self.route.append(r)

# ----------------------------
# Simulation runner
# ----------------------------
def run_with_metrics(reqs: List[Request], alg, name:str):
    """Return dict of curves: route length, lateness p50/p95/p99, compute time, final summary."""
    route_lengths, p50s, p95s, p99s = [], [], [], []
    start = time.perf_counter()
    for r in reqs:
        alg.insert(r)
        # recompute metrics after each insert
        route = alg.route
        if route:
            # route length as time from depot (distance proxy)
            rl = 0.0
            if route:
                rl = euclid(DEPOT, (route[0].x, route[0].y))
                for i in range(len(route)-1):
                    rl += euclid((route[i].x, route[i].y), (route[i+1].x, route[i+1].y))
            etas = compute_etas(route, 0.0)
            p50, p95, p99 = lateness_metrics(route, etas)
        else:
            rl, p50, p95, p99 = 0.0, 0.0, 0.0, 0.0
        route_lengths.append(rl)
        p50s.append(p50); p95s.append(p95); p99s.append(p99)
    elapsed_ms = (time.perf_counter() - start) * 1000.0
    return {
        "name": name,
        "route_lengths": route_lengths,
        "p50": p50s, "p95": p95s, "p99": p99s,
        "time_ms": elapsed_ms,
        "final_length": route_lengths[-1],
        "final_p50": p50s[-1], "final_p95": p95s[-1], "final_p99": p99s[-1],
    }

# Generate data
N_REQ = 100
requests = generate_requests_with_deadlines(N_REQ, interarrival=0.5)

# Instantiate algorithms
h3 = H3Route(freeze_F=6, k=16, order_pow2=14, warp_alpha=2)
gni = GreedyNearestInsertion(freeze_F=6)
fifo = FifoAppend()

# Run
res_h3 = run_with_metrics(requests, h3, "H3-Route (priority-aware)")
res_gn = run_with_metrics(requests, gni, "Greedy Nearest (distance-only)")
res_ff = run_with_metrics(requests, fifo, "FIFO Append (arrival order)")

# Build summary table
summary = pd.DataFrame({
    "Algorithm": [res_h3["name"], res_gn["name"], res_ff["name"]],
    "Final route length": [res_h3["final_length"], res_gn["final_length"], res_ff["final_length"]],
    "Final priority-weighted lateness P50": [res_h3["final_p50"], res_gn["final_p50"], res_ff["final_p50"]],
    "Final priority-weighted lateness P95": [res_h3["final_p95"], res_gn["final_p95"], res_ff["final_p95"]],
    "Final priority-weighted lateness P99": [res_h3["final_p99"], res_gn["final_p99"], res_ff["final_p99"]],
    "Total compute time (ms)": [res_h3["time_ms"], res_gn["time_ms"], res_ff["time_ms"]]
})

display_dataframe_to_user("Routing with Deadlines — Summary", summary)

# ----------------------------
# Plot: P50, P95, P99 lateness
# ----------------------------
x = np.arange(1, N_REQ+1)

plt.figure()
plt.plot(x, res_h3["p50"], label=res_h3["name"])
plt.plot(x, res_gn["p50"], label=res_gn["name"])
plt.plot(x, res_ff["p50"], label=res_ff["name"])
plt.xlabel("Number of requests inserted")
plt.ylabel("Priority-weighted lateness P50")
plt.title("P50 lateness vs Online Insertions")
plt.legend()
plt.show()

plt.figure()
plt.plot(x, res_h3["p95"], label=res_h3["name"])
plt.plot(x, res_gn["p95"], label=res_gn["name"])
plt.plot(x, res_ff["p95"], label=res_ff["name"])
plt.xlabel("Number of requests inserted")
plt.ylabel("Priority-weighted lateness P95")
plt.title("P95 lateness vs Online Insertions")
plt.legend()
plt.show()

plt.figure()
plt.plot(x, res_h3["p99"], label=res_h3["name"])
plt.plot(x, res_gn["p99"], label=res_gn["name"])
plt.plot(x, res_ff["p99"], label=res_ff["name"])
plt.xlabel("Number of requests inserted")
plt.ylabel("Priority-weighted lateness P99")
plt.title("P99 lateness vs Online Insertions")
plt.legend()
plt.show()

# Save summary
summary.to_csv("/mnt/data/routing_deadlines_summary.csv", index=False)
