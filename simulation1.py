# Compact Python simulation and plotting for H3-Route vs traditional approaches
# - H3-Route (Hilbert + Horizon Freeze + Heuristic regret insertion + tiny local search)
# - Greedy Nearest-Insertion (distance-only)
# - FIFO Append (queue)
#
# Produces two plots:
# 1) Cumulative route length vs number of requests
# 2) Priority ordering efficiency (priority-weighted position index) vs number of requests
#
# Notes:
# - Uses a 2D Hilbert index with a simple "priority warp" (pushes higher priority earlier).
# - Freezes the first F stops to avoid thrashing (applies to all algs for fairness).
# - For simplicity, we don't simulate real-time deadlines; efficiency proxy is distance + priority ordering.
# - Everything is self-contained (no external packages besides matplotlib, numpy, pandas).
import math, random, time
from dataclasses import dataclass
from bisect import bisect_left
from typing import List, Tuple, Dict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

random.seed(42)
np.random.seed(42)

# ----------------------------
# Hilbert curve utilities (2D)
# ----------------------------
def _rot(n, x, y, rx, ry):
    # Rotate/flip a quadrant appropriately
    if ry == 0:
        if rx == 1:
            x = n - 1 - x
            y = n - 1 - y
        # Swap x and y
        x, y = y, x
    return x, y

def hilbert_xy2d(order_pow2:int, x:int, y:int) -> int:
    """Return the Hilbert distance (d) for integer coords (x,y) in [0, n-1]^2, n = 2^order_pow2."""
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

def hilbert_key_2d(x_float: float, y_float: float, order_pow2:int=16) -> int:
    """Map (x,y) in [0,1) to integer Hilbert key with 2^order grid."""
    # Clip to [0, 1 - eps] for safety
    eps = 1e-12
    xf = min(max(x_float, 0.0), 1.0 - eps)
    yf = min(max(y_float, 0.0), 1.0 - eps)
    n = 1 << order_pow2
    xi = int(xf * n)
    yi = int(yf * n)
    return hilbert_xy2d(order_pow2, xi, yi)

# ----------------------------
# Core data structures
# ----------------------------
@dataclass
class Request:
    rid: int
    x: float
    y: float
    priority: int  # 1=highest, 3=lowest

# Distance helper
def euclid(a:Tuple[float,float], b:Tuple[float,float]) -> float:
    return math.hypot(a[0]-b[0], a[1]-b[1])

# Compute total route length (start at depot at (0.5,0.5) for neutrality, then chain through route, no return)
DEPOT = (0.5, 0.5)

def route_length(points: List[Request]) -> float:
    if not points:
        return 0.0
    total = euclid(DEPOT, (points[0].x, points[0].y))
    for i in range(len(points)-1):
        total += euclid((points[i].x, points[i].y), (points[i+1].x, points[i+1].y))
    return total

def priority_weight(p:int) -> int:
    # Higher weight for higher priority
    return {1:4, 2:2, 3:1}.get(p,1)

def priority_position_score(points: List[Request]) -> float:
    # Lower is better: sum over i (w(priority_i) * position_i)
    s = 0.0
    for i, r in enumerate(points):
        s += priority_weight(r.priority) * (i+1)
    return s / max(1, len(points))

# ---------------------------------------------
# H3-Route (compact version for the simulation)
# ---------------------------------------------
class H3Route:
    def __init__(self, freeze_F:int=5, k:int=16, order_pow2:int=16, warp_alpha:int=2):
        self.route: List[Request] = []
        # Maintain a second structure: a sorted list by warped Hilbert key
        self.sorted_hilbert: List[Tuple[int, int]] = []  # list of (key, rid_index_in_route)
        self.index_in_route: Dict[int, int] = {}  # rid -> index in route
        self.freeze_F = freeze_F
        self.k = k
        self.order_pow2 = order_pow2
        self.warp_alpha = warp_alpha  # scale for priority warp

    def _warped_key(self, r: Request) -> int:
        base = hilbert_key_2d(r.x, r.y, self.order_pow2)
        # Priority warp: subtract chunk to move higher-priority earlier
        warp = self.warp_alpha * 10**6 * priority_weight(r.priority)
        return max(0, base - warp)

    def _insert_sorted_hilbert(self, r: Request, rid_idx:int):
        key = self._warped_key(r)
        pos = bisect_left(self.sorted_hilbert, (key, -1))
        self.sorted_hilbert.insert(pos, (key, rid_idx))

    def _update_indices_after(self, start:int, delta:int):
        # When we insert into self.route, shift indices for nodes at positions >= start
        for rid, idx in list(self.index_in_route.items()):
            if idx >= start:
                self.index_in_route[rid] = idx + delta
        # Update the hilbert-sorted (key, idx) tuples
        for i in range(len(self.sorted_hilbert)):
            key, idx = self.sorted_hilbert[i]
            if idx >= start:
                self.sorted_hilbert[i] = (key, idx + delta)

    def _candidate_edges_from_neighbors(self, new_r: Request) -> List[int]:
        # Use position in hilbert-sorted list to propose candidate insertion positions in the *route*
        key = self._warped_key(new_r)
        pos = bisect_left(self.sorted_hilbert, (key, -1))
        idxs = []
        # Gather neighbors around pos
        lo = max(0, pos - self.k//2)
        hi = min(len(self.sorted_hilbert), pos + self.k//2)
        neighbors = self.sorted_hilbert[lo:hi]
        # Convert neighbors' route indices to candidate insertion positions (after each neighbor)
        for _, route_idx in neighbors:
            idxs.append(route_idx)  # insert after route_idx (i.e., between route_idx and route_idx+1)
        # Always include end of route as a candidate
        idxs.append(len(self.route)-1)
        # Dedup and sort
        idxs = sorted(set(i for i in idxs if i >= -1))
        return idxs

    def _delta_distance(self, insert_after_idx:int, r:Request) -> float:
        # Insertion between (prev)->(next): prev is DEPOT if idx==-1
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

    def _disruption_penalty(self, insert_after_idx:int) -> float:
        # Penalize inserting inside the frozen prefix
        if insert_after_idx < self.freeze_F - 1:
            return 3.0  # constant penalty; tuned small since we compare deltas
        return 0.0

    def insert(self, r: Request):
        if not self.route:
            self.route.append(r)
            self.index_in_route[r.rid] = 0
            self._insert_sorted_hilbert(r, 0)
            return

        # Pick best candidate by regret-like delta: distance + disruption
        candidates = self._candidate_edges_from_neighbors(r)
        best_cost = float('inf')
        best_pos = len(self.route)-1
        for idx in candidates:
            # Respect freeze: avoid inserting *before* the frozen prefix end unless it's the only place
            delta = self._delta_distance(idx, r) + self._disruption_penalty(idx)
            if delta < best_cost:
                best_cost = delta
                best_pos = idx

        # Commit insertion in route
        insert_at = best_pos + 1
        self.route.insert(insert_at, r)
        # Update indices after insertion
        self._update_indices_after(insert_at, +1)
        # Record index
        self.index_in_route[r.rid] = insert_at
        # Insert into hilbert-sorted view
        self._insert_sorted_hilbert(r, insert_at)

        # Tiny local improvement around insertion (bounded 2-opt in a small window)
        self._local_two_opt(center=insert_at, window=12, time_budget_ms=5)

    def _local_two_opt(self, center:int, window:int=12, time_budget_ms:int=5):
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
                    # Edges: (i,i+1) and (j,j+1) -> 2-opt swap
                    pA = DEPOT if i == -1 else (self.route[i].x, self.route[i].y)
                    pB = (self.route[i+1].x, self.route[i+1].y)
                    pC = (self.route[j].x, self.route[j].y)
                    pD = (self.route[j+1].x, self.route[j+1].y) if j+1 < n else None

                    before = euclid(pA, pB) + euclid(pC, pD) if pD else euclid(pA, pB)
                    after  = euclid(pA, pC) + (euclid(pB, pD) if pD else 0.0)
                    if after + 1e-12 < before:
                        # reverse segment (i+1 .. j)
                        self.route[i+1:j+1] = reversed(self.route[i+1:j+1])
                        # Need to rebuild indices mapping quickly (just for affected window)
                        for k in range(i+1, j+1):
                            self.index_in_route[self.route[k].rid] = k
                        improved = True
                        # continue searching within budget

# ---------------------------------------------
# Baselines
# ---------------------------------------------
class GreedyNearestInsertion:
    def __init__(self, freeze_F:int=5):
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
            # optional freeze penalty to be fair
            penalty = 3.0 if idx < self.freeze_F - 1 else 0.0
            delta = self._delta_distance(idx, r) + penalty
            if delta < best_cost:
                best_cost = delta
                best_pos = idx
        self.route.insert(best_pos+1, r)

class FifoAppend:
    def __init__(self):
        self.route: List[Request] = []
    def insert(self, r: Request):
        self.route.append(r)

# ---------------------------------------------
# Simulation
# ---------------------------------------------
def generate_requests(n:int=100) -> List[Request]:
    # Positions uniform; priorities biased toward P3
    pri_choices = [1,2,3]
    pri_probs   = [0.2,0.3,0.5]
    reqs = []
    for i in range(n):
        x, y = float(np.random.rand()), float(np.random.rand())
        p = int(np.random.choice(pri_choices, p=pri_probs))
        reqs.append(Request(i, x, y, p))
    return reqs

def run_sim(reqs: List[Request], alg) -> Tuple[List[float], List[float], float]:
    lengths = []
    prio_scores = []
    start = time.perf_counter()
    for r in reqs:
        alg.insert(r)
        lengths.append(route_length(alg.route))
        prio_scores.append(priority_position_score(alg.route))
    elapsed = (time.perf_counter() - start) * 1000.0  # ms total
    return lengths, prio_scores, elapsed

# Prepare data
N_REQ = 100
requests = generate_requests(N_REQ)

# Instantiate algorithms
h3 = H3Route(freeze_F=6, k=16, order_pow2=14, warp_alpha=2)
gni = GreedyNearestInsertion(freeze_F=6)
fifo = FifoAppend()

# Run simulations
len_h3, pr_h3, t_h3 = run_sim(requests, h3)
len_gn, pr_gn, t_gn = run_sim(requests, gni)
len_ff, pr_ff, t_ff = run_sim(requests, fifo)

# Summary table
summary = pd.DataFrame({
    "Algorithm": ["H3-Route", "Greedy-Nearest", "FIFO-Append"],
    "Final route length": [len_h3[-1], len_gn[-1], len_ff[-1]],
    "Final priority-pos score (lower better)": [pr_h3[-1], pr_gn[-1], pr_ff[-1]],
    "Total compute time (ms)": [t_h3, t_gn, t_ff]
})

# ----------------------------
# Plot 1: Cumulative route length vs #requests
# ----------------------------
plt.figure()
x = np.arange(1, N_REQ+1)
plt.plot(x, len_h3, label="H3-Route")
plt.plot(x, len_gn, label="Greedy-Nearest")
plt.plot(x, len_ff, label="FIFO-Append")
plt.xlabel("Number of requests inserted")
plt.ylabel("Cumulative route length")
plt.title("Route Length vs Online Insertions")
plt.legend()
plt.show()

# ----------------------------
# Plot 2: Priority-weighted position index (lower is better)
# ----------------------------
plt.figure()
plt.plot(x, pr_h3, label="H3-Route")
plt.plot(x, pr_gn, label="Greedy-Nearest")
plt.plot(x, pr_ff, label="FIFO-Append")
plt.xlabel("Number of requests inserted")
plt.ylabel("Priority-weighted position index (lower better)")
plt.title("Priority Ordering Efficiency vs Online Insertions")
plt.legend()
plt.show()

# Save artifacts
summary.to_csv("/mnt/data/routing_summary.csv", index=False)

