# H3-Route: A novel priority-aware online routing algorithm
## for hazard zone relief Robots and Drones

## Background

In this project, I invented a novel, priority-aware online routing algorithm specifically designed for hazard zone relief operations using autonomous robots and drones. This heuristic addresses the challenges of dynamic traffic requests by aiming for bounded latency, route stability, and near-optimal total travel distance, which are crucial for battery- or fuel-limited vehicles. The core of H3-Route lies in creating a 3-D "priority-time warp" using a Hilbert space-filling curve to combine geographical coordinates with urgency metadata for efficient index-based ordering and quick online insertion of new requests. The algorithm also incorporates a frozen horizon of imminent stops and periodic local optimization to ensure route stability and overall efficiency, outperforming traditional greedy and FIFO approaches in simulation experiments.

## The Idea

The H3-Route algorithm (Hilbert, Horizon, Heuristic) defines its online routing strategy for hazard zone relief robots and drones through three integrated technical innovations:

### 1. 3-D “Priority-Time Warp” using the Hilbert Space-Filling Curve (H)

This innovation creates a locality-preserving backbone for managing dynamic requests by turning the multi-dimensional routing complexity (travel cost, priority, and time-window penalties) into a single ordered index thus able to do quick range query or updates.

- **Multi-Dimensional Key**: The algorithm utilizes the Hilbert space-filling curve (SFV) to preserve locality. It normalizes geo-information (X, Y) and converts priority and deadlines into a scalar "urgency" (U).
- **Total Ordering**: These components are combined to compute a single 3D Hilbert index (D). This key provides a total order that is simultaneously geometrically local (via X,Y) and urgency-aware (via U).
- **Dynamic Sequencing**: Sorting requests by D creates the "priority-time warp" that allows urgent nodes to "bubble up" in the sequence without requiring a global re-solve.

### 2. Efficient Online Insertion with Priority-Aware Regret (H)

H3-Route utilizes a fast insertion heuristic to handle spikes of arriving requests with sublinear complexity.

- **Data Structure and Query**: The current route order is stored in a balanced order tree (e.g., augmented skip-list or red-black tree), keyed by the Hilbert index (h). This structure allows the proximal neighbor of a new request to be found quickly in O(logn).
- **Regret Scoring**: The insertion decision is made using a priority-aware regret score over a small set (k) of nearby candidate edges (u,v). This objective function balances three factors:

   1. **1. Distance/Energy Cost**: The added travel cost.
   1. **2. Priority-Weighted Lateness**: A penalty weighted by the new request's priority (w(p)) if the predicted arrival time (ETA) exceeds the deadline (dr).
   1. **3. Stability Penalty**: A term that discourages changes inside the frozen prefix of length F.

- **Complexity**: The time complexity per insertion is efficient, at O(\log(n)+k)), which is effectively near constant when k (the number of scored candidates) is small.

### 3. Periodic Receding Horizon Optimization and Stability Prefix (H)

This approach maintains route stability while achieving near-optimal performance by limiting route modification.

- **Stability Mechanism**: A frozen prefix of length F (the next F imminent stops) is maintained. This prevents the route from "thrashing geographically" during request spikes and ensures stability, which is a practical requirement for safe operation.
- **Local Optimization**: Full, costly global optimization is avoided. Instead, a novel local optimization heuristic is run occasionally (e.g., every τ inserts).
- **Horizon Window**: Optimization is strictly bounded to a horizon window defined immediately following the frozen prefix.
- **Methodology**: Within this window, bounded local search (such as 2-opt, 3-opt, or Or-opt) is executed with a small time cap to remove detours. The window "recedes" forward for the next optimization. This periodical cleanup helps achieve near-static-TSP quality without needing a full re-solve.

## The impact

If the H3-Route algorithm is proven successful and deployed, its efficiency, stability, and priority awareness can significantly escalate the social impact of disaster response and urgent logistics.

### 1. Saving Lives Through Timely and Efficient Response

The algorithm directly confronts the problem of people being left without timely assistance following disasters.

- **Priority-Aware Urgency**: H3-Route utilizes a 3-D "priority-time warp" and a soft partial order mechanism to ensure urgent stops rise in the sequence without requiring a computationally heavy global re-solve. This ability to prioritize based on urgency and delivery deadlines is paramount in life-saving missions.
- **Fuel/Battery Conservation**: By delivering near shortest routing performance (even outperforming greedy methods in some simulations), H3-Route optimizes geographical efficiency. This is critical for the limited service route of a robot/drone restricted by fuel or battery. Conserving energy allows the vehicles to assist more victims or stay operational longer in the hazard zone.

### 2. Ensuring Stability for Reliable Operations

The stability features of H3-Route directly translate into operational reliability, which is essential when lives are at stake.

- **Reliable Control**: The heuristic maintains a frozen prefix of imminent stops (F). This feature limits route churn and prevents the system from "thrashing" during spikes in request arrivals, which is a practical requirement for safe flight and drive controllers during real applications. Reliable operation enhances trust and safety in critical relief missions.
- **Bounded Latency**: The algorithm operates with O(\log(n)+k) efficiency per insert. Since decisions are made quickly and efficiently, aid is not delayed by computational bottlenecks, ensuring that the service remains bounded-latency.

### 3. Scaling Impact Through Fleet Coordination

The social impact can be vastly scaled up by extending the algorithm to handle multiple vehicles, a planned area for future work.

- **Coordinated Fleet Management**: Future implementation aims to extend the method to a coordinated fleet of UAVs or drones. This extension is vital for large-scale disaster scenarios, allowing resources to be efficiently distributed across a wider geographical area.
- **Optimized Assignment and Handoffs**: The multi-vehicle heuristic will manage assignments using a partition method in the Hilbert index space and study cross-vehicle handoffs. This level of coordination will ensure that relief efforts are highly optimized, guaranteeing that aid reaches its destination via the quickest and most cost-effective vehicle, ultimately commanding a full fleet ready for hazard zone recovery and relief efforts

