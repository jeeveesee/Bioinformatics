from collections import defaultdict, Counter, deque

# Function to read the adjacency list and construct an Eulerian cycle
def eulerian_cycle(adj_list):
    # -------------------------
    # 1) Parse adjacency safely
    # -------------------------
    graph = defaultdict(list)
    nodes = set()

    for line in adj_list:
        raw = line.strip()
        if not raw:
            continue
        if ":" not in raw:
            raise ValueError(f"Bad line (missing ':'): {raw}")
        left, right = raw.split(":", 1)
        u = left.strip()                     # keep nodes as strings (like your code)
        right = right.strip()
        if right:
            vs = right.split()               # flexible to multiple spaces
            graph[u].extend(vs)
            nodes.add(u)
            nodes.update(vs)
        else:
            nodes.add(u)
            graph[u]                         # ensure key exists

    # Ensure all nodes show up as keys (even if outdegree 0)
    for v in nodes:
        graph[v]  # touch to create empty list if missing

    # -----------------------------------------
    # 2) Degree checks for Eulerian *cycle*
    #    (in-degree == out-degree for every node)
    # -----------------------------------------
    indeg = Counter()
    outdeg = Counter()
    total_edges = 0

    for u, nbrs in graph.items():
        outdeg[u] += len(nbrs)
        total_edges += len(nbrs)
        for v in nbrs:
            indeg[v] += 1

    if total_edges == 0:
        raise ValueError("No edges: a nontrivial Eulerian cycle does not exist.")

    unbalanced = [u for u in graph if indeg[u] != outdeg[u]]
    if unbalanced:
        details = ", ".join(f"{u}(out={outdeg[u]}, in={indeg[u]})" for u in unbalanced)
        raise ValueError(f"Indegree/outdegree mismatch (cycle requires equality): {details}")

    # ---------------------------------------------------
    # 3) Weak connectivity on the underlying undirected graph
    #    (among nodes that have degree > 0)
    # ---------------------------------------------------
    active = {u for u in graph if outdeg[u] > 0 or indeg[u] > 0}
    # Build undirected view
    undirected = defaultdict(set)
    for u, nbrs in graph.items():
        for v in nbrs:
            undirected[u].add(v)
            undirected[v].add(u)

    # BFS over the undirected view starting from any active node
    start_for_bfs = next(iter(active))
    seen = set()
    q = deque([start_for_bfs])
    while q:
        x = q.popleft()
        if x in seen:
            continue
        seen.add(x)
        for y in undirected[x]:
            if y not in seen:
                q.append(y)

    if seen != active:
        missing = sorted(active - seen)
        raise ValueError(f"Graph is not (weakly) edge-connected. Unreachable nodes: {missing}")

    # ---------------------------------------------------
    # 4) Hierholzer's algorithm (your original structure)
    # ---------------------------------------------------
    # Work on a copy we can mutate; reverse lists so .pop() uses left-to-right input order
    G = {u: list(nbrs)[::-1] for u, nbrs in graph.items()}

    # Start from any node with outgoing edges (valid for Eulerian cycles)
    start_node = next((u for u in graph if G[u]), next(iter(graph)))

    cycle = []
    stack = []
    current_node = start_node  # Start from a valid node with outgoing edges

    # Main loop: walk edges, never reusing an edge (we pop it when used)
    while stack or G[current_node]:
        if not G[current_node]:
            cycle.append(current_node)  # dead-end: output node and backtrack
            current_node = stack.pop()
        else:
            stack.append(current_node)  # remember where we came from
            next_node = G[current_node].pop()  # consume one unused edge
            current_node = next_node

    cycle.append(current_node)
    cycle.reverse()  # Reverse to get the correct direction

    # 5) Safety: verify all edges were used
    if len(cycle) != total_edges + 1:
        raise ValueError("Not all edges were used—graph might not be properly connected.")

    return "->".join(cycle)
