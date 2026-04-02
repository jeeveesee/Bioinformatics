# from collections import Counter

# # Define the counters
# out_degrees = Counter({'3': 2, '6': 2, '9': 2, '0': 1, '1': 1, '2': 1, '7': 1, '8': 1})
# in_degrees = Counter({'3': 2, '6': 2, '2': 1, '1': 1, '0': 1, '4': 1, '7': 1, '8': 1, '9': 1})

# # Get all unique keys from both counters to ensure no node is missed.
# # A Counter will return 0 for a key that is not explicitly present,
# # which simplifies checking the difference.
# # all_keys = set(out_degrees.keys()) | set(in_degrees.keys())

# # List to store the keys that meet the condition
# keys_with_out_greater = []

# for key in out_degrees.keys():
#     # Check if the out-degree value is exactly 1 greater than the in-degree value
#     if out_degrees[key] == in_degrees[key] + 1:
#         keys_with_out_greater.append(key)

# print(f"Keys where out-degrees is 1 greater than in-degrees: {keys_with_out_greater}")

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

from collections import defaultdict, Counter

adj = defaultdict(list)
out_degrees = Counter()
in_degrees = Counter()
all_nodes = set()


adj_list = """
0: 2
1: 3
2: 1
3: 0 4
6: 3 7
7: 8
8: 9
9: 6
"""
champaran = adj_list.strip().split("\n")
champaran_dict = {}
for line in champaran:
    node, edges = line.split(": ")
    champaran_dict[node] = edges.split(" ")
print(f"\n\nChamparan Dict: {champaran_dict}\n\n")

# Process input for graph and degrees
for u, neighbors in champaran_dict.items():
    # Using reversed list for adj to make pop() easier to manage in the recursive step,
    # but for a true Hierholzer's, order doesn't strictly matter. 
    # Here, we'll keep it as is, and use the simpler approach: pop() from the end.
    print(f"\n{u=} and {neighbors=}")
    print(adj[u].extend(neighbors))
    print(f"adj[u]: {adj[u]}")
    out_degrees[u] += len(neighbors)
    print(f"out_degrees[u]: {out_degrees[u]}")
    all_nodes.add(u)
    print(f"add[u]: {all_nodes}")
    
    for v in neighbors:
        in_degrees[v] += 1
        all_nodes.add(v)

print("\n\n")
print(f"{all_nodes=}")
print(f"{in_degrees=}")
print(f"{out_degrees=}")