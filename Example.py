import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from PathLengthWeightdistanceAlg_Class import Graph, calc_distance_matrices


'''
The following 6-node network is created by adding vertices with weights.
'''
g = Graph(6, weight=1)

g.add_edge(1, 0, 4)
g.add_edge(1, 2, 2)
g.add_edge(2, 5, 4)
g.add_edge(3, 0, 3)
g.add_edge(3, 4, 1)
g.add_edge(4, 1, 1)
g.add_edge(4, 5, 3)
g.add_edge(5, 0, 1)


'''
# If you have the proximity matrix, you can create the network as follows
M = np.array([[0, np.inf, np.inf, np.inf, np.inf, np.inf],
 [4, 0, 2, np.inf, np.inf, np.inf],
 [np.inf, np.inf, 0, np.inf, np.inf, 4],
 [ 3, np.inf, np.inf, 0, 1, np.inf],
 [np.inf, 1, np.inf, np.inf, 0, 3],
 [1, np.inf, np.inf, np.inf, np.inf, 0]])

g = Graph(int(M.shape[0]), weight=1)
g.add_edge_matrix(M)

'''


'''
Graph proximity matrix.
'''
print(g.edges)


'''
Weights associated with path lengths for path-length-weights distance calculation (W).
'''
print(g.weights)


'''
The graph is displayed using the networkx library.
'''
G = nx.DiGraph()

n_nodes = g.V
G.add_nodes_from(range(n_nodes))
for i in range(n_nodes):
    for j in range(n_nodes):
        w = g.edges[i][j]
        if w != 0 and w != np.inf:
            G.add_edge(i, j, weight=w)

fig, ax = plt.subplots(figsize=(10,10))

pos = nx.kamada_kawai_layout(G, weight='weight') # pos = nx.spring_layout(G), etc.
labels = {i: f'$v_{{{i}}}$' for i in G.nodes}
node_colors = 'lightblue' 

nx.draw(G, pos, labels=labels, with_labels=True, ax=ax, node_color=node_colors, node_size=2000, font_weight='bold', font_size=20)
weights = {(u, v): d['weight'] for u, v, d in G.edges(data=True)}
nx.draw_networkx_edge_labels(G, pos, edge_labels=weights, font_size=16)
plt.title('Graph with the weigths')
plt.show()


'''
The path-length-weights distances and dijkstra distances the of all the nodes of
 the graph are calculated using the function 'calc_distance_matrices'.
'''
plw_distances, dijkstra_distances = calc_distance_matrices(g, filt=0)

print(plw_distances)
print(dijkstra_distances)


'''
The path-length-weights distances of the graph is displayed.
'''
G_plw = nx.DiGraph()

n_nodes = g.V
G_plw.add_nodes_from(range(n_nodes))
for i in range(n_nodes):
    for j in range(n_nodes):
        w = plw_distances[i][j]
        if w != 0 and w != np.inf:
            G_plw.add_edge(i, j, weight=round(w,3))

fig, ax = plt.subplots(figsize=(10,10))

pos = nx.kamada_kawai_layout(G_plw, weight='weight') # pos = nx.spring_layout(G), etc.
labels = {i: f'$v_{{{i}}}$' for i in G_plw.nodes}
node_colors = 'lightblue' 

nx.draw(G_plw, pos, labels=labels, with_labels=True, ax=ax, node_color=node_colors, node_size=2000, font_weight='bold', font_size=20)
weights = {(u, v): d['weight'] for u, v, d in G_plw.edges(data=True)}
nx.draw_networkx_edge_labels(G_plw, pos, edge_labels=weights, font_size=16)
plt.title('path-length-weights distances of the graph')
plt.show()