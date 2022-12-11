#!/usr/bin/env python3

import networkx as nx

def calculate_average_betweenness_node(G):
	sum_betweenness = 0.0
	betweenness_dict = nx.betweenness_centrality(G, normalized=False)
	return sum(betweenness_dict.values()) / len(G.nodes())

G = nx.Graph()



G.add_node(0)
G.add_node(1)
G.add_node(2)
G.add_node(3)
G.add_node(4)



G.add_edge(0, 1)
G.add_edge(1, 2)
G.add_edge(1, 3)
G.add_edge(2, 3)
G.add_edge(2, 4)
G.add_edge(1, 4)

print(calculate_average_betweenness_node(G))


