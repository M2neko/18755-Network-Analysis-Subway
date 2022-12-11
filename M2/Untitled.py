#!/usr/bin/env python3

import networkx as nx


def calculate_average_betweenness_node(G):
	sum_betweenness = 0.0
	betweenness_dict = nx.betweenness_centrality(G, normalized=False)
	return sum(betweenness_dict.values()) / len(G.nodes())


G_1 = nx.read_gml("M2_London.gml")
G = nx.read_graphml("M2_London_re.graphml")
print(calculate_average_betweenness_node(G_1))
print(calculate_average_betweenness_node(G))