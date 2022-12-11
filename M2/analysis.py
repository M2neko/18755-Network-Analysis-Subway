import json
import copy
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt



def calculate_efficiency(G):
	sum_efficiency = 0.0
	for i in G.nodes():
		for j in G.nodes():
			if i != j:
				sum_efficiency += nx.efficiency(G, i, j)
				
	sum_efficiency /= 2
	
	N = len(G.nodes())
	
	sum_efficiency /= N * (N - 1)
	
	return sum_efficiency

def calculate_average_betweenness_node(G):
	betweenness_dict = nx.betweenness_centrality(G, normalized=True)
	return sum(betweenness_dict.values()) / len(G.nodes())

def calculate_average_betweenness_edge(G):
	edge_betweenness_centrality_dict = nx.edge_betweenness_centrality(G, normalized=True)
	return sum(edge_betweenness_centrality_dict.values()) / len(edge_betweenness_centrality_dict)
	
def calculate_network_size(G):
	return G.size()

def calculate_LCC(G, s0):
	s = -1
	for connected_components in nx.connected_components(G):
		s = max(s, len(connected_components))
	return s / s0

# section 4
def remove_nodes_random_FLN(G, p):
	FLN = 0
	f_map = {}
	fl_map = {}
	for node in G.nodes():
		f_map[node] = 1
		fl_map[node] = 0

	N = len(G.nodes())
	num_removed = int(N * p)
	removed_list = np.random.choice(G.nodes(), num_removed, replace=False)
	for i in removed_list:
		for j in G.nodes():
			if i == j or not G.has_edge(i, j):
				continue
			original_fl = f_map[j]
			new_fl = f_map[j] - f_map[j] / (nx.shortest_path_length(G, source=i, target=j) ** 2 * G.degree[j])
			f_map[j] = new_fl
			fl_map[j] += original_fl - new_fl

		del f_map[i]
		del fl_map[i]
		G.remove_node(i)
	
	for _, val in fl_map.items():
		FLN += val	

	return FLN

def remove_nodes_degree_FLN(G, p):
	FLN = 0
	f_map = {}
	fl_map = {}
	for node in G.nodes():
		f_map[node] = 1
		fl_map[node] = 0
	N = len(G.nodes())
	num_removed = int(N * p)
	degree_list = sorted(G.degree, key=lambda x: x[1], reverse=True)
	removed_list = degree_list[0 : num_removed]
	for [node, val] in removed_list:
		i = node
		for j in G.nodes():
			if i == j or not G.has_edge(i, j):
				continue
			original_fl = f_map[j]
			new_fl = f_map[j] - f_map[j] / (nx.shortest_path_length(G, source=i, target=j) ** 2 * G.degree[j])
			f_map[j] = new_fl
			fl_map[j] += original_fl - new_fl

		del f_map[i]
		del fl_map[i]
		G.remove_node(i)
	
	for _, val in fl_map.items():
		FLN += val	
		
	return FLN


def remove_nodes_betweenness_FLN(G, p):
	FLN = 0
	f_map = {}
	fl_map = {}
	for node in G.nodes():
		f_map[node] = 1
		fl_map[node] = 0

	N = len(G.nodes())
	num_removed = int(N * p)

	betweenness_dict = nx.betweenness_centrality(G)
	betweenness_list = sorted(betweenness_dict.items(), key=lambda x:x[1], reverse=True)
	removed_list = betweenness_list[0 : num_removed]
	for [node, val] in removed_list:
		i = node
		for j in G.nodes():
			if i == j or not G.has_edge(i, j):
				continue
			original_fl = f_map[j]
			new_fl = f_map[j] - f_map[j] / (nx.shortest_path_length(G, source=i, target=j) ** 2 * G.degree[j])
			f_map[j] = new_fl
			fl_map[j] += original_fl - new_fl

		del f_map[i]
		del fl_map[i]
		G.remove_node(i)
	
	for _, val in fl_map.items():
		FLN += val	
	return FLN


#############


def remove_nodes_random(G, r):
	removed_list = np.random.choice(G.nodes(), r, replace=False)
	for i in removed_list:
		G.remove_node(i)
	return G

def remove_nodes_degree(G, p):
	N = len(G.nodes())
	num_removed = int(N * p)
	degree_list = sorted(G.degree, key=lambda x: x[1], reverse=True)
	removed_list = degree_list[0 : num_removed]
	for [node, val] in removed_list:
		G.remove_node(node)
	return G

def remove_nodes_betweenness(G, p):
	N = len(G.nodes())
	num_removed = int(N * p)
	betweenness_dict = nx.betweenness_centrality(G)
	betweenness_list = sorted(betweenness_dict.items(), key=lambda x:x[1], reverse=True)
	removed_list = betweenness_list[0 : num_removed]
	for [node, val] in removed_list:
		G.remove_node(node)
	return G


def remove_nodes_articulation(G, r):
	
	removed = list(nx.articulation_points(G))
	
	if len(removed) >= r:
		removed_list = np.random.choice(removed, r, replace=False)
		for i in removed_list:
			G.remove_node(i)
	else:
		remaining = r - len(removed)
		for i in removed:
			G.remove_node(i)
		
		random_removed_list = np.random.choice(G.nodes(), remaining, replace=False)
		
		for i in random_removed_list:
			G.remove_node(i)
		
	return G


def remove_nodes_collective(G, r):
	
	l = 3
	
	CI = {}
	for i in G.nodes():
		CI[i] = 0
	
	paths = dict(nx.all_pairs_shortest_path(G))

	for i in paths.keys():
		for j in paths[i].keys():
			if i != j and len(paths[i][j]) <= l:
				CI[i] += (G.degree[i] - 1) * (G.degree[j] - 1)
	
	CI_list = sorted(CI.items(), key=lambda x:x[1], reverse=True)
	removed_list = CI_list[0 : r]
	for [node, val] in removed_list:
		G.remove_node(node)
		
	return G
	





##############################################

def generate_removed_graph(G, p, G_random, G_articulation, G_colletive, r):
	G_removed_betweenness = remove_nodes_betweenness(copy.deepcopy(G), p)
	G_removed_degree = remove_nodes_degree(copy.deepcopy(G), p)
	G_removed_random = remove_nodes_random(G_random, r)
	G_removed_articulation = remove_nodes_articulation(G_articulation, r)
	G_removed_collective = remove_nodes_collective(G_colletive, r)
	return [G_removed_betweenness, G_removed_degree, G_removed_random, G_removed_articulation, G_removed_collective]


def generate_efficiency(G, p, G_random, G_articulation, G_collective,r):
	
	[G_removed_betweenness, G_removed_degree, G_removed_random, G_removed_articulation, G_removed_collective] = generate_removed_graph(G, p, G_random, G_articulation, G_collective, r)
	
	e_betweenness = calculate_efficiency(G_removed_betweenness)
	e_degree = calculate_efficiency(G_removed_degree)
	e_random = calculate_efficiency(G_removed_random)
	e_articulation = calculate_efficiency(G_removed_articulation)
	e_collective = calculate_efficiency(G_removed_collective)
	return [e_betweenness, e_degree, e_random, e_articulation, e_collective]



def generate_betweenness_node(G, p, G_random, G_articulation, G_collective, r):
	
	[G_removed_betweenness, G_removed_degree, G_removed_random, G_removed_articulation, G_removed_collective] = generate_removed_graph(G, p, G_random, G_articulation, G_collective, r)
	
	
	e_betweenness = calculate_average_betweenness_node(G_removed_betweenness)
	e_degree = calculate_average_betweenness_node(G_removed_degree)
	e_random = calculate_average_betweenness_node(G_removed_random)
	e_articulation = calculate_average_betweenness_node(G_removed_articulation)
	e_collective = calculate_average_betweenness_node(G_removed_collective)
	return [e_betweenness, e_degree, e_random, e_articulation, e_collective]


def generate_betweenness_edge(G, p, G_random, G_articulation, G_collective, r):
	
	[G_removed_betweenness, G_removed_degree, G_removed_random, G_removed_articulation, G_removed_collective] = generate_removed_graph(G, p, G_random, G_articulation, G_collective, r)
	
	
	e_betweenness = calculate_average_betweenness_edge(G_removed_betweenness)
	e_degree = calculate_average_betweenness_edge(G_removed_degree)
	e_random = calculate_average_betweenness_edge(G_removed_random)
	e_articulation = calculate_average_betweenness_edge(G_removed_articulation)
	e_collective = calculate_average_betweenness_node(G_removed_collective)
	return [e_betweenness, e_degree, e_random, e_articulation, e_collective]


def generate_network_size(G, p, network_size, G_random, G_articulation, G_collective, r):
	
	[G_removed_betweenness, G_removed_degree, G_removed_random, G_removed_articulation, G_removed_collective] = generate_removed_graph(G, p, G_random, G_articulation, G_collective, r)
	
	
	e_betweenness = calculate_network_size(G_removed_betweenness) / network_size
	e_degree = calculate_network_size(G_removed_degree) / network_size
	e_random = calculate_network_size(G_removed_random) / network_size
	e_articulation = calculate_network_size(G_removed_articulation) / network_size
	e_collective = calculate_network_size(G_removed_collective) / network_size
	
	
	return [e_betweenness, e_degree, e_random, e_articulation, e_collective]

def generate_LCC(G, p, s0, G_random, G_articulation, G_collective, r):
	
	[G_removed_betweenness, G_removed_degree, G_removed_random, G_removed_articulation, G_removed_collective] = generate_removed_graph(G, p, G_random, G_articulation, G_collective, r)
	
	
	e_betweenness = calculate_LCC(G_removed_betweenness, s0)
	e_degree = calculate_LCC(G_removed_degree, s0)
	e_random = calculate_LCC(G_removed_random, s0)
	e_articulation = calculate_LCC(G_removed_articulation, s0)
	e_collective = calculate_LCC(G_removed_collective, s0)
	return [e_betweenness, e_degree, e_random, e_articulation, e_collective]

def generate_FLN(G, p):
	
	e_betweenness = remove_nodes_betweenness_FLN(copy.deepcopy(G), p)
	e_degree = remove_nodes_degree_FLN(copy.deepcopy(G), p)
	e_random = remove_nodes_random_FLN(copy.deepcopy(G), p)
	return [e_betweenness, e_degree, e_random]
	

def generate(G, p, method, G_random, G_articulation, G_collective):
	
	if p == 0.00:
		r = 0
	else:
		r = int(len(G.nodes()) * STEP)
	
	if method == "efficiency":
		[e_betweenness, e_degree, e_random, e_articulation, e_collective] = generate_efficiency(G, p, G_random, G_articulation, G_collective, r)

	elif method == "betweenness_node":
		[e_betweenness, e_degree, e_random, e_articulation, e_collective] = generate_betweenness_node(G, p, G_random, G_articulation, G_collective, r)
		
	elif method == "betweenness_edge":
		[e_betweenness, e_degree, e_random, e_articulation, e_collective] = generate_betweenness_edge(G, p, G_random, G_articulation, G_collective, r)
		
	elif method == "network_size":
		network_size = calculate_network_size(G)
		[e_betweenness, e_degree, e_random, e_articulation, e_collective] = generate_network_size(G, p, network_size, G_random, G_articulation, G_collective, r)
		
	elif method == "LCC":
		s0 = len(G.nodes())
		[e_betweenness, e_degree, e_random, e_articulation, e_collective] = generate_LCC(G, p, s0, G_random, G_articulation, G_collective, r)
		
	elif method == "FLN":
		[e_betweenness, e_degree, e_random] = generate_FLN(G, p)
		e_articulation = 0
		e_collective = 0
		
	return [e_betweenness, e_degree, e_random, e_articulation, e_collective]

def plot_graph(p_list, e_list_random, e_list_degree, e_list_betweenness, e_list_articulation, e_list_collective, method):
	plt.plot(p_list, e_list_random, label = "Random", c = 'b', marker='^', ls='--')
	plt.plot(p_list, e_list_degree, label = "Largest Degree", c = 'r', marker = 'o', ls='-')
	plt.plot(p_list, e_list_betweenness, label = "Highest Betweenness", c = 'g', marker = 'v', ls=':')
	plt.plot(p_list, e_list_articulation, label = "Articulation Targeted", c = 'y', marker = 's', ls=':')
	plt.plot(p_list, e_list_collective, label = "Collective Influence", c = 'c', marker = 'p', ls='--')
	
	plt.title(method)
	
	plt.xlabel("Fraction of removed nodes")
	plt.ylabel("Value")
	plt.legend()
	plt.show()
#	plt.savefig("graph.png")


cities = ["M2_Shanghai.graphml", "M2_NYC.graphml", "M2_London.graphml", "M2_Chicago.graphml"]

G = nx.read_graphml(cities[3])

MAX_P = 0.50
STEP = 0.01
p = 0.00

p_list = []
e_list_random = []
e_list_degree = []
e_list_betweenness = []
e_list_articulation = []
e_list_collective = []

G_random = copy.deepcopy(G)
G_articulation = copy.deepcopy(G)
G_collective = copy.deepcopy(G)


methods = ["efficiency", "betweenness_node", "betweenness_edge", "network_size", "LCC", "FLN"]

method = methods[4]

while p <= MAX_P:
	[e_betweenness, e_degree, e_random, e_articulation, e_collective] = generate(G, p, method, G_random, G_articulation, G_collective)
	e_list_betweenness.append(e_betweenness)
	e_list_degree.append(e_degree)
	e_list_random.append(e_random)
	e_list_articulation.append(e_articulation)
	e_list_collective.append(e_collective)
	
	p_list.append(p)
	print("Complete: p =", f'{p:.{2}f}')
	p += STEP
	
plot_graph(p_list, e_list_random, e_list_degree, e_list_betweenness, e_list_articulation, e_list_collective, method)
