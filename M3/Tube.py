import json
import networkx as nx

G = nx.MultiGraph()

f = open('TFL-Tube-Data.json')

data = json.load(f)

f.close()

# add all station nodes to the graph
for station in data['stations']:
	G.add_node(station['tla'])
	
# add edge from the previous station to this station to the graph
for line in data['lines']:
	for i in range(len(line['stations']) - 1):
		G.add_edge(line['stations'][i], line['stations'][i + 1])

G.remove_nodes_from(list(nx.isolates(G)))
	
nx.write_graphml(G, "M3_London.graphml")