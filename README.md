# Project Documentation

## Complie Instructions
Users can simply use `python analysis.py` to generate the graph of the city subway network under a certain attack.

Users can change the city by selecting the cities array index in our code.

```python
cities = ["M3_Shanghai.graphml", "M3_NYC.graphml", "M3_London.graphml", "M3_Chicago.graphml"]
G = nx.read_graphml(cities[0])
```

Users can change the property method by selecting the methods array index in our code.

```python
methods = ["efficiency", "betweenness_node", "betweenness_edge", "network_size", "LCC", "FLN"]
method = methods[0]
```

TFL-Tube-Data.json and Tube.py are the example network generation data of the city London. Users can use `python Tube.py` to generate the graphml file and import it into Gephi for further analysis.

## Source Files
In this project, we have the following python files, graphml files, and gephi files.

### TFL-Tube-Data.json (subway example)
This is the London Underground json file we fetch from Transportation For London website and [GitHub repo](https://gist.github.com/JhumanJ/59b785798248598a0e95c3e25eaad346). We use this file to generate the network of London underground.

### Tube.py (subway example)
This is the python script we use to generate the network file (graphml).

### analysis.py (main script)
This is the main python file to analyze the robutness property under five attacking propotols. Users can change the setting as discussed in Complie Instructions.

### M3_city.graphml
These four graphml files can be used both by Python and Gephi.

### M3_city.gephi
These four gephi files are the Gephi files for network property analyzation.

## Important Functions

### calculate_method(G)
calculate the property method value of the graph G.

### remove_nodes_protocol(G, p)
Remove the nodes of graph G based on the fraction of removed nodes p and the specific protocol. (Remove size(G.nodes) * p nodes from G)

### generate_method()
Generate the property method values of five attacking protocols (Random, Largest Degree, Highest Betweennes, Articulation Points, Collective Influence).

### plot_graph()
Plot the fraction of removed nodes vs. property method value graph based on setting.

## Implementation

Cities: Shanghai, New York City, London, Chicago.

Methods: Efficiency, Betweenness Node, Betweenness Edge, Network Size, Largest Clustering Component, Functionality Loss.

Protocols: Random, Largest Degree, Highest Betweennes, Articulation Points, Collective Influence.
