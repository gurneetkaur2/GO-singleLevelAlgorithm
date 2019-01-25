import networkx as nx
import matplotlib.pyplot as plt

g=nx.read_edgelist('test', create_using=nx.Graph(),nodetype=int)

print (nx.info(g))

nx.draw(g)

plt.show()
