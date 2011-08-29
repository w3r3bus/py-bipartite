# ===================================================================
# Import mrinfo graph file
#
# (C) 2011, Bruno Quoitin (bruno.quoitin@umons.ac.be)
# ===================================================================


import re
import igraph


# -----[ usage ]-----------------------------------------------------
#
# -------------------------------------------------------------------
def usage():
    print "Usage: python mrinfo.py FILENAME"
    print


class MRINFO:

    # -----[ load ]--------------------------------------------------
    @classmethod
    def load(cls, filename):
        nodes= {}
        next_node_id= 0
        types= []

        edges= []
        num_ignored= 0

        prog_switch= re.compile("ID([0-9]+).+")

        # Each line in an mrinfo file declares a single adjacency.
        # An adjacency is <SRC>: <DST>: <NUMBER>
        # if <SRC> or <DST> has format "ID[0-9]+\[IP\]" it is a switch (L2)
        # otherwise, it is a router
        # Header lines start with '#' and are ignored
        with open(filename, 'r') as f:
            for line in f:
                if line[0] == '#': continue
                fields= line.split(':')
                (src, dst)= [name.strip() for name in fields[0:2]]
                for name in [src, dst]:
                    # L2/L3 node ?
                    m= prog_switch.match(name)
                    # Register node, get node id
                    if not(name in nodes):
                        # node: (id, type=False/True for L2/L3)
                        nodes[name]= next_node_id
                        types.append(m == None)
                        next_node_id+= 1
                # L2/L3 edges are added as is.
                if (types[nodes[src]] != types[nodes[dst]]):
                    edges.append((nodes[src], nodes[dst]))
                # L3/L3 edges are converted to a new L2 vertex +
                # two L2/L3 edges.
                else:
                    types.append(False)
                    edges.append((nodes[src], next_node_id))
                    edges.append((next_node_id, nodes[dst]))
                    next_node_id+= 1

        # To create an igraph bipartite graph we need to provide the
        # types of vertices (False for L2 nodes, True for L3 nodes)
        # and the edges
        g= igraph.Graph.Bipartite(types, edges)

        return g


# -----[ main ]------------------------------------------------------
#
# -------------------------------------------------------------------
def main():
    if len(sys.argv) != 2:
        usage()
        sys.exit(-1)

    g= MRINFO.load(sys.argv[1])

    l2_vertices= [idx for idx in range(g.vcount()) if not(g.vs["type"][idx])]
    l3_vertices= [idx for idx in range(g.vcount()) if g.vs["type"][idx]]
    
    print "Number of vertices: %d" % (g.vcount())
    print "Number of L2 vertices: %d" % (len(l2_vertices))
    print "Number of L3 vertices: %d" % (len(l3_vertices))
    print "Number of edges: %d" % (g.ecount())
    clusters= g.clusters()
    print "Number of connected components: %d" % (len(clusters))
    print "Size of largest connected component: %d" % (clusters.giant().vcount())

    degrees= g.degree()
    l2_degrees= g.degree(l2_vertices)
    l3_degrees= g.degree(l3_vertices)
    max_degree= g.maxdegree()
    print "Highest degree: %d" % (max_degree)
    degrees_hist= numpy.histogram(degrees, max_degree)
    l2_degrees_hist= numpy.histogram(l2_degrees, max_degree)
    l3_degrees_hist= numpy.histogram(l3_degrees, max_degree)


    print "MLE (overall): %.20f" % (Pareto.mle(degrees))
    print "MLE (L2): %.20f" % (Pareto.mle(l2_degrees))
    print "MLE (L2): %.20f" % (Pareto.mle(l2_degrees, 2.0))
    print "MLE (L3): %.20f" % (Pareto.mle(l3_degrees))


    plt.figure()
    plt.suptitle(sys.argv[1])
    plt.subplot(3, 2, 1)
    plt.title("Overall degree distribution", fontsize=10)
    plt.hist(degrees, bins=max_degree, range=(0, max_degree), cumulative=False, normed=False, histtype='bar')
    plt.subplot(3, 2, 2)
    plt.title("(log-log)", fontsize=10)
    plt.grid()
    plt.loglog(range(0, max_degree), degrees_hist[0], 'o')
    plt.subplot(3, 2, 3)
    plt.title("L2 degree distribution", fontsize=10)
    plt.hist(l2_degrees, bins=max_degree, range=(0, max_degree), cumulative=False, normed=False, histtype='bar')
    plt.subplot(3, 2, 4)
    plt.title("(log-log)", fontsize=10)
    plt.grid()
    plt.loglog(range(0, max_degree), l2_degrees_hist[0], 'o')
    plt.subplot(3, 2, 5)
    plt.title("L3 degree distribution", fontsize=10)
    plt.hist(l3_degrees, bins=max_degree, range=(0, max_degree), cumulative=False, normed=False, histtype='bar')
    plt.subplot(3, 2, 6)
    plt.title("(log-log)", fontsize=10)
    plt.grid()
    plt.loglog(range(0, max_degree), l3_degrees_hist[0], 'o')
    plt.show()
    
    #l2_degrees= g.degree([id for id where g.types[id] == False])

    return


if __name__ == "__main__":
    from igraph import igraph
    import sys
    import numpy
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    from pareto import Pareto
    main()
