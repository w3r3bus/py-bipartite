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


class MRINFO_VERTEX:

    L2= False
    L3= True

    def __init__(self, name, id, type):
        self.__name= name
        self.__id= id
        self.__type= type

    def get_id(self):
        return self.__id

    def get_name(self):
        return self.__name

    def get_type(self):
        return self.__type


class MRINFO:

    def __init__(self, filename):
        self.__vertices= {}
        self.__id2vertex= {}
        self.__edges= {}
        self.__next_vertex_id= 0
        self.__load(filename)

    # -----[ __add_vertex ]------------------------------------------
    # Add a vertex to the database if it does not already exist.
    # Generate a unique ID for the new vertex (or retrieve the ID if
    # vertex already exists)
    def __add_vertex(self, name, type):
        if name in self.__vertices:
            v= self.__vertices[name]
        else:
            v= MRINFO_VERTEX(name, self.__next_vertex_id, type)
            self.__vertices[name]= v
            self.__id2vertex[v.get_id()]= v
            self.__next_vertex_id+= 1
        return v.get_id()

    # -----[ __add_edge ]--------------------------------------------
    # Add an edge between two vertices. The head and tail vertices
    # must already exist.
    def __add_edge(self, head_name, tail_name):
        if not(head_name in self.__edges):
            self.__edges[head_name]= {}
        self.__edges[head_name][tail_name]= 1

    # -----[ __get_vertex_by_id ]------------------------------------
    # Obtain a vertex from its ID
    def __get_vertex_by_id(self, id):
        return self.__id2vertex[id]

    # -----[ __load ]------------------------------------------------
    # Load a file in MRINFO format
    def __load(self, filename):
        # Each line in an mrinfo file declares a single adjacency.
        # An adjacency is <SRC>: <DST>: <NUMBER>
        # if <SRC> or <DST> has format "ID[0-9]+\[IP\]" it is a switch (L2)
        # otherwise, it is a router
        # Header lines start with '#' and are ignored
        prog_switch= re.compile("ID([0-9]+).+")
        with open(filename, 'r') as f:
            for line in f:
                if line[0] == '#': continue
                fields= line.split(' ')
                (src, dst)= [name.strip(': ') for name in fields[0:2]]
                for name in [src, dst]:
                    # L2/L3 node ?
                    m= prog_switch.match(name)
                    # Register node, get node id
                    if (m == None):
                        type= MRINFO_VERTEX.L3
                    else:
                        type= MRINFO_VERTEX.L2
                    self.__add_vertex(name, type)
                self.__add_edge(src, dst)
        return

    # -----[ vertex_count ]------------------------------------------
    # Return number of vertices. If a type is provided, return the
    # number of vertices of that type.
    def vertex_count(self, typ=None):
        if typ != None:
            return len([n for n in self.__vertices \
                        if self.__vertices[n].get_type() == typ])            
        return len(self.__vertices)

    # -----[ edge_count ]--------------------------------------------
    # Return number of edges
    def edge_count(self):
        count= 0
        for head in self.__edges:
            count+= len(self.__edges[head])
        return count

    # -----[ to_graph ]----------------------------------------------
    # Return a flat graph with all vertices and edges
    def to_graph(self):
        edges= []
        for head in self.__edges:
            for tail in self.__edges[head]:
                edges.append((self.__vertices[head].get_id(),
                              self.__vertices[tail].get_id()))
        return igraph.Graph(self.vertex_count(), edges)

    # -----[ to_l3l3_graph ]-----------------------------------------
    # Return a flat graph with only the L3-L3 edges
    def to_l3l3_graph(self):
        edges= []
        for head in self.__edges:
            for tail in self.__edges[head]:
                head_n= self.__vertices[head]
                tail_n= self.__vertices[tail]
                if head_n.get_type() == MRINFO_VERTEX.L3 and \
                   tail_n.get_type() == MRINFO_VERTEX.L3:
                    edges.append((head_n.get_id(), tail_n.get_id()))
        return igraph.Graph(self.vertex_count(), edges)

    # -----[ to_biparite ]-------------------------------------------
    # Return a bipartite graph with L3-L3 edges converted to L3-L2-L3
    def to_bipartite(self):
        next_virtual_vertex_id= self.__next_vertex_id
        types= []
        for i in range(self.vertex_count()):
            typ= self.__get_vertex_by_id(i).get_type()
            if typ == MRINFO_VERTEX.L2:
                types.append(False)
            else:
                types.append(True)
        edges= []
        for head in self.__edges:
            for tail in self.__edges[head]:
                head_n= self.__vertices[head]
                tail_n= self.__vertices[tail]
                if head_n.get_type() == MRINFO_VERTEX.L3 and \
                   tail_n.get_type() == MRINFO_VERTEX.L3:
                    virtual_l2_id= next_virtual_vertex_id
                    types.append(False)
                    edges.append((head_n.get_id(), virtual_l2_id))
                    edges.append((virtual_l2_id, tail_n.get_id()))
                    next_virtual_vertex_id+= 1
                else:
                    edges.append((head_n.get_id(), tail_n.get_id()))
        return igraph.Graph.Bipartite(types, edges)

    # -----[ to_pure_bipartite ]-------------------------------------
    def to_pure_bipartite(self):
        types= []
        for i in range(self.vertex_count()):
            typ= self.__get_vertex_by_id(i).get_type()
            if typ == MRINFO_VERTEX.L2:
                types.append(False)
            else:
                types.append(True)
        edges= []
        for head in self.__edges:
            for tail in self.__edges[head]:
                head_n= self.__vertices[head]
                tail_n= self.__vertices[tail]
                if head_n.get_type() == MRINFO_VERTEX.L3 and \
                   tail_n.get_type() == MRINFO_VERTEX.L3:
                    continue
                edges.append((head_n.get_id(), tail_n.get_id()))
        return igraph.Graph.Bipartite(types, edges)

# -----[ test1 ]-----------------------------------------------------
def test1(filename):
    g= MRINFO.load(filename)

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
    plt.suptitle(filename)
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

# -----[ test2 ]-----------------------------------------------------
def test2(filename):
    mi= MRINFO(filename)
    print "Number of nodes:", mi.vertex_count()
    print "Number of L2 nodes: ", mi.vertex_count(MRINFO_VERTEX.L2)
    print "Number of L3 nodes: ", mi.vertex_count(MRINFO_VERTEX.L3)
    print "Number of edges: ", mi.edge_count()
    print "> complete flat graph"
    g= mi.to_graph()
    print g
    print "> pure L3/L3 flat graph"
    g= mi.to_l3l3_graph()
    print g
    print "> bipartite graph"
    g= mi.to_bipartite()
    print g
    print "> pure bipartite graph"
    g= mi.to_pure_bipartite()
    print g
    return

# -----[ main ]------------------------------------------------------
#
# -------------------------------------------------------------------
def main():
    if len(sys.argv) != 2:
        usage()
        sys.exit(-1)
    #test1(sys.argv[1])
    test2(sys.argv[1])
    return


if __name__ == "__main__":
    import igraph
    import sys
    import numpy
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    from pareto import Pareto
    main()
