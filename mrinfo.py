#!/usr/bin/env python2.7
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

    TYPE_FLAT    = "flat"
    TYPE_PURE_L3 = "pure-l3"
    TYPE_BIP     = "bipartite"
    TYPE_PURE_BIP= "pure-bipartite"

    def __init__(self, filename):
        self.__vertices= {}
        self.__id2vertex= {}
        self.__edges= {}
        self.__next_vertex_id= 0
        self.__load(filename)
        self.__types= {MRINFO.TYPE_FLAT:self.__to_flat_graph,
                       MRINFO.TYPE_PURE_L3:self.__to_l3_graph,
                       MRINFO.TYPE_BIP:self.__to_bipartite,
                       MRINFO.TYPE_PURE_BIP:self.__to_pure_bipartite}
        return
    

    @classmethod
    def get_types(cls):
        return [MRINFO.TYPE_FLAT,
                MRINFO.TYPE_PURE_L3,
                MRINFO.TYPE_BIP,
                MRINFO.TYPE_PURE_BIP]

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
    def edge_count(self, typ=None):
        count= 0
        for head in self.__edges:
            count+= len(self.__edges[head])
        return count

    # -----[ to_graph ]----------------------------------------------
    def to_graph(self, typ):
        if not(typ in self.__types):
            return None
        return self.__types[typ]()

    # -----[ __to_flat_graph ]---------------------------------------
    # Return a flat graph with all vertices and edges
    def __to_flat_graph(self):
        edges= []
        for head in self.__edges:
            for tail in self.__edges[head]:
                edges.append((self.__vertices[head].get_id(),
                              self.__vertices[tail].get_id()))
        g= igraph.Graph(self.vertex_count(), edges)
        for v in g.vs:
            if self.__id2vertex[v.index].get_type() == MRINFO_VERTEX.L2:
                v["type"]= False
            else:
                v["type"]= True
        return g

    # -----[ __to_l3_graph ]-----------------------------------------
    # Return a flat graph with only the L3-L3 edges
    def __to_l3_graph(self):
        next_id= 0
        idmap= {}
        edges= []
        for i in range(self.vertex_count()):
            v= self.__get_vertex_by_id(i)
            if v.get_type() == MRINFO_VERTEX.L3:
                idx= v.get_id()
                if not(idx in idmap):
                    idmap[idx]= next_id
                    next_id+= 1
        for head in self.__edges:
            for tail in self.__edges[head]:
                head_n= self.__vertices[head]
                tail_n= self.__vertices[tail]
                if head_n.get_type() == MRINFO_VERTEX.L3 and \
                   tail_n.get_type() == MRINFO_VERTEX.L3:
                    edges.append((idmap[head_n.get_id()], \
                                  idmap[tail_n.get_id()]))
        g= igraph.Graph(len(idmap), edges)
        for v in g.vs:
            v["type"]= True
        return g

    # -----[ __to_biparite ]-----------------------------------------
    # Return a bipartite graph with L3-L3 edges converted to L3-L2-L3
    def __to_bipartite(self):
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

    # -----[ __to_pure_bipartite ]-----------------------------------
    def __to_pure_bipartite(self):
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
