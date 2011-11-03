#!/usr/bin/env python2.7
# ===================================================================
# Random bipartite graph generation with random or prescribed degrees
#
# (C) 2011, Bruno Quoitin (bruno.quoitin@umons.ac.be)
# ===================================================================
#
# This program is organized as follows:
#
#   The main entry is in function 'main' which is responsible for
#   parsing the command line.
#
#   The user provides a variable number of operations along with
#   their arguments. More details about the supported operations
#   and their arguments can be obtained by running the program
#   without argument (see function 'usage'). All operations are
#   defined in the data structure 'operations').
#
#   The most important operations are named 'genpareto' and 'genzipf'
#   (see functions 'bipartite_op_gen_pareto' and 'bipartite_op_gen_zipf').
#   These operations generate a bipartite graph using L2/L3 degree sequences
#   obtained from 2 Pareto/Zipf distributions. These operations require
#   the user to provide the number of L2/L3 vertices and the parameters of
#   the power-law distribution (such as the exponent).
#
#     python bipartite.py genzipf:2.5:20:3:-1:1.5:10:2:100 project plot
#
#   The above invocation specifies that the generator must produce
#   a bipartite graph with 100 routers. The L2 degrees distribution
#   must follow a Zipf with exponent 2.5 and only contain degrees in
#   interval [3,20]. The L3 degrees distribution must follow a Zipf
#   distribution with exponent 1.5 and only contain degrees in interval
#   [2,10]. The number of L2 devices (switches) is left undefined (-1)
#   to allow the generator to compute it so as to satisfy the following
#   basic property of bipartite graphs: the sums of L2 degrees and L3
#   degrees must match.
#
#   The bipartite graph is then L3-projected and displayed on the screen.
#
#   Other important operations are
#     'giant' which extracts the largest connected component (see
#       function 'bipartite_op_giant'),
#     'save' which allows to save the generated graph into a file
#       (see function 'bipartite_op_save')
#     'project' which projects the bipartite graph onto the L3 level
#       ( see function 'bipartite_op_project').
#
# ===================================================================


import igraph
import igraph.drawing
import numpy
import random
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

#from pareto import Pareto
from mrinfo import MRINFO
from zipf import Zipf


# -----[ Constants ]-------------------------------------------------

TYPE_L2= False
TYPE_L3= True

COLOR_DICT= {TYPE_L2:"blue", TYPE_L3:"green"}
SHAPE_DICT= {TYPE_L2:"rectangle", TYPE_L3:"circle"}
STR_DICT= {TYPE_L2:"L2", TYPE_L3:"L3"}
FONT_SIZE= 10


# -----[ Global options ]--------------------------------------------
verbosity = 0
max_redraws = 1000000


# -----[ log ]-------------------------------------------------------
# Display log message if current verbosity level is higher or equel
# than level in argument.
# -------------------------------------------------------------------
def log(level, msg):
    if (verbosity >= level):
        print msg
    return


# -----[ error ]-----------------------------------------------------
# Display error message and exit with error
# -------------------------------------------------------------------
def error(msg):
    print >> sys.stderr, "Error: %s" % (msg)
    sys.exit(-1)
    return


# -----[ info ]------------------------------------------------------
# Display information message (on stderr)
# -------------------------------------------------------------------
def info(msg):
    print >> sys.stderr, msg
    return


# -----[ cm ]--------------------------------------------------------
# Flat configuration model (Molloy-Reed)
# -------------------------------------------------------------------
def cm(degrees):
    # Generate stubs
    stubs= []
    for i in range(len(degrees)):
        for j in range(degrees[i]):
            stubs.append(i)

    # Connect stubs uniformly randomly
    edges= []
    while (len(stubs) >= 2):
        i= random.randint(0, len(stubs)-1)
        j= random.randint(0, len(stubs)-1)
        if (i != j):
            if (i > j):
                tmp= i
                i= j
                j= tmp
            edges.append((stubs[i], stubs[j]))
            stubs= stubs[0:i] + stubs[i+1:j] + stubs[j+1:]
    
    g= igraph.Graph(len(degrees), edges)
    return g


# -----[ bipartite_cm ]----------------------------------------------
# Bipartite configuration model. Generate a bipartite graph using the
# provided L2 and L3 degree sequences. The sums of L2 and L3 degrees
# MUST match.
# -------------------------------------------------------------------
def bipartite_cm(degrees2, degrees3):
    vertices_types= []
    edges= []

    td2= 0
    td3= 0

    stubs2= []
    stubs3= []
    id= 0

    # Generate L2/L3 stubs according to the respective prescribed
    # degree sequences
    for i in range(len(degrees2)):
        vertices_types.append(TYPE_L2)
        td2+= degrees2[i]
        for j in range(degrees2[i]):
            stubs2.append((id, j))
        id+= 1
        
    for i in range(len(degrees3)):
        vertices_types.append(TYPE_L3)
        td3+= degrees3[i]
        for j in range(degrees3[i]):
            stubs3.append((id, j))
        id+= 1

    # Check that total L2/L3 degrees match
    if (td2 != td3):
        error("Sum of L2 degrees does not match sum of L3 degrees")

    # Shuffle L3 stubs
    random.shuffle(stubs3)
    log(1, "lengths: %d ; %d" % (len(stubs2), len(stubs3)))

    # Connect L2 and L3 stubs
    for i in range(len(stubs2)):
        log(2, "connect L2/L3 stubs (%d, %d)" % (stubs2[i][0], stubs3[i][0]))
        edges.append((stubs2[i][0], stubs3[i][0]))

    # Generate Bipartite graph
    g= igraph.Graph.Bipartite(vertices_types, edges)

    return g


# -----[ random_matching_degrees ]-----------------------------------
# Generate matching L2/L3 degree sequences using provided
# probability distributions (e.g. Pareto, Zipf).
# -------------------------------------------------------------------
def random_matching_degrees(l2_dist, n2, l3_dist, n3):
    degrees2= []
    degrees3= []
    td2= 0
    td3= 0

    # Generate L2 degree sequence
    for i in range(n2):
        r= int(round(l2_dist.random()))
        degrees2.append(r)
        td2+= r

    # Generate L3 degree sequence
    for i in range(n3):
        r= int(round(l3_dist.random()))
        degrees3.append(r)
        td3+= r

    log(1, "# Total degrees: td2=%d  vs  td3=%d" % (td2, td3))

    # Try to ensure that L2/L3 degrees sums match by re-drawing the
    # degree for randomly picked pairs of L2/L3 nodes
    #
    # Note: another strategy might be to try to add nodes for the
    #       smallest degree distribution until a match is found.
    #       This part needs to be investigated further...
    num_redraws= 0
    while (td2 != td3):
        i= random.randint(0, n2-1)
        j= random.randint(0, n3-1)
        td2-= degrees2[i]
        td3-= degrees3[j]
        degrees2[i]= int(round(l2_dist.random()))
        degrees3[j]= int(round(l3_dist.random()))
        td2+= degrees2[i]
        td3+= degrees3[j]
        log(2, "#i=%d,j=%d   td2=%d  vs  td3=%d" % (i, j, td2, td3))
        num_redraws+= 1
        if (num_redraws > max_redraws):
            error("maximum redraws limit was reached")
            return None

    info("  Number of redraws required = %d" % (num_redraws))
    return (degrees2, degrees3)


# -------------------------------------------------------------------
# This class holds the paramaters required to generate a random
# degree sequence, based on a Pareto distribution. This class allows
# to automatically compute one of the parameters if it was left
# undefined.
# -------------------------------------------------------------------
class DegreeSeqParams:

    def __init__(self, n, alpha, beta):
        self.__n= n
        self.__alpha= alpha
        self.__beta= beta
        return

    def get_alpha(self):
        return self.__alpha

    def get_beta(self):
        return self.__beta

    def get_n(self):
        return self.__n

    def num_undefined(self):
        n= 0
        if (self.__n < 0): n+= 1
        if (self.__alpha < 0): n+= 1
        if (self.__beta < 0): n+= 1
        return n

    def solve_undefined(self, other):
        if (self.__n < 0):
            alpha2= other.get_alpha()
            beta2= other.get_beta()
            n2= other.get_n()
            self.__n= int((1.0 * alpha2) / self.__alpha \
                          * (1.0 * beta2) / self.__beta \
                          * n2 \
                          * (self.__alpha - 1.0) / (alpha2 - 1.0))
                      
        elif (self.__alpha < 0):
            n2= other.get_n()
            alpha2= other.get_alpha()
            beta2= other.get_beta()
            self.__alpha= (1.0 * alpha2) \
                          /(alpha2 - ((1.0 * self.__n) / n2) \
                            * (self.__beta/beta2) * (alpha2 - 1.0))

        elif (self.__beta < 0):
            alpha2= other.get_alpha()
            beta2= other.get_beta()
            n2= other.get_n()
            self.__beta= (1.0 * alpha2) / self.__alpha \
                         * (1.0 * n2) / self.__n \
                         * beta2 \
                         * (self.__alpha - 1.0) / (alpha2 - 1.0)

        return

    def exp_tot_degree(self):
        return self.__n * (self.__alpha * self.__beta)/(self.__alpha - 1)

    def dump(self):
        info("    number of nodes (N) = %d" % (self.__n))
        if (self.__alpha > 0):
            info("    exponent (A)        = %.20f" % (self.__alpha))
        else:
            info("print     exponent (A)        = undefined")
        if (self.__beta > 0):
            info("    scale (B)           = %.20f" % (self.__beta))
        else:
            info("    scale (B)           = undefined")
        return


#####################################################################
#
# Implementation of the graph operations/transformations
#
#####################################################################


# -----[ bipartite_op_gen_pareto ]-----------------------------------
# Generate a bipartite graph using random degrees drawn from
# two Pareto distributions.
#
# Note: the degrees are obtained by rounding the pseudo-random Pareto
# distributed numbers. This introduces a bias and this bias grows
# with the power-law exponent => prefer the Zipf distribution (see
# function 'bipartite_op_gen_zipf').
# -------------------------------------------------------------------
def bipartite_op_gen_pareto(g, params):
    if len(params) != 6:
        error("invalid number of arguments for \"genpareto\"")

    # Obtain degree distribution parameters
    l2_params= DegreeSeqParams(int(params[0]), float(params[1]), float(params[2]))
    l3_params= DegreeSeqParams(int(params[3]), float(params[4]), float(params[5]))

    # Check there is at most one undefined parameter
    if (l2_params.num_undefined() + l3_params.num_undefined() > 1):
        error("too many (%d) undefined parameters (1 max)" % (num_undefined))

    # Determine L2/L3 exponent automatically from other parameters
    l3_params.solve_undefined(l2_params)
    l2_params.solve_undefined(l3_params)

    info("  L2 Parameters:")
    l2_params.dump()
    info("  L3 Parameters:")
    l3_params.dump()

    # Check expected number of edges at level 2
    # against expected number of edges at level 3
    info("  Expected total L2 degree = %.20f" % (l2_params.exp_tot_degree()))
    info("  Expected total L3 degree = %.20f" % (l3_params.exp_tot_degree()))

    l2_dist= Pareto2(l2_params.get_alpha(), l2_params.get_beta())
    l3_dist= Pareto2(l3_params.get_alpha(), l3_params.get_beta())
    (degrees2, degrees3)= \
        random_matching_degrees(l2_dist, l2_params.get_n(), \
                                l3_dist, l3_params.get_n())

    return bipartite_cm(degrees2, degrees3)


# -----[ bipartite_op_gen_zipf ]-------------------------------------
# Generate a bipartite graph using random degrees drawn from
# two Zipf distributions.
# -------------------------------------------------------------------
def bipartite_op_gen_zipf(g, params):
    if len(params) != 8:
        error("invalid number of arguments for \"genzipf\"")

    zipf_alpha2= float(params[0])
    zipf_N2= int(params[1])
    zipf_xmin2= int(params[2])
    n2= int(params[3])
    zipf_alpha3= float(params[4])
    zipf_N3= int(params[5])
    zipf_xmin3= int(params[6])    
    n3= int(params[7])

    z2= Zipf(zipf_alpha2, zipf_N2, zipf_xmin2)
    z3= Zipf(zipf_alpha3, zipf_N3, zipf_xmin3)
    
    z2_exp= z2.expectation()
    z3_exp= z3.expectation()
    info("  L2 distrib expectation = %.20f" % (z2_exp))
    info("  L3 distrib expectation = %.20f" % (z3_exp))

    if (n2 <= 0) and (n3 > 0):
        n2= int(round((n3*z3_exp)/z2_exp))
        info("  n2 auto-computed = %d" % (n2))
    elif (n2 > 0) and (n3 <= 0):
        n3= int(round((n2*z2_exp)/z3_exp))
        info("  n3 auto-computed = %d" % (n3))
    else:
        error("  n2 and n3 cannot both be undefined");

    info("  a2=%f, N2=%d, xmin2=%d, n2=%d" % \
         (zipf_alpha2, zipf_N2, zipf_xmin2, n2))
    info("  a3=%f, N3=%d, xmin3=%d, n3=%d" % \
         (zipf_alpha3, zipf_N3, zipf_xmin3, n3))

    # Check expected number of edges at level 2
    # against expected number of edges at level 3
    l2_exp= z2_exp * n2
    l3_exp= z3_exp * n3
    info("  Expected total L2 degree = %.20f" % (l2_exp))
    info("  Expected total L3 degree = %.20f" % (l3_exp))

    (degrees2, degrees3)= random_matching_degrees(z2, n2, z3, n3)

    log(1, "  MLE L2 = %f" % (Zipf.mle(degrees2, zipf_xmin2)))
    log(1, "  MLE L3 = %f" % (Zipf.mle(degrees3, zipf_xmin3)))

    return bipartite_cm(degrees2, degrees3)


# -----[ bipartite_op_prescribed ]-----------------------------------
# Generate a bipartite graph using prescribed degree sequences from
# either a previously loaded bipartite graph or obtained from files
# -------------------------------------------------------------------
def bipartite_op_prescribed(g, params):
    if len(params) == 2:
        l2_filename= params[0]
        l3_filename=params[1]
        info("  Loading L2 degrees from \"%s\"" % (l2_filename))
        info("  Loading L3 degrees from \"%s\"" % (l3_filename))
        l2_degrees= []
        l3_degrees= []
        with open(l2_filename, "r") as f:
            for l in f:
                l2_degrees.append(int(l))
        with open(l3_filename, "r") as f:
            for l in f:
                l3_degrees.append(int(l))
        if sum(l2_degrees) != sum(l3_degrees):
            error("sums of L2/L3 degrees do not match")

    elif len(params) == 0:
        info("  Obtaining degrees from bipartite graph")
        if (g == None):
            error("no graph to determine prescribed degree sequences")
        if not(g.is_bipartite()):
            error("command \"prescribed\" can only be used with bipartite graphs")
        l2_vertices= [idx for idx in range(g.vcount()) if not(g.vs["type"][idx])]
        l3_vertices= [idx for idx in range(g.vcount()) if g.vs["type"][idx]]
        l2_degrees= g.degree(l2_vertices)
        l3_degrees= g.degree(l3_vertices)
        
    else:
        error("invalid number of arguments for \"prescribed\"")

    info("  Running bipartite configuration model")
    return bipartite_cm(l2_degrees, l3_degrees)


# -----[ bipartite_op_plot ]-----------------------------------------
#
# -------------------------------------------------------------------
def bipartite_op_plot(g, params):
    if (len(params) < 1):
        filename= None
    else:
        filename= params[0]
    layout= g.layout("fr")
    #layout= g.layout("kk")
    if "type" in g.vs.attribute_names():
        g.vs["color"]= [COLOR_DICT[type] for type in g.vs["type"]]
        g.vs["shape"]= [SHAPE_DICT[type] for type in g.vs["type"]]
    else:
        g.vs["color"]= [COLOR_DICT[TYPE_L3] for idx in range(g.vcount())]
        g.vs["shape"]= [SHAPE_DICT[TYPE_L3] for idx in range(g.vcount())]
    if filename:
        igraph.drawing.plot(g, layout=layout, target=filename)
    else:
        igraph.drawing.plot(g, layout=layout)
    return g


def save_bipartite_edgelist(g, filename):
    types= g.vs["type"]
    with open(filename, "w") as f:
        for v in range(g.vcount()):
            print >> f, v, STR_DICT[types[v]]
        print >> f
        for e in g.get_edgelist():
            print >> f, e[0], e[1]
    return


def load_bipartite_edgelist(filename):
    edge_types= {}
    edges= []
    state= 0
    with open(filename, "r") as f:
        for line in f:
            line= line.rstrip("\r\n")
            assert(state in [0, 1])
            if len(line) == 0:
                state+= 1
                continue
            fields= line.split()
            assert(len(fields) == 2)
            if (state == 0):
                vertex= int(fields[0])
                assert(fields[1] in ["L2", "L3"])
                if fields[1] == "L2":
                    level= TYPE_L2
                else:
                    level= TYPE_L3
                edge_types[vertex]= level
            elif (state == 1):
                source= int(fields[0])
                target= int(fields[1])
                edges.append((source, target))
                
    types= [edge_types[i] for i in sorted(edge_types)]
    g= igraph.Graph.Bipartite(types, edges)
    return g


# -----[ bipartite_op_load ]-----------------------------------------
#
# -------------------------------------------------------------------
def bipartite_op_load(g, params):
    filename= params[0]
    format= params[1]
    info("  Load graph from file \"%s\" with format \"%s\"" \
          % (filename, format))
    if (format == "gml"):
        g= igraph.Graph.Read_GML(filename)
        for v in g.vs:
            v["type"]= TYPE_L3
    elif (format == "edgelist"):
        g= igraph.Graph.Read_Edgelist(filename)
        for v in g.vs:
            v["type"]= TYPE_L3
    elif (format == "edgelistb"):
        g= load_bipartite_edgelist(filename)
    elif (format == "graphml"):
        f= igraph.Graph.Read_GraphML(filename)
        for v in g.vs:
            v["type"]= TYPE_L3
    elif (format == "graphmlz"):
        f= igraph.Graph.Read_GraphMLz(filename)
        for v in g.vs:
            v["type"]= TYPE_L3
    else:
        error("unknown file format \"%s\"" % (format))
        return None
    return g


# -----[ bipartite_op_save ]-----------------------------------------
#
# -------------------------------------------------------------------
def bipartite_op_save(g, params):
    if (g == None):
        error("no graph to save")
    filename= params[0]
    format= params[1]
    info("  Save graph to file \"%s\" with format \"%s\"" \
          % (filename, format))
    if (format == "gml"):
        g.write_gml(filename)
    elif (format == "graphml"):
        g.write_graphml(filename)
    elif (format == "graphmlz"):
        g.write_graphmlz(filename)
    elif (format == "dot"):
        g.write_dot(filename)
    elif (format == "edgelist"):
        g.write_edgelist(filename)
    elif (format == "edgelistb"):
        save_bipartite_edgelist(g, filename)
    else:
        error("unknown file format \"%s\"" % (format))
    return g


def __degrees_summary(degrees, prefix):
    degrees= sorted(degrees)
    if len(degrees) > 0:
        print "%s.degree.min\t%d" % (prefix, degrees[0])
        print "%s.degree.max\t%d" % (prefix, degrees[-1])
        s= sum(degrees)
        print "%s.degree.avg\t%f" % (prefix, float(s)/len(degrees))
        print "%s.degree.total\t%d" % (prefix, s)
    else:
        print "%s.degree.min\t%d" % (prefix, 0)
        print "%s.degree.max\t%d" % (prefix, 0)
        print "%s.degree.avg\t%f" % (prefix, 0)
        print "%s.degree.total\t%d" % (prefix, 0)
    return


# -----[ bipartite_op_summary ]--------------------------------------
#
# -------------------------------------------------------------------
def bipartite_op_summary(g, params):
    types= g.vs["type"]
    l2_vertices= [idx for idx in range(g.vcount()) \
                  if types[idx] == False]
    l3_vertices= [idx for idx in range(g.vcount()) \
                  if types[idx] == True]
    l2l3_edges= [e for e in g.get_edgelist() \
                 if not(types[e[0]] and types[e[1]])]
    l3l3_edges= [e for e in g.get_edgelist() \
                 if (types[e[0]] and types[e[1]])]
    print "vcount\t%d" % (g.vcount())
    print "l2.vcount\t%d" % (len(l2_vertices))
    print "l3.vcount\t%d" % (len(l3_vertices))
    print "ecount\t%d" % (g.ecount())
    print "l2l3.ecount\t%d" % (len(l2l3_edges))
    print "l3l3.ecount\t%d" % (len(l3l3_edges))
    clusters= g.clusters()
    print "ccount\t%d" % (len(clusters))
    print "giant.vcount\t%d" % (clusters.giant().vcount())
    print "giant.ecount\t%d" % (clusters.giant().ecount())
    __degrees_summary(g.degree(l2_vertices), "l2")
    __degrees_summary(g.degree(l3_vertices), "l3")
    l3l2_degrees= {}
    for e in l2l3_edges:
        if types[e[0]]:
            v= e[0]
        else:
            v= e[1]
        if not(v in l3l2_degrees):
            l3l2_degrees[v]= 1
        else:
            l3l2_degrees[v]+= 1
    __degrees_summary(l3l2_degrees.values(), "l3l2")
    l3l3_degrees= {}
    for e in l3l3_edges:
        if not(e[0] in l3l3_degrees):
            l3l3_degrees[e[0]]= 1
        else:
            l3l3_degrees[e[0]]+= 1
        if not(e[1] in l3l3_degrees):
            l3l3_degrees[e[1]]= 1
        else:
            l3l3_degrees[e[1]]+= 1
    __degrees_summary(l3l3_degrees.values(), "l3l3")
    return g


# -----[ bipartite_op_giant ]----------------------------------------
# Extract largest connected component
# -------------------------------------------------------------------
def bipartite_op_giant(g, params):
    if (g == None):
        error("no graph to extract largest component")
    print "  Extract largest connected component"
    return g.clusters().giant()


# -----[ bipartite_op_project ]--------------------------------------
# Project bipartite graph onto L3
# -------------------------------------------------------------------
def bipartite_op_project(g, params):
    if (g == None):
        error("no graph to project")
    print "  Project onto L3"
    (g_l2, g_l3)= g.bipartite_projection()
    for v in g_l3.vs:
        v["type"]= TYPE_L3
    return g_l3


# -----[ bipartite_op_orphans ]--------------------------------------
# Remove orphans (o-degree) nodes
# -------------------------------------------------------------------
def bipartite_op_orphans(g, params):
    if (g == None):
        error("no graph to remove orphans from")
    orphans= [idx for idx in range(g.vcount()) if g.degree(idx)==0]
    g.delete_vertices(orphans)
    print "  Number of orphan nodes removed: %d" % (len(orphans))
    return g


# -----[ _estimate_pareto_params ]-----------------------------------
def _estimate_pareto_params(degrees, beta_range):
    alphas= {}
    for beta in beta_range:
        #alphas[beta]= Pareto.mle(degrees, beta)
        alphas[beta]= Zipf.mle(degrees, beta)
    return alphas

# -----[ _plot_pareto_pdf ]------------------------------------------
def _plot_pareto_pdf(ax, alpha, xmin, degrees, max_degree, hist=False):
    degree_freq= []
    #p= Pareto2(alpha, xmin)
    p= Zipf(alpha, max_degree, xmin)
    n= 0
    for k in degrees:
        if k >= xmin:
            n+= 1
    total= 0
    xrange= range(xmin, max_degree, 1)
    for k in xrange:
        pdf= p.pdf(k)
        if not(hist):
            degree_freq.append(pdf)
        else:
            #degree_freq.append(total+pdf)
            degree_freq.append(p.cdf(k))
        total+= pdf
    degree_freq= [i * n / total for i in degree_freq]
        
    if not(hist):
        ax.plot(xrange, degree_freq)
    else:
        ax.plot([i - 0.5 for i in xrange], degree_freq, 'v')
    return

# -----[ _plot_degree_hist ]-----------------------------------------
def _plot_degree_hist(ax, degrees, max_degree):
    ax.set_title("Degree distribution", fontsize=FONT_SIZE)
    ax.hist(degrees, bins=max_degree, range=(0, max_degree),
            cumulative=False, normed=False, histtype='bar')
    return

# -----[ _plot_degree_dist_loglog ]----------------------------------
def _plot_degree_dist_loglog(ax, degrees, max_degree):
    degrees_hist= numpy.histogram(degrees, bins=max_degree,
                                  range=(0, max_degree))
    ax.set_title("(log-log)", fontsize=FONT_SIZE)
    ax.grid()
    ax.loglog(range(0, max_degree), degrees_hist[0], 'o')    
    return

# -----[ _dump_degrees_to_file ]-------------------------------------
def _dump_degrees_to_file(degrees, filename):
    file= open(filename, 'w')
    for i in range(len(degrees)):
        file.write("%d\n" % (degrees[i]))
    file.close()
    return

# -----[ bipartite_op_stats_flat ]-----------------------------------
def _bipartite_op_stats_flat(g, output=None):
    # Vertex/edge count
    info("  Number of vertices: %d" % (g.vcount()))
    info("  Number of edges: %d" % (g.ecount()))
    degrees= g.degree()
    max_degree= g.maxdegree()
    info("  Highest degree: %d" % (max_degree))

    # Connectedness
    if (g.is_connected()):
        info("  Graph is connected")
    else:
        info("  Graph is not connected")
    clusters= g.clusters()
    info("  Number of connected components = %d" % (len(clusters)))
    giant= clusters.giant()
    info("  Size of largest component = %d" % (giant.vcount()))

    # Power law exponent estimators
    #mle= Pareto.mle(degrees)
    mle= Zipf.mle(degrees)
    info("  Maximum Likelihood Estimator:")
    info("    alpha (xm=1): %.20f" % (mle))

    # Node degrees distributions
    if (output):
        fig=plt.figure()
        ax= plt.subplot(1, 1, 1)
        _plot_degree_hist(ax, degrees, max_degree)
        filename= "degree-hist-%s" % (output)
        info("  Create \"%s\"" % (filename))
        fig.savefig(filename)
        fig=plt.figure()
        ax= plt.subplot(1, 1, 1)
        _plot_degree_dist_loglog(ax, degrees, max_degree)
        filename= "degree-dist-loglog-%s" % (output)
        info("  Create \"%s\"" % (filename))
        fig.savefig(filename)
    else:
        fig= plt.figure()
        ax= plt.subplot(1, 2, 1)
        _plot_degree_hist(ax, degrees, max_degree)
        ax= plt.subplot(1, 2, 2)
        _plot_degree_dist_loglog(ax, degrees, max_degree)
        plt.show()

    return g

# -----[ bipartite_op_stats_bipartite ]------------------------------
def _bipartite_op_stats_bipartite(g, output=None):
    # Vertex/edge count
    l2_vertices= [idx for idx in range(g.vcount()) if not(g.vs["type"][idx])]
    l3_vertices= [idx for idx in range(g.vcount()) if g.vs["type"][idx]]
    info("  Number of L2 vertices: %d" % (len(l2_vertices)))
    info("  Number of L3 vertices: %d" % (len(l3_vertices)))
    info("  Number of edges: %d" % (g.ecount()))

    l2_degrees= g.degree(l2_vertices)
    l3_degrees= g.degree(l3_vertices)

    l3_degrees= [i for i in l3_degrees if i > 0]
    
    max_l2_degree= g.maxdegree(l2_vertices)
    max_l3_degree= g.maxdegree(l3_vertices)
    if (max_l2_degree > max_l3_degree):
        max_degree= max_l2_degree
    else:
        max_degree= max_l3_degree
    info("  Highest L2 degree: %d" % (max_l2_degree))
    info("  Highest L3 degree: %d" % (max_l3_degree))
    l2_degrees_hist= numpy.histogram(l2_degrees, bins=max_degree,
                                     range=(0, max_degree))
    l3_degrees_hist= numpy.histogram(l3_degrees, bins=max_degree,
                                     range=(0, max_degree))

    # Connectedness
    if (g.is_connected()):
        info("  Graph is connected")
    else:
        info("  Graph is not connected")
    clusters= g.clusters()
    info("  Number of connected components = %d" % (len(clusters)))
    giant= clusters.giant()
    info("  Size of largest component = %d" % (giant.vcount()))

    # L2/L3 power law exponents estimators
    info("  Maximum Likelihood Estimator:")
    l2_xmin_range= range(1, 4)
    l2_mle= _estimate_pareto_params(l2_degrees, l2_xmin_range)
    for xmin in l2_xmin_range:
        info("    alpha (L2, xm=%d): %.20f" % (xmin, l2_mle[xmin]))
    #l3_mle= Pareto.mle(l3_degrees)
    l3_mle= Zipf.mle(l3_degrees)
    info("    alpha (L3, xm=1): %.20f" % (l3_mle))

    # L2/L3 degrees distributions
    if (output):
        fig= plt.figure()
        ax= plt.subplot(1, 1, 1)
        _plot_degree_hist(ax, l2_degrees, max_degree)
        for xmin in l2_xmin_range:
            _plot_pareto_pdf(ax, l2_mle[xmin], xmin, l2_degrees, max_degree, hist=False)
        filename= "l2-degree-hist-%s" % (output)
        info("  Create \"%s\"" % (filename))
        fig.savefig(filename)
        
        fig= plt.figure()
        ax= plt.subplot(1, 1, 1)
        _plot_degree_dist_loglog(ax, l2_degrees, max_degree)
        for xmin in l2_xmin_range:
            _plot_pareto_pdf(ax, l2_mle[xmin], xmin, l2_degrees, max_degree)
        filename= "l2-degree-dist-loglog-%s" % (output)
        info("  Create \"%s\"" % (filename))
        fig.savefig(filename)
        
        fig= plt.figure()
        ax= plt.subplot(1, 1, 1)
        _plot_degree_hist(ax, l3_degrees, max_degree)
        _plot_pareto_pdf(ax, l3_mle, 1, l3_degrees, max_degree, hist=False)
        filename= "l3-degree-hist-%s" % (output)
        info("  Create \"%s\"" % (filename))
        fig.savefig(filename)
        
        fig= plt.figure()
        ax= plt.subplot(1, 1, 1)
        _plot_degree_dist_loglog(ax, l3_degrees, max_degree)
        _plot_pareto_pdf(ax, l3_mle, 1, l3_degrees, max_degree)
        filename= "l3-degree-dist-loglog-%s" % (output)
        info("  Create \"%s\"" % (filename))
        fig.savefig(filename)
    else:
        fig= plt.figure()
        ax= plt.subplot(2, 2, 1)
        _plot_degree_hist(ax, l2_degrees, max_degree)
        ax= plt.subplot(2, 2, 2)
        _plot_degree_dist_loglog(ax, l2_degrees, max_degree)
        for xmin in l2_xmin_range:
            _plot_pareto_pdf(ax, l2_mle[xmin], xmin, l2_degrees, max_degree)

        ax= plt.subplot(2, 2, 3)
        _plot_degree_hist(ax, l3_degrees, max_degree)
        ax= plt.subplot(2, 2, 4)
        _plot_degree_dist_loglog(ax, l3_degrees, max_degree)
        _plot_pareto_pdf(ax, l3_mle, 1, l3_degrees, max_degree)
        plt.show()

    return g

# -----[ bipartite_op_stats ]----------------------------------------
# Compute some statistics:
# - number of vertices/edges
# - connectedness (number of connected component + size of giant component)
# - L2/L3 degree frequency distribution
# - Pareto maximum likelyhood estimators
# -------------------------------------------------------------------
def bipartite_op_stats(g, params):
    if (g == None):
        error("no graph to compute statistics")
    output= None
    if (len(params) > 0):
        output= params[0]

    # Is this a multi-graph ? Number of parallel edges ?
    edict= {}
    for e in g.get_edgelist():
        eid= "%d_%d" % (e[0], e[1])
        if not(eid in edict):
            edict[eid]= 1
        else:
            edict[eid]+= 1
    mecount= 0
    max_emultiplicity= 0
    for eid in edict:
        if edict[eid] > 1:
            mecount+= 1
            if edict[eid] > max_emultiplicity:
                max_emultiplicity+= 1
    info("  Number of multi-edges: %d" % (mecount))
    info("  Maximum edge multiplicity: %d" % (max_emultiplicity))

    # Existence of loops ?
    loops= 0
    for e in g.get_edgelist():
        if e[0] == e[1]:
            loops+= 1
    info("  Number of self-loops: %d" % (loops))
            
    bipartite= g.is_bipartite()
    if (bipartite):
        return _bipartite_op_stats_bipartite(g, output)
    else:
        return _bipartite_op_stats_flat(g, output)


# -----[ bipartite_op_mrinfo ]---------------------------------------
# Load an 'mrinfo' file
# -------------------------------------------------------------------
def bipartite_op_mrinfo(g, params):
    filename= params[0]
    typ= params[1]
    info("  Load graph \"%s\"" % (filename))
    info("  as \"%s\"" % (typ))
    mi= MRINFO(filename)
    g= mi.to_graph(typ)
    if g == None:
        error("unknown graph type \"%s\"" % (typ))
    return g


# -----[ bipartite_op_toy ]------------------------------------------
# Generates a toy bipartite graph.
# This is mainly used for demo/debugging.
# -------------------------------------------------------------------
def bipartite_op_toy(g, params):
    info("  Generate toy bipartite graph")
    return bipartite_cm([2, 2, 3, 3, 2], [1, 1, 2, 2, 1, 2, 3])


# -----[ bipartite_op_degrees ]--------------------------------------
# Extract the degrees (L2/L3/*) from the current graph
# -------------------------------------------------------------------
def bipartite_op_degrees(g, params):
    if g.is_bipartite():
        level= params[0]
        if (level == "2"):
            level= TYPE_L2
        elif (level == "3"):
            level= TYPE_L3
        else:
            error("unknown level \"%s\"" % (level))
        info("  Extract L%s degrees from bipartite graph" % (params[0]))
        vertices= [idx for idx in range(g.vcount()) \
                   if (g.vs["type"][idx] == level)]
        degrees= g.degree(vertices)
    else:
        info("  Extract degrees from flat graph")
        degrees= g.degree()

    for d in degrees:
        print "%d" % (d)
    return g


#####################################################################
#
# Main program: definition of supported operations, command-line
# parsing, help.
#
#####################################################################


# -------------------------------------------------------------------
# Define supported graph transformations/operations
#   1st param = function of (graph, typle of params), must return
#               new graph or None
#   2nd param = number of arguments (min, max)
#   3rd param = synopsis of parameters
#   4th param = description
# -------------------------------------------------------------------
operations= {
    "degrees":(
        bipartite_op_degrees, (0, 1),
        "[:LEVEL]",
        [
            "Extract the degrees from the current graph. If the graph is",
            "bipartite, the LEVEL parameter must be either 2 or 3 to",
            "specify which degrees are extracted."
        ]
        ),
    "genpareto":(
        bipartite_op_gen_pareto, (6, 6),
        ":n2:A2:B2:n3:A3:B3",
        [
            "Produces a bipartite graph with the L2/L3 degrees obtained",
            "by sampling 2 Pareto distributions.",
            "  n2: number of L2 vertices (integer, > 0)",
            "  (A2,B2): parameters of the L2 Pareto distribution,",
            "           A2 float >1, B2 float >=1",
            "  n3: number of L3 vertices (integer, > 0)",
            "  (A3,B3): parameters of the L3 Pareto distribution,",
            "           A3 float >1, B2 float >=1",
            "one of the above parameters can be left undefined (-1)"
        ]
        ),
    "genzipf":(
        bipartite_op_gen_zipf, (8, 8),
        ":A2:N2:XM2:n2:A3:N3:XM3:n3",
        [
            "Produces a bipartite graph with L2/L3 degrees obtained",
            "by sampling 2 Zipf distributions.",
            "  (A2,N2,XM2): parameters of the L2 Zipf distribution,",
            "               A2 float >1, N2 integer >0, XM2 integer >0",
            "  n2: number of L2 vertices (integer > 0)",
            "  (A3,N3,XM3): parameters of the L3 Zipf distribution",
            "               A3 float >1, N3 integer >0, XM3 integer >0",
            "  n3: number of L3 vertices (integer, > 0)",
            "the n2 or n3 parameter can be left undefined (<=0)"
        ]
        ),
    "giant":(
        bipartite_op_giant, (0, 0),
        "",
        [
            "Extract the largest connected component of a pre-computed graph."
        ]
        ),
    "load":(
        bipartite_op_load, (2, 2),
        ":FILENAME:FORMAT",
        [
            "Load a graph from a file. The file FORMAT can be one of ",
            "gml, graphml, graphmlz, edgelist, edgelistb."
        ]
        ),
    "mrinfo":(
        bipartite_op_mrinfo, (2, 2),
        ":FILENAME:TYPE",
        [
            "Load a graph in the MRINFO format. The TYPE parameter",
            "selects how the MRINFO must be loaded. The TYPE parameter can be",
            "one of %s." % (", ".join( MRINFO.get_types()))
        ]
        ),
    "orphans":(
        bipartite_op_orphans, (0, 0),
        "",
        [
            "Remove orphan vertices (O-degree vertices) from the current",
            "graph."
        ]
        ),
    "plot":(
        bipartite_op_plot, (0, 1),
        "[:FILENAME]",
        [
            "Plot the current graph. If no FILENAME is provided, the graph",
            "is displayed on the screen. If a FILENAME is provided, the graph",
            "is plotted in a file. The file format is derived from the",
            "FILENAME extension (currently, igraph supports the following",
            "extensions: .pdf, .png and .svg - see the igraph documentation",
            "for more details)"
        ]
        ),
    "prescribed":(
        bipartite_op_prescribed, (0, 2),
        "[:L2_FILENAME:L3_FILENAME]",
        [
            "Generates a bipartite graph using the configuration model",
            "with the L2/L3 degrees obtained from a previously loaded",
            "bipartite graph or from 2 files."
        ]
        ),
    "project":(
        bipartite_op_project, (0, 0),
        "",
        [
            "Produces the L3 projection of a pre-computed bipartite graph.",
        ]
        ),
    "save":(
        bipartite_op_save, (2, 2),
        ":FILENAME:FORMAT",
        [
            "Save the current graph into a file. The file FORMAT can be",
            "one of dot, gml, graphml, graphmlz, edgelist, edgelistb"
        ]
        ),
    "stats":(bipartite_op_stats, (0, 1)),
    "summary":(
        bipartite_op_summary, (0,0),
        "",
        [
            "Show basic statistics about the current graph: number of L2/L3",
            "vertices, number of edges, min/max/avg L2/L3 degrees, ...",
        ]
        ),
    "toy":(
        bipartite_op_toy, (0, 0),
        "",
        [
            "Generates a toy bipartite graph used for demo/debugging."
        ]
        ),
    }


# -----[ usage ]-----------------------------------------------------
# Display program usage help
# -------------------------------------------------------------------
def usage():
    info("Usage: python bipartite.py [ OPERATION [ OPERATION [ ... ] ] ]")
    info("")
    info("  where OPERATION is one of")
    info("")

    for op in sorted(operations):
        op_spec= operations[op]
        if len(op_spec) == 4:
            param_doc= op_spec[2]
            doc= op_spec[3]
            info("  - %s%s" % (op, param_doc))
            info("")
            for l in doc:
                info("     %s" % (l))
        info("")
    return


# -----[ header ]----------------------------------------------------
# Display program purpose, copyright, license and version
# -------------------------------------------------------------------
def header():
    info("Power-law based random bipartite graph generation")
    info("(c) 2011, B. Quoitin, University of Mons, Belgium")
    info("This program is released under the LGPL license")
    info("Version: 0.1")
    info("Note: this program heavily relies on 'igraph'")
    info("")
    return
             

# -----[ main ]-----------------------------------------------------
# This is where everything starts...
# -------------------------------------------------------------------
def main():
    global verbosity
    g= None

    if len(sys.argv) <= 1:
        header()
        usage()
        sys.exit(0)

    # Execute, in sequence, each graph transformation/operation
    # requested by the user on the command-line.
    for op in sys.argv[1:]:
        op_fields= op.split(":")
        op= op_fields[0]
        
        found= False
        for sop in operations.keys():
            if sop == op[0:len(sop)]:
                ndashes= ''.join(['-' for i in range(79-6-len(sop))])
                info("=[ %s ]=%s" % (sop.upper(), ndashes))
                func= operations[sop][0]
                min_args= operations[sop][1][0]
                max_args= operations[sop][1][1]
                if (len(op_fields[1:]) < min_args):
                    error("operation \"%s\" requires %d argument(s)" % \
                          (sop, min_args))
                elif (len(op_fields[1:]) > max_args):
                    error("operation \"%s\" takes max. %d argument(s)" % (sop, max_args))
                g= func(g, op_fields[1:])
                found= True
                break
        if not(found):
            error("unknown operation \"%s\"" % (op))
        info("")
                
    return


if __name__ == "__main__":
    import sys
    main()
