# Copyright ESIEE (2015) 
# 
# benjamin.perret@esiee.fr
# 
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 

# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 

# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 

# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

'''
Created on 10 juin 2015


@author: perretb
'''

from HiPy.IO import * #@UnusedWildImport
from HiPy.Hierarchies.DirectedComponentHierarchy import * #@UnusedWildImport


# for drawing graphs in a pdf file
try:
    import pydot 
    pydotAvailable = True
except ImportError:
    pydotAvailable = False


# Test of the algorithm "as close as possible" of the one given in the article 
# Input: 
#    - file fig5-Graph.list : plain text file encoding the graph of Fig. 5
#    - file fig5-Image.png: png image encoding the weights of the graph of Fig. 5
# Output:
#    - file "DCC-Hierarchy-Fig5.pdf": pdf plot of the graph Fig. 6 c
def testAlgoArticle():
    print("Test of the paper algorithm for vertex weighted graphs")
    graph = readGraph('../samples/DCC/fig5-Graph.list')
    image=readImage("../samples/DCC/fig5-Image.png")
    image.adjacency=graph
    parent, adj= DirectedComponentHierarchy(image)
    parents = [-1]*graph.nbPoints
    parents.extend(parent)
    levels = computeLevelsAsDepth(parents)

    dccTree = buildFinalDCCTree(graph.nbPoints, parents, adj, levels, image)
    drawGraphVizAttr("Results/DCC-Hierarchy-Fig5.pdf", dccTree, ["level"],False)
    print("=> Result saved in file " + " 'DCC-Hierarchy-Fig5.pdf'")
    print("Done\n\n")


# Test of the naive algorithm for stack 
# Input: 
#    - file fig5-Graph.list: plain text file encoding the graph of Fig. 5
#    - file fig5-Image.png: png image encoding the weights of the graph of Fig. 5
# Output:
#    - file "DCC-Hierarchy-Fig5-Stack.pdf": pdf plot of the graph Fig. 6 c
def testAlgoStackGraphArticle():
    print("Test of the general algorithm for stacks of graphs")
    adj = readGraph("../samples/DCC/fig5-Graph.list")
    image= readImage("../samples/DCC/fig5-Image.png")
    stack = createGraphStackFromVertexWeightedGraph(adj,image)
    
    parent, DAGsAdj, Lvls, _ = directedComponentHierarchyStack(stack)
    ensureTree(parent, Lvls)
    dccTree = buildFinalDCCTree(stack.nbPoints, parent, DAGsAdj, Lvls, image)
    
    drawGraphVizAttr("Results/DCC-Hierarchy-Fig5-Stack.pdf", dccTree, ["level"],False)
    print("=> Result saved in file " + " 'DCC-Hierarchy-Fig5-Stack.pdf'")
    print("Done\n\n")


# Test of the semi-naive algorithm for stack 
# Input: 
#    - file fig5-Graph.list : plain text file encoding the graph of Fig. 5
#    - file fig5-Image.png: png image encoding the weights of the graph of Fig. 5
# Output:
#    - file "DCC-Hierarchy-Fig5-Stack-Fast.pdf": pdf plot of the graph Fig. 6 c
def testAlgoStackGraphArticleFast():
    print("Test of the general algorithm for stacks of graphs with optimization")
    adj = readGraph("../samples/DCC/fig5-Graph.list")
    image = readImage("../samples/DCC/fig5-Image.png")
    stack = createGraphStackFromVertexWeightedGraph(adj,image)

    parent, DAGsAdj, Lvls = directedComponentHierarchyStackFast(stack)
    ensureTree(parent, Lvls)
    dccTree = buildFinalDCCTree(stack.nbPoints, parent, DAGsAdj, Lvls, image)
    
    drawGraphVizAttr("Results/DCC-Hierarchy-Fig5-Stack-Fast.pdf", dccTree, ["level"],False)
    print("=> Result saved in file " + " 'DCC-Hierarchy-Fig5-Stack-Fast.pdf'")
    print("Done\n\n")
    

# draw directed component hierarchy in a pdf with the values of the specified attributes
# "name" : filename of the produced pdf
# "graph" : graph to plot
# "attributes" : list of attributes values to print
# "leaves" : if true leaves are ploted

def drawGraphVizAttr(name, dccTree, attributes=[], leaves=False):
    if not pydotAvailable:
        print ("Error: this function requires the library pydot (which requires the software graphviz.")
        return
    
    #n = len(parent)
    G = pydot.Dot(graph_type='graph')

    nbLeaves = dccTree.nbLeaves
    
    for i in dccTree.iterateFromLeavesToRoot():
        if(i>=nbLeaves or leaves):
            node = pydot.Node(i)
            if(i < nbLeaves):
                node.set_shape("square")
            label = str(i)
            if attributes != []:
                labels = ",".join([str(dccTree.getAttribute(a)[i]) for a in attributes])
                label += "(" + labels + ")"
            node.set_label(label)
            G.add_node(node)
       
    for i in dccTree.iterateFromLeavesToRoot():
        if(i>=nbLeaves or leaves):
            par = dccTree[i]
            if(par != -1):
                edge = pydot.Edge(par, i)
                edge.set_dir("forward")
                G.add_edge(edge)
    
            for d in dccTree.getAttribute("sucs")[i]:
                edge = pydot.Edge(i, d)
                edge.set_dir("forward")
                edge.set_style("dotted")
                edge.set_constraint("false")
                G.add_edge(edge)
       

    G.write_pdf(name)
  
    
def main():
    testAlgoArticle()
    testAlgoStackGraphArticle() 
    testAlgoStackGraphArticleFast()    
    
if __name__ == '__main__':
    main()
    