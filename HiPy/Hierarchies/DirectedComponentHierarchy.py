# Copyright ESIEE (2015) 
# 
# benjamin.perret@esiee.fr
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

from HiPy.Structures  import Tree, AdjacencyEdgeWeightedGraph


from HiPy.Util.UnionFind import *
from HiPy.Hierarchies.ComponentTree import addAttributeChildren, addAttributeArea, addAttributeElongationOrientation2d, addAttributeSimpleMoments2d, addAttributeInertia2d
from math import * #@UnusedWildImport



# Computes the directed component hierarchy of a vertex weighted graph
# ############################################################################# 

# Returns :
#    - the parent relation between the strong connected components of the different upper level sets of the graph
#    - the adjacency relation between the strong connected components of each upper level sets of the graph
def DirectedComponentHierarchy(image):

    graph=image.adjacency
    
    # first label of strong components 
    currentLabel=graph.nbPoints
    # parent relation
    parent=[]
    # adjacency relations
    adj=[]   
    # max weight
    lvlMax=max(image.data)
    # variables related to the previous level of the hierarchy
    graphPred, labelsPred, currentLabelPred, numLabelsPred = None, None, None, None
    
    # Let's rock !
    for k in range(lvlMax,-1,-1):
        
        # graph induced by vertices of weight >= k
        gk=builgGraph(graph,image,k)
        
        # strong component labelling, labels start at currentLabel
        # lastLabel is the number of the last label used
        labels,lastLabel=strongComponentLabelings(gk, currentLabel)
        # number of labels in the labelling of the level k
        numLabels = lastLabel - currentLabel
        
        # create adjacency relation of the DAG os strong component
        adjk=adjacencyDAG(gk,labels,currentLabel,numLabels)
        adj.extend(adjk)
        
        if k!= lvlMax:
            # create parent relation between current and previous level
            parentk=parentRelation(graphPred, labelsPred, currentLabelPred, numLabelsPred, labels)
            parent.extend(parentk)
        
        # save current level information for future use
        graphPred, labelsPred, currentLabelPred, numLabelsPred = gk, labels, currentLabel, numLabels
        currentLabel = lastLabel
    
    # parents of the last level (if the stack is connected then numLabels==1: the last label is the root of the tree)
    parent.extend([-1]*numLabels)
    
    return parent, adj
    
# Build the graph induced by vertices of "graph" of weight (given by "image") greater than or equal to "level"
# Returns a graph according to the internal representation of this program
def builgGraph(graph, image, level):
    dim=graph.nbPoints
    
    # create a graph with dim vertices
    # some vertices will certainly be unused.
    # using this trick we are sure that a same vertice get the same index in all levels
    g=AdjacencyEdgeWeightedGraph(dim) 
    g.enabled=[False]*dim
    
    for i in range(dim):
        # if pixel grey level is greater than or equal to level
        if image[i]>=level:
            # mark vertice as used
            g.enabled[i]=True
            # get through all edges going out from vertice i in graph
            for j in graph.getSuccesors(i):
                if image[j]>=level:
                    g.createEdge(i,j)
    return g
                
# Iterative version of Tarjan's algorithm.
# The directed graph "graph" is labelled into strong components.
# Labels are consecutive integers starting at "firstLabel".
# Returns an array containing the strong component label assigned to each vertex of the graph
def strongComponentLabelings(graph, firstLabel=0):
    
    
    # inner variables
    dim=graph.nbPoints
    Index=[-1]*dim
    LowLink=[-1]*dim
    IsInStack=[False]*dim

    ind = 0
    stack = []
    derecStack = []
    
    # results
    Labels=[-1]*dim
    curLabel=firstLabel
    
    # recursive sub-procedure of tarjan strong comp construction              
    def strongConnectComp(base):
        nonlocal ind, curLabel
        derecStack.append([base,graph.out[base]])
        Index[base] = ind
        LowLink[base] = ind
        ind = ind + 1
        IsInStack[base]=True
        stack.append(base)
        while len(derecStack)!=0:
            [i,curAdj]=derecStack.pop() 
            while curAdj!=-1:
                vc = graph.edges[curAdj][1]
                if Index[vc]==-1 : 
                    derecStack.append([i,curAdj])
                    derecStack.append([vc,graph.out[vc]])
                    Index[vc] = ind
                    LowLink[vc] = ind
                    ind = ind + 1
                    IsInStack[vc]=True
                    stack.append(vc)
                    break
                elif IsInStack[vc]:
                    LowLink[i] = min(LowLink[i], Index[vc])
                curAdj=graph.nextEdge[curAdj]

            if curAdj==-1: # all adjacent nodes have been processed
                for vc in graph.getSuccesors(i):
                    if Index[vc] > Index[i] and IsInStack[vc]: 
                        LowLink[i] = min(LowLink[i], LowLink[vc])
                    
                # if true i is the "root" of a strong comp, pop the stack until finding i to construct the strong comp
                if(LowLink[i] == Index[i]):
                    w = None 
                    while w != i:
                        w = stack.pop()
                        IsInStack[w]=False
                        Labels[w]=curLabel
                    curLabel=curLabel+1

       
    for i in range(dim):
        if(graph.enabled[i] and Index[i] == -1):
            strongConnectComp(i)
                     
    return Labels,curLabel

# Construct the adjacency relation on the strong components of the graph "graph"
# "labels" is an array containing the strong component label assigned to each vertex of the graph(see function strongComponentLabelings)
# It is assumed that labels are a serie of "numLalbels" consecutive integers starting at "firstLabel"
# Returns a list of arcs such that [l1 l2] is in  the list iff an arc from l1 to l2 exists in the adjacency relation of strong components
def adjacencyDAG(graph,labels,firstLabel,numLabels):
    dim=graph.nbPoints
    flags=[-1]*numLabels
    sccs=[]
    for i in range(numLabels):
        sccs.append([])
    adj=[]
    for i in range(dim):
        if graph.enabled[i]:
            sccs[labels[i]-firstLabel].append(i)
            
    for l in range(numLabels):
        l1 = l + firstLabel
        for x in sccs[l]:
            for y in graph.getSuccesors(x):
                l2 = labels[y]
                if l1 != l2 and flags[l2 - firstLabel] != l1:
                    adj.append([l1,l2])
                    flags[l2 - firstLabel] = l1
    return adj
        
# Construct the parent relation between two consecutive levels of the stack
# "GraphPred" is the graph of the previous level
# "LabelsPred" is the strong component labelling of "GraphPred"
# It is assumed that labels in "labelsPred" are a serie of "numLabelsPred" consecutive integers starting at "firstLabelPred"
# "Labels" is the strong component labelling of the vertices of "GraphPred" according to the next graph in the stack
# Returns an array containing for each strong component of label in "labelsPred" the label of its parent strong component.
def parentRelation(graphPred, labelsPred, firstLabelPred, numLabelsPred, labels):
    parent = [-1]*numLabelsPred
    dim=graphPred.nbPoints
    for i in range(dim):
        if graphPred.enabled[i]:
            lPar = labels[i]
            l = labelsPred[i]
            parent[l - firstLabelPred] = lPar
    return parent
    
    

# Directed component hierarchy construction from a stack using tarjan union find
############################################################################## 
# "stack" stack of graphs : see functions createGraphStackFromEdgeWeightedGraph and createGraphStackFromVertexWeightedGraph
# "completeGraph" if False then only the strong component tree is built.
#                 if True then the whole directed component hierarchy is built.
# Returns :
#    - the parent relation between the strong connected components of the different upper level sets of the graph
#    - the adjacency relation between the strong connected components of each upper level sets of the graph

def directedComponentHierarchyStackFast(stack, completeGraph=True):
    completeGraphEdges = []
    nbPoints = stack.nbPoints
    
    # pre-allocations of lots of things
    
    # union find makeset
    Par = [i for i in range(nbPoints)]
    parent = [-1] * nbPoints
    Root = [i for i in range(nbPoints)]
    Rank = [0] * nbPoints
    Lvls = [-1] * nbPoints
    # strong comp detection
    Index = [-1] * nbPoints
    LowLink = [-1] * nbPoints
    IsInStack = [False] * nbPoints
    # list of canonical vertices in union find
    activeVertices = []  # nodes
    edges = stack.edges
    i = 0
    c = 0
    while i < len(edges):
        lvl = stack.lvlInsertionSequence[c]
        # add nodes at lvl c
        newnodes = stack.pointInsertionSequence[c]
        c += 1
        for n in newnodes:
            Lvls[n] = lvl
        activeVertices.extend(newnodes)
    
        # add all edges at level lvl (assuming edges are sorted)              
        while i < len(edges) and edges[i][2] == lvl:
            u = edges[i][0]
            v = edges[i][1]
            uc = findTarjan(u, Par)
            vc = findTarjan(v, Par)
            if(uc!=vc):
                stack.setEdge(uc, vc, i)
            i = i + 1
        activeVertices = strongConnectDerecFast(stack, activeVertices, Par, Rank, Root, parent, Index, LowLink, IsInStack, completeGraph, Lvls, lvl)
        if completeGraph:
            saveCompleteGraphEdges(stack, activeVertices, Par, Root, completeGraphEdges)
    return parent, completeGraphEdges, Lvls


# performs a contraction into a DAG of strong connected comps and maintain the strong comp tree,
# returns all the list of canonical nodes representing the strong connected components
def strongConnectDerecFast(graph, activeVertices, Par, Rank, Root, parent, Index, LowLink, IsInStack, completeGraph, Lvls, lvl):
    ind = 0
    stack = []
    derecStack = []
    NActiveVertices = []
    for i in activeVertices:
        Index[i] = -1
        IsInStack[i] = False
        graph.flags[i] = -1
    
    # recursive sub-procedure of tarjan strong comp construction              
    def strongConnectComp(base):
        nonlocal ind, completeGraph
        derecStack.append([base,graph.out[base]])

        Index[base] = ind
        LowLink[base] = ind
        IsInStack[base]=True
        stack.append(base)
        ind = ind + 1

        while len(derecStack)!=0:
            [i,curAdj]=derecStack.pop()

            while curAdj != -1:
                vc = findTarjan(graph.edges[curAdj][1], Par)
                if Index[vc]==-1: # first touch for this adjacent vertex
                    derecStack.append([i,curAdj])
                    derecStack.append([vc,graph.out[vc]])
                    Index[vc] = ind
                    LowLink[vc] = ind
                    ind = ind + 1
                    IsInStack[vc]=True
                    stack.append(vc)
                    break
                elif  IsInStack[vc]: # this node has already been touched by someone else and is not in another scc
                    LowLink[i] = min(LowLink[i], Index[vc])
                
                curAdj=graph.nextEdge[curAdj]
    
            if curAdj==-1: # all adjacent nodes have been processed
                e = graph.out[i]
                # compute lowlink
                while e != -1:
                    vc = findTarjan(graph.edges[e][1], Par)
                    # adjacent node that were not already processed during first pass and that have not been added in another scc
                    if Index[vc] > Index[i] and IsInStack[vc]: 
                        LowLink[i] = min(LowLink[i], LowLink[vc])

                    e = graph.nextEdge[e]
                
                # if true i is the "root" of a strong comp, pop the stack until finding i to construct the strong comp
                if(LowLink[i] == Index[i]):
                    tmp = i
                    w = stack.pop()
                    IsInStack[w]=False
                    if(w != i) or completeGraph:  # else scc has not changed
                        newNode = len(parent)
                        parent.append(-1)
                        Lvls.append(lvl)
                        Root.append(newNode)
                        parent[Root[i]] = newNode
                        while w != i:
                            tmp,tmp2 = unionTarjan(tmp, w, Par, Rank)
                            graph.concatEdgeOut(tmp, tmp2)  # perhaps we can be more fancy here
                            parent[Root[w]] = newNode
                            w = stack.pop()
                            IsInStack[w]=False
                        Root[tmp] = newNode
                        cleanOutEdges(graph, tmp, Par)       
                    NActiveVertices.append(tmp)  # tmp is the canonical node, representative of the strong comp in the DAG

       
    for i in activeVertices:
        if(Index[i] == -1):
            strongConnectComp(i)
                     
    return NActiveVertices




# remove reflexive and redundant out edges from vertice i                     
def cleanOutEdges(graph, i, Par):
    e = graph.out[i]
    while e != -1:
        nexte = graph.nextEdge[e]
        vc = findTarjan(graph.edges[e][1], Par)
        if vc == i or graph.flags[vc] == i:
            graph.removeEdgeOut(i, e)
        else:
            graph.flags[vc] = i
        e = nexte


# save the edges of the dag of strong connected components at the current level in the list completeGraphEdges
def saveCompleteGraphEdges(graph, activeVertices, Par, Root, completeGraphEdges):
    for i in activeVertices:
        source = Root[i]
        e = graph.out[i]
        while e != -1:
            vc = findTarjan(graph.edges[e][1], Par)
            dest = Root[vc]
            completeGraphEdges.append([source, dest])
            e = graph.nextEdge[e]


############################################################################## 
############################################################################## 
############################################################################## 
############################################################################## 
############################################################################## 

### Completely Naive algorithm for any stack of graphs

############################################################################## 
############################################################################## 
############################################################################## 
############################################################################## 
############################################################################## 


# Directed component hierarchy construction from a stack with no optimization
############################################################################## 
# "stack" stack of graphs : see functions createGraphStackFromEdgeWeightedGraph and createGraphStackFromVertexWeightedGraph
# Returns :
#    - the parent relation between the strong connected components of the different upper level sets of the graph
#    - the adjacency relation between the strong connected components of each upper level sets of the graph
#    - the number of arcs going out of each strong connected component

def directedComponentHierarchyStack(stack):
    nbPoints = stack.nbPoints
    activeVertices = []  # nodes
    parent = [-1] * nbPoints
    Lvls = [-1] * nbPoints
    nbOut =[0] * nbPoints
    
    Index = [-1] * nbPoints
    LowLink = [-1] * nbPoints
    IsInStack = [0] * nbPoints
    strongCompLvL = [i for i in range(nbPoints)]
    strongCompLvLPrec = [i for i in range(nbPoints)]

    dagEdges =[]
    
    edges = stack.edges
    i = 0
    c = 0

    while i < len(edges):
        lvl = stack.lvlInsertionSequence[c]
        # add nodes at lvl c
        newnodes = stack.pointInsertionSequence[c]
        c += 1
        for n in newnodes:
            Lvls[n] = lvl
        activeVertices.extend(newnodes)

        # add all edges at level lvl (assuming edges are sorted)
              
        while i < len(edges) and edges[i][2] == lvl:
            u = edges[i][0]
            v = edges[i][1]
            stack.setEdge( u, v, i)
            i = i + 1
        
        strongConnectDerec(stack,activeVertices, Index, LowLink, IsInStack, strongCompLvL,dagEdges,parent,Lvls,lvl,nbOut)
        
        
        for j in range(nbPoints):
            parent[strongCompLvLPrec[j]]=strongCompLvL[j]
        
        strongCompLvLPrec, strongCompLvL =  strongCompLvL, strongCompLvLPrec         

       
    return parent, dagEdges, Lvls, nbOut 


def strongConnectDerec(graph,activeVertices, Index, LowLink, IsInStack, strongCompLvL,dagEdges,parent,Lvls,lvl,nbOut):
    ind = 0
    stack = []
    derecStack = []
    
    for i in activeVertices:
        Index[i] = -1
        IsInStack[i] = False  
        graph.flags[i] = -1
                 
    def strongConnectComp(base):
        nonlocal ind
        derecStack.append([base,graph.out[base]])

        Index[base] = ind
        LowLink[base] = ind
        IsInStack[base]=True
        stack.append(base)
        ind = ind + 1

        def addOutEdges(lst, n):
            oe = graph.out[n]
            while oe != -1:
                v = graph.edges[oe][1]
                lst.append(v)
                oe = graph.nextEdge[oe]
        
        def cleanOutEdgesNaif(i, out):
            marked=[True]*len(out)
            for j in range(len(out)):
                vc = strongCompLvL[out[j]]
                if vc == i:
                    marked[j]=False
                elif graph.flags[vc] == i:
                    nbOut[i] = nbOut[i] +1
                    marked[j]=False
                else:
                    nbOut[i] = nbOut[i] +1
                    graph.flags[vc] = i
            nout = [strongCompLvL[out[j]] for j in range(len(out)) if marked[j] ]
            return nout

        while len(derecStack)!=0:
            [i, curAdj]=derecStack.pop()

                     
            while curAdj != -1:
                vc = graph.edges[curAdj][1]
                if Index[vc]==-1: # first touch for this adjacent vertex
                    derecStack.append([i,curAdj])
                    derecStack.append([vc,graph.out[vc]])
                    Index[vc] = ind
                    LowLink[vc] = ind
                    ind = ind + 1
                    IsInStack[vc]=True
                    stack.append(vc)
                    break

                elif  IsInStack[vc]: # this node has already been touched by someone else and is not in another scc
                    LowLink[i] = min(LowLink[i], Index[vc])
                curAdj=graph.nextEdge[curAdj]

            if curAdj==-1: # all adjacent nodes have been processed

                e = graph.out[i]
                # compute lowlink
                while e != -1:
                    vc = graph.edges[e][1]
                    # adjacent node that were not already processed during first pass and that have not been added in another scc
                    if Index[vc] > Index[i] and IsInStack[vc]: 
                        LowLink[i] = min(LowLink[i], LowLink[vc])
                    e = graph.nextEdge[e]
                
                # if true i is the "root" of a strong comp, pop the stack until finding i to construct the strong comp
                if(LowLink[i] == Index[i]):
                    w = None
                    comp=len(parent)
                    parent.append(-1)
                    graph.flags.append(-1)
                    strongCompLvL[i]=comp
                    Lvls.append(lvl)
                    nbOut.append(0)
                    out = []
           
                    while w != i:
                        w = stack.pop()
                        strongCompLvL[w]=comp
                        IsInStack[w]=False
                        addOutEdges(out,w)
                
                    out = cleanOutEdgesNaif(comp, out)
            
                    for i in out:
                        dagEdges.append([comp, i])
       
    for i in activeVertices:
        if(Index[i] == -1):
            strongConnectComp(i)




# Create a stack from the edge weighted graph 
def createGraphStackFromEdgeWeightedGraph(graph):
    dim=graph.nbPoints
    stack = AdjacencyEdgeWeightedGraph(dim)
    edges = graph.edges
    nbe = len(edges)
    #stack["edges"].extend(edges)
    
    el = [i for i in range(nbe)]
    el.sort(key=lambda x:edges[x][2])
    
    added = [False]*dim
    addedNodes = []
    addedLvls =[]
    stack.pointInsertionSequence = addedNodes
    stack.lvlInsertionSequence = addedLvls
    stackEdges = stack.edges
    
    c=0
    
    while c<len(edges):
        e=edges[el[c]]
        lvl=e[2]
        addedLvls.append(lvl)
        nnodes = []
        while(c<len(edges) and edges[el[c]][2]==lvl):
            e=edges[el[c]]
            stackEdges.append(e)
            if not added[e[0]]:
                nnodes.append(e[0])
                added[e[0]]=True
            if not added[e[1]]:
                nnodes.append(e[1])
                added[e[1]]=True
            c = c + 1
        addedNodes.append(nnodes)
    
    stack.nextEdge = [-1] * len(edges)
    stack.prevEdge = [-1] * len(edges)
    stack.flags = [-1] * dim
    return stack

# Construct a stack for a vertex weighted graph.   
# Weights are given by the image.
def createGraphStackFromVertexWeightedGraph(graph,image):
    dim = len(image)
    stack = AdjacencyEdgeWeightedGraph(dim)
    
    pix = [i for i in range(dim)]
 
    pix.sort(key=lambda x:image[x]) 

    addedNodes = []
    addedLvls =[]
    stack.pointInsertionSequence = addedNodes
    stack.lvlInsertionSequence = addedLvls
    edges = stack.edges
   
    for e in graph.edges:
        x=e[0]
        y=e[1]
        edges.append([x, y, min(image[x],image[y])])
        

    i = dim - 1
    while i >= 0:
        j = pix[i]
        lvl = image[j]
        
        nnodes = []
        while(i >= 0 and image[j] == lvl):
            nnodes.append(j)
            i -= 1
            j = pix[i]
        addedLvls.append(lvl)
        addedNodes.append(nnodes)

    edges.sort(key=lambda x:-x[2])  # lambda x,y: cmp(x[2],y[2]) )

    stack.nextEdge = [-1] * len(edges)
    stack.prevEdge = [-1] * len(edges)
    stack.flags = [-1] * dim
    return stack


############################################################################## 
############################################################################## 
############################################################################## 
############################################################################## 
############################################################################## 

### Usefull functions to finalyze the structure

############################################################################## 
############################################################################## 
############################################################################## 
############################################################################## 
############################################################################## 

# takes a parent relation and the Lvls of each node. 
# If the parent relation does not define a tree, a root is created 
# and each tree of the forest is rooted in this new node.
# The new root is assigned to level -1
def ensureTree(parent, Lvls=[]):
    if parent.count(-1) > 1:
        nel = len(parent)
        for i in range(len(parent)):
            if(parent[i] == -1):
                parent[i] = nel
        parent.append(-1)
        Lvls.append(-1)  # yeark


# Construct a more usable structure from the result of the algorithm directedComponentHierarchyXXXXXX
def buildFinalDCCTree(nbPoints, parent, completeGraphEdges, Lvls, image, nbOut=[] ):
    tree = Tree(parent,Lvls,image)
    tree.addAttribute("sucs", [])
    tree.addAttribute("preds", [])
    tree.addAttribute("nbOut", -1)
    tree.nbOut.data = nbOut
    
    for e in completeGraphEdges:
        tree.sucs[e[0]].append(e[1])
        tree.preds[e[1]].append(e[0])
        
    addAttributeChildren(tree)
    return tree


# Compute the level of each node as its depth in the tree given by parent relation
def computeLevelsAsDepth(parent):
    levels = [0]*len(parent)
    for j in range(len(levels)-1,-1,-1):
        if parent[j]!=-1:
            levels[j] = levels[parent[j]]+1
    return levels           







############################################################################## 
############################################################################## 
############################################################################## 
############################################################################## 
############################################################################## 

### Attributes

############################################################################## 
############################################################################## 
############################################################################## 
############################################################################## 
############################################################################## 


def computeAttributeDirectedComponent(dccTree, attributeName, baseAttributeName, accumulator):
    attr=dccTree.addAttribute(attributeName)
    if attr==None:
        return
    baseAttr=dccTree.getAttribute(baseAttributeName)
    succesors=dccTree.sucs
    visit=[None]*len(dccTree.data)
    
    def computeRec(i,ref):
        visit[i]=ref
        tmp=baseAttr[i]
        for j in succesors[i]:
            if(visit[j]!=ref):
                tmp=accumulator(tmp,computeRec(j,ref))
        return tmp
    
    for i in dccTree.iterateFromLeavesToRoot():
        attr[i]=computeRec(i,i)
            
# Compute the area (number of nodes) of each SCC (attribute "area") and each DCC (attribute "area_directed")                   
def addAttributeAreaDirectedComponent(dccTree):
    addAttributeArea(dccTree)
    computeAttributeDirectedComponent(dccTree,"area_directed","area",lambda x, y : x + y)
    

# Compute the moments [M00 M10 M01 M11 M20 M02] of each SCC (attribute "moments") and each DCC  (attribute "moments_semi")
def addAttributeSimpleMomentsDirectedComponent(dccTree):
    addAttributeSimpleMoments2d(dccTree)
    computeAttributeDirectedComponent(dccTree,"moments_directed","moments",lambda x, y : [ x[i] + y[i] for i in range(len(x))])
    

# Compute the moment of inertia of each SCC (attribute "inertia") and each DCC  (attribute "inertia_semi")
def addAttributeInertiaDirectedComponent(dccTree):
    addAttributeInertia2d(dccTree)
    addAttributeInertia2d(dccTree,"inertia_directed","moments_directed")
    
                           
# Compute the elongation and orientation of each SCC (attributes "elongation" and "orientation") and each DCC  (attribute "elongation_semi" and "orientation_semi")
def addAttributeElongationOrientationDirectedComponent(dccTree):
    addAttributeElongationOrientation2d(dccTree)
    addAttributeElongationOrientation2d(dccTree,"elongation_directed","orientation_directed","moments_directed")

# Compute the attribute value "attributeName_directed" on the graph.
# For each node, the attribute is set to true if the directed connected 
# component rooted int it contains a node having its attribute "attributeName" set to true.
# see function addAttributeMarker
def addAttributeMarkerDirectedComponent(dccTree,attributeName):
    attributeNameDirected = attributeName +"_directed"
    attr=dccTree.getAttribute(attributeName)
    attrDirected=dccTree.addAttribute(attributeNameDirected,False)
    if attrDirected==None:
        return
    
    predecesors = dccTree.preds

    def propagateSemiPred(i):
        for p in predecesors[i]:
            if attrDirected[p]==False:
                attrDirected[p]=True
                propagateSemiPred(p)
                         
    for i in dccTree.iterateFromLeavesToRoot(): 
        if attr[i] and not attrDirected[i]:
            attrDirected[i]=True
            propagateSemiPred(i)

# regularize the result of a selection criterion using the select max strategy
def regularizeSelectMax(dccTree):
    deleted=dccTree.getAttribute("deleted")
    succesors=dccTree.getAttribute("sucs")
    def keepDCC(i):
        deleted[i] = False
        for c in succesors[i]:
            if deleted[c]:
                keepDCC(c)

    for i in dccTree.iterateFromLeavesToRoot():
        if not deleted[i]:   
            keepDCC(i)



# regularize the result of a selection criterion using the discard min strategy              
def regularizeDiscardMin(dccTree):
    deleted=dccTree.getAttribute("deleted")
    succesors=dccTree.getAttribute("sucs")

    def removeDCC(i):
        deleted[i] = True
        for c in succesors[i]:
            removeDCC(c)

    for i in dccTree.iterateFromLeavesToRoot():
        if deleted[i]:   
            removeDCC(i)
