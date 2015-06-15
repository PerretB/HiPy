# Copyright ESIEE (2015) 
# 
# benjamin.perret@esiee.fr
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
from HiPy.Util import UnionFind

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
Created on 15 juin 2015

@author: perretb
'''

import HiPy.Util.UnionFind as UnionFind
from HiPy.Structures import Tree, AdjacencyEdgeWeightedGraph
from HiPy.Hierarchies.ComponentTree import addAttributeChildren

def constructAltitudeBPT(image, verbose=False):
    
    if verbose:
        print("Sorting") 
    edges = image.adjacency.edges
    edges = sorted(edges, key=lambda x:x[2])
    
    if verbose:
        print("Kruskaling") 
    MST,parent=computeMSTBPT(image,edges)
    
    if verbose:
        print("Stuffing")
    levels=[-1]*len(image)
    adjMST=AdjacencyEdgeWeightedGraph(len(image))
    for e in MST:
        levels.append(e[2])
        adjMST.createEdge(*e)
    
    if verbose:
        print("Finalizing")
    tree = Tree(parent, levels, image)
    tree.leavesAdjacency=adjMST
    return tree

def transformAltitudeBPTtoComponentTree(bpt):
    '''
    Delete node n such that bpt.level[n]=bpt.level[bpt[n]] (and update the parent relation accordingly...
    
    Warning: this function returns nothing, the argument is modified!
    '''
    nbLeaves = bpt.nbLeaves
    nbNodes = len(bpt)
    addAttributeChildren(bpt)
    children=bpt.children
    level=bpt.level
    
    count=0
    deleted=[False]*nbNodes
    dmap=[0]*nbNodes
    
    #from root to leaves, compute the new parent relation, don't care of the holes in the parent tab
    for i in range(nbNodes-2,nbLeaves-1,-1):
        p = bpt[i]
        if level[i]==level[p]:
            for c in children[i]:
                bpt[c]=p
            deleted[i]=True
            count=count+1
        #inverse of what we want: number of deleted nodes after node i
        dmap[i] = count
    
    #correct the mapping 
    for i in bpt.iterateFromLeavesToRoot(False):
        dmap[i] = count - dmap[i]
    
    #new relations with correct size  
    nparent = [0]*(nbNodes-count)
    nlevel = [-1]*(nbNodes-count) 
    
    count=0
    for i in range(0,nbNodes-1):
        if not deleted[i]:
            p = bpt[i]
            np = p - dmap[p] 
            nparent[count]=np
            nlevel[count]=level[i]
            count = count + 1
    nparent[count]=-1
    
    del bpt[:]
    del level[:]
    for i in range(len(nparent)):
        bpt.append(nparent[i])
        level(nlevel[i])

def computeMSTBPT(image,sortedEdgeList):
    '''
    precondition: The edge list must be sorted
    '''
    nbEdgeMST=len(image)-1
    mst=[]
    parent = [-1]*len(image)
    ufParent = [None]*len(image)
    ufRank = [0]*len(image)
    root = [i for i in range(len(image))]
    nbEdgeFound=0
    i=0
    while nbEdgeFound<nbEdgeMST :
        e = sortedEdgeList[i]
        c1 = UnionFind.findTarjan(e[0],ufParent)
        c2 = UnionFind.findTarjan(e[1],ufParent)
        if c1!=c2:
            newParent=len(parent)
            parent.append(-1)
            parent[root[c1]]=newParent
            parent[root[c2]]=newParent
            nr,_ = UnionFind.unionTarjan(c1, c2, ufParent, ufRank)
            root[nr]=newParent
            mst.append(e)
            nbEdgeFound = nbEdgeFound+1
        i=i+1
    
    return mst, parent

