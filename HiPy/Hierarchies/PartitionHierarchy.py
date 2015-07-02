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

'''
Created on 15 juin 2015

@author: perretb
'''

import HiPy.Util.UnionFind as UnionFind
from HiPy.Structures import Tree, Adjacency2d4, Embedding2dGrid, Image, TreeType, WeightedAdjacency
from HiPy.Processing.Attributes import addAttributeChildren,\
    addAttributeDepth, addAttributeRank
from HiPy.Util.Histogram import imageMap, rescaleGray, normalizeToByte

def constructAltitudeBPT(adjacency, verbose=False):
    '''
    Construct the binary partition tree by altitude of the given image.
    
    It assumes that adjacency has weighted edges.
    
    The corresponding MST is stored in the attribute leavesAdjacency of the returned tree
    '''
    if verbose:
        print("Sorting") 
    
    nbPoints = adjacency.nbPoints
    edgesI = sorted(range(len(adjacency)), key=lambda x:adjacency[x])
    edges = [[adjacency.source[i], adjacency.target[i], adjacency[i]]  for i in edgesI]
    if verbose:
        print("Kruskaling") 
    MST,parent=computeMSTBPT(nbPoints,edges)
    
    if verbose:
        print("Stuffing")
    levels=[0]*nbPoints 
    adjMST=WeightedAdjacency(nbPoints)

    for e in MST:
        levels.append(e[2])
        adjMST.createEdge(*e)

    if verbose:
        print("Finalizing")
    tree = Tree(TreeType.PartitionHierarchy,parent, levels)
    tree.leavesAdjacency=adjMST
    return tree

def transformAltitudeBPTtoComponentTree(bpt):
    '''
    Copy bpt and delete the nodes n such that bpt.level[n]=bpt.level[bpt[n]] 
    (and update the parent relation accordingly...
    '''
    nbLeaves = bpt.nbPixels
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
    for i in bpt.iteratorFromPixelsToRoot(False):
        dmap[i] = count - dmap[i]
    
    #new relations with correct size  
    nparent = [-1]*(nbNodes-count)
    nlevel = [0]*(nbNodes-count) 
    
    count=0
    for i in range(0,nbNodes-1):
        if not deleted[i]:
            p = bpt[i]
            np = p - dmap[p] 
            nparent[count]=np
            nlevel[count]=level[i]
            count = count + 1
    nparent[count]=-1
    nlevel[count]=level[-1]
    
    ntree = Tree(TreeType.PartitionHierarchy,nparent, nlevel)
    ntree.leavesAdjacency=bpt.leavesAdjacency
    
    return ntree

def transformBPTtoAttributeHierarchy(bpt, attributeName):
    adj=reweightMSTByAttribute(bpt, attributeName)
    return constructAltitudeBPT(adj)

def transformAltitudeBPTtoWatershedHierarchy(bpt):
    nadj=extractWatershedEdges(bpt) 
    return constructAltitudeBPT(nadj)
    
def computeMSTBPT(nbPoints,sortedEdgeList):
    '''
    precondition: The edge list must be sorted
    '''
    nbEdgeMST=nbPoints-1
    
    mst=[]
    parent = [-1]*nbPoints
    ufParent = [i for i in range(nbPoints)]
    ufRank = [0]*nbPoints
    root = [i for i in range(nbPoints)]
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

def extractWatershedEdges(bpt):
    '''
    Return a copy of bpt.leavesAdjacency (ie the MST) where the weight of the non-watershed edges are set to 0
    '''
    
    nadj=bpt.leavesAdjacency.copy()
    nbLeaves=bpt.nbPixels
    nbNodes=len(bpt)
    addAttributeChildren(bpt)
    level=bpt.level
    children=bpt.children
    
    minima=[0]*nbNodes
    for i in bpt.iteratorFromPixelsToRoot(False):
        cl = children[i]
        m1=minima[cl[0]]
        m2=minima[cl[1]]
        nb=m1+m2

        if m1==0 or m2==0:
            nadj[i-nbLeaves]=0
       
        if nb!=0:
            minima[i]=nb
        else:
            p=bpt[i]
            #if p==-1 there is a unique minima in the image
            if p==-1 or level[i]!=level[p]:
                minima[i]=1
            #else minima[i]=0
    
    return nadj


def correctAttributeValueBPT(bpt,attributeName):
    '''
    Most attributes (@todo precision needed!!!) computed using fonctions in module HiPy.Component 
    will be incorrect for the BPT. However their values can be corrected a posteriori
    using this function !
    ''' 
    addAttributeChildren(bpt)
    children=bpt.children        
    attr=bpt.getAttribute(attributeName)
    level=bpt.level
    for i in bpt.iterateOnLeaves():
        attr[i]=0
    
    for i in bpt.iteratorFromPixelsToRoot(False):
        p = bpt[i]
        if p!=-1 and level[i]==level[p]:
            vc1=attr[children[i][0]]
            vc2=attr[children[i][1]]
            attr[i]=max(vc1,vc2)

def reweightMSTByAttribute(bpt, attributeName, extinctionValue=True): 
    '''
    Reweights the MST associated to a binary partition tree according to a given attribute.
    
    If the attribute represents an extinctionValue, then the new weight of an edge is the 
    minimum of the attribute values of the two regions separated by this edge.
    
    Otherwise the new value of the edge is simply the attribut value of the node associated
    to this edge.
    
    The function returns a new weighted adjacency containing with the same edges as the MST but different weights.
    '''
    nbLeaves=bpt.nbPixels
    nadj=bpt.leavesAdjacency.copy()
    addAttributeChildren(bpt)
    children=bpt.children
    attr=bpt.getAttribute(attributeName)
    if extinctionValue:
        for i in bpt.iteratorFromPixelsToRoot(False):
            nadj[i-nbLeaves]=min(attr[children[i][0]],attr[children[i][1]])
    else:
        for i in bpt.iteratorFromPixelsToRoot(False):
            nadj[i-nbLeaves]=attr[i]
    return nadj

def computeSaliencyMap(partitionTree,adjacency, attribute="level"):
    '''
    Compute the saliency values of the edges of the given adjacency w.r.t the given partition tree.
    
    A new adjacency is created during the process.
    '''
    addAttributeDepth(partitionTree)
    attr=partitionTree.getAttribute(attribute)
    lca=partitionTree.lca
    
    # could be shortened to
    # WeightedAdjacency.createAdjacency(adjacency, lambda i,j:level[lca(i,j)])
    # but doubles the cost due to symmetric weights
#     adj=WeightedAdjacency(adjacency.nbPoints)
#     for i in range(adj.nbPoints):
#         for j in adjacency.getSuccesors(i):
#             if j>i:
#                 w=attr[lca(i,j)]
#                 adj.createEdge(i, j, w)
#                 #adj.createEdge(j, i, w)
#      
#     return adj 

    return WeightedAdjacency.createAdjacency(adjacency, lambda i,j:attr[lca(i,j)])

def drawSaliencyMap(size,saliency):
    '''
    Represent a saliency map as a contour image.
    Size is the size [width, height] of the image => result size [2*width-1,2*height-1].
    Saliency must represent a 4 adjacency on the 2d grid (typically an Adjacency2d4), results are unpredictable otherwise
    '''
    w=size[0]
    h=size[1]
    grid = Embedding2dGrid(w,h)
    rw=w*2-1
    rh=h*2-1
    grid2 = Embedding2dGrid(rw,rh)
    
    res = Image(rw*rh,0,Adjacency2d4([rw,rh]),grid2)
    for y in range(h):
        for x in range(w):
            p=grid.getLinearCoordinate(x,y)
            for n in saliency.getOutEdges(p):
                if n[1]>p:
                    ng=grid.fromLinearCoordinate(n[1])
                    res.setPixelWCS(n[2],x+ng[0],y+ng[1])
    
    for y in range(1,rh-1,2):
        for x in range(1,rw-1,2):
            p=grid2.getLinearCoordinate(x,y)
            vmax=-999999
            for n in res.getNeighbours(p):
                vn=res[n]
                if vn>vmax:
                    vmax=vn
            res[p]=vmax
    
    return res


def drawSaliencyForVizu(tree,image, attr="level", gammaFactor=0.33333):
    '''
    Draw the saliency map associated to the given partition tree in an image.
    
    The tree is assumed to represent the given image.
    
    Only saliceny associated to a 4 adjacency can be drawn consistently.
    '''
    adj4 = Adjacency2d4(image.embedding.size)
    addAttributeRank(tree)
    saliency=computeSaliencyMap(tree,adj4,attr)
    sal=drawSaliencyMap(image.embedding.size,saliency)
    
    sal = rescaleGray(sal, 0, 1)
    sal = imageMap(sal, lambda x:x**gammaFactor)
    sal = normalizeToByte(sal)
    return sal