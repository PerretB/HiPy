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
Created on 3 juin 2015

@author: perretb
'''

from HiPy.Structures import Image, Tree, enum, TreeType
from HiPy.Processing.Attributes import addAttributeChildren
from HiPy.Util.UnionFind import findTarjan
from math import * #@UnusedWildImport




ComponentTreeType  = enum(MinTree = 1, MaxTree = 2)



def constructComponentTree(image,treeType = ComponentTreeType.MaxTree, verbose=False):
    '''
    Construct the min or max tree of the given image
    '''
    # sorting
    if verbose:
        print("Sorting")  
    if treeType == ComponentTreeType.MaxTree:
        sortedPixels = sorted(range(len(image)), key=lambda x:image[x])
    elif treeType == ComponentTreeType.MinTree:
        sortedPixels = sorted(range(len(image)), key=lambda x:-image[x])
    else:
        raise Exception("Component Tree construction error : unknown tree type")
        
    # generic tree construction
    if verbose:
        print("Union Find")
    parent=preTreeConstruction(image,sortedPixels)
     
    # canonize tree
    if verbose:
        print("Canonize")
    canonizeTree(parent,image,sortedPixels)
    parent,levels = expandCanonizedParentRelation(parent,image,sortedPixels)
    
    # encapsulated tree for ease of manipulation
    if verbose:
        print("Tree finalization")
    return Tree(TreeType.ComponentTree,parent,levels,image)


def preTreeConstruction(image, sortedPixels):
    '''
    Generic tree construction using ordered pixels and union find
    '''
    parent = image.copy(False)
    ufParent = Image(len(image),None)
    ufRank = Image(len(image))
    reprez = Image(len(image))
    for i in range(len(image)-1,-1,-1):
        currentPoint = sortedPixels[i]
        parent[currentPoint]=currentPoint
        ufParent[currentPoint]=currentPoint
        reprez[currentPoint]=currentPoint
        cpReprez = currentPoint
        for neighbour in image.getNeighbours(currentPoint):
            if ufParent[neighbour]!=None:
                cpNeigh = findTarjan(neighbour, ufParent)
                if cpNeigh != cpReprez: 
                    parent[reprez[cpNeigh]] = currentPoint
                    if ufRank[cpReprez]<ufRank[cpNeigh]:
                        cpNeigh,cpReprez=cpReprez,cpNeigh
                    ufParent[cpNeigh]=cpReprez
                    reprez[cpReprez]=currentPoint
                    if ufRank[cpReprez]==ufRank[cpNeigh]:
                        ufRank[cpReprez]=ufRank[cpReprez]+1
                    
    return parent


def canonizeTree(parent,enqueuedLevels,sortedPixels):
    '''
    Parent relation "canonization" (path compression) after preTreeConstruction
    '''
    for p in sortedPixels:
        q=parent[p]
        if enqueuedLevels[parent[q]]==enqueuedLevels[q]:
            parent[p]=parent[q]


def expandCanonizedParentRelation(canonizedTreeImage,nodeLevels,sortedPixels):
    '''
    Expand a canonized parent relation to a regular parent relation (each node is represented individually)
    '''
    levels=nodeLevels.copy(True)
    data = [None]*len(canonizedTreeImage)
    for j in range(len(canonizedTreeImage)-1,-1,-1):
        i=sortedPixels[j]
        if(nodeLevels[i]!=nodeLevels[canonizedTreeImage[i]]):
            parent=i
        else:
            parent=canonizedTreeImage[i]
        if(data[parent]==None):
            data.append(-1)
            data[parent]=len(data)-1;
            levels.append(nodeLevels[parent])
        data[i]=data[parent]
        
    for j in range(len(canonizedTreeImage)-1,-1,-1):
        i=sortedPixels[j]
        if(nodeLevels[i]!=nodeLevels[canonizedTreeImage[i]]):
            parent=i
            pparent=canonizedTreeImage[parent]
            data[data[parent]]=data[pparent]
    #data[-1]=-1
    return data,levels            

            
def regularizeSelectMaxTree(tree):
    '''
    Regularize the result of a selection criterion using the select max hierarchical strategy.
    
    A node is mark deleted if and only if it and all his ancestors were marked deleted. 
    '''
    deleted=tree.getAttribute("deleted")
    addAttributeChildren(tree)
    children=tree.getAttribute("children")
    def keepBranch(i):
        deleted[i] = False
        for c in children[i]:
            if deleted[c]:
                keepBranch(c)

    for i in tree.iteratorFromPixelsToRoot():
        if not deleted[i]:   
            keepBranch(i)

