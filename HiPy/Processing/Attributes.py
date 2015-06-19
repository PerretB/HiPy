# Copyright ESIEE (2015) 
# 
# benjamin.perret@esiee.fr
# 
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
Created on 19 juin 2015

@author: perretb
'''

import random
from math import * #@UnusedWildImport
import HiPy.Structures

def updateAttributeAfterFiltering(tree, attributeName, defaultValue=0):
    '''
    Update attribute values according to a filtering.
    
    The attribute value of each node marked as deleted is set to the attribute value of its closest non deleted ancestor.
    
    If no such ancestor exist, its attribute value is set to "defaultValue".
    '''
    deleted=tree.deleted
    attr = tree.getAttribute(attributeName)
    for i in tree.iteratorFromRootToPixels():
        if deleted[i]:
            par=tree[i]
            if par!=-1:
                attr[i]=attr[par]
            else:
                attr[i]=defaultValue
              

            
def addAttributeRandomColor(tree, name="randomColor"):
    attr=tree.addAttribute(name)
    if attr==None:
        return
    for i in tree.iteratorFromPixelsToRoot(False):
        attr[i]=(random.randint(0,255),random.randint(0,255),random.randint(0,255))

def addAttributeArea(tree, name="area"):
    attr=tree.addAttribute(name,0)
    if attr==None:
        return
    for i in tree.iteratorFromPixelsToRoot(True):
        if i<tree.nbPixels:
            attr[i]=1
        
        par=tree[i]
        if(par!=-1):
            attr[par]=attr[par]+attr[i]

           
def addAttributeChildren(tree, name="children"):
    attr=tree.addAttribute(name,[])
    if attr==None:
        return
    for i in tree.iteratorFromPixelsToRoot(True):
        par = tree[i]
        if(par!=-1):
            attr[par].append(i)
    

        
def addAttributeVolume(tree, name="volume"):
    addAttributeArea(tree)
    attr=tree.addAttribute(name,0)
    if attr==None:
        return
    for i in tree.iteratorFromPixelsToRoot(True):
        par = tree[i]
        if(par!=-1):
            lvl = tree.levels[i]
            lvlp = tree.levels[par]
            v=abs(lvl-lvlp)*tree.area[i]
            attr[i]=attr[i]+v
            attr[par]=attr[par]+v
            

def addAttributeMarker(tree, markerImage, markerValue, attributeName):
    '''
    Compute the attribute "attributeName" on the graph.
    For each node, the attribute is set to true if the  connected 
    component of the node contains a point of the image "markerImage" 
    having the value "markerValue"
    '''
    addAttributeChildren(tree)
    
    attr=tree.addAttribute(attributeName)
    if attr==None:
        return
   
    for i in range(tree.nbPixels):
        attr[i]=(markerImage[i]==markerValue)
    
    children=tree.getAttribute("children")
    for i in tree.iteratorFromPixelsToRoot(False):
        flag=False
        for j in children[i]:
            if attr[j]:
                flag = True
                break
        attr[i]=flag


def addAttributeSimpleMoments2d(tree, name="moments"):
    '''
    Compute the moments [M00 M10 M01 M11 M20 M02] of each component
    '''
    attr=tree.addAttribute(name)
    if attr==None:
        return
    nbLeaves=tree.nbPixels
    embedding=tree.leavesEmbedding
    
    for i in tree.iteratorFromPixelsToRoot():
        if i<nbLeaves:
            c=embedding.fromLinearCoordinate(i)
            x = c[0]
            y = c[1]
            attr[i]=[1, x, y, x * y, x * x, y * y]
        else:
            attr[i]=[0, 0, 0, 0, 0, 0]
    
    for i in tree.iteratorFromPixelsToRoot():
        par=tree[i]
        if(par != -1):
            m=attr[i]
            mp=attr[par]
            for j in range(len(m)):
                mp[j] += m[j]



def addAttributeInertia2d(tree,name="inertia",momentAttributeName="moments"):
    '''
    Compute the moment of inertia of each component 
    '''
    attr=tree.addAttribute(name)
    if attr==None:
        return
    attrMoment=tree.getAttribute(momentAttributeName)
   
    
    def computeInitertia(m):
        xmean = m[1] / m[0]
        ymean = m[2] / m[0]
        u20 = m[4]-xmean*m[1]
        u02 = m[5]-ymean*m[2]
        return (u20+u02)/(m[0]*m[0])
    
    for i in tree.iteratorFromPixelsToRoot():
        attr[i]=computeInitertia(attrMoment[i])
        

def addAttributeElongationOrientation2d(tree,nameElongation="elongation",nameOrientation="orientation",momentAttributeName="moments"):
    '''
    Compute the elongation and orientation of each component
    '''
    attrElongation=tree.addAttribute(nameElongation)
    attrOrientation=tree.addAttribute(nameOrientation)
    if attrElongation==None or attrOrientation==None:
        return
    attrMoment=tree.getAttribute(momentAttributeName)
    
    
    
    # Compute the elongation and the orientation from an array of moments [M00 M10 M01 M11 M20 M02]    
    def elongationOrientationFromMoments(m):
        sign = lambda x: copysign(1, x)
        xmean = m[1] / m[0]
        ymean = m[2] / m[0]
        xvar = m[4] / m[0] - xmean * xmean
        yvar = m[5] / m[0] - ymean * ymean
        xycovar = m[3] / m[0] - xmean * ymean
        if(xvar - yvar != 0):
            a = 0.5 * atan(2 * xycovar / (xvar - yvar))
            if(sign(a) * sign(xycovar) < 0.0):
                a += pi / 2.0;
            alpha = a
        else:
            if xycovar == 0:
                alpha = 0
            else:
                alpha = pi / 4.0 
        lambda1 = max(0,0.5 * (xvar + yvar + sqrt(4 * xycovar * xycovar + (xvar - yvar) * (xvar - yvar))))
        lambda2 = max(0,0.5 * (xvar + yvar - sqrt(4 * xycovar * xycovar + (xvar - yvar) * (xvar - yvar))))
    
        
        if ( lambda1 == 0): # ill posed case
            el = 1;
        elif(lambda2==0): # ill posed case
            el = sqrt((lambda2+1) / (lambda1+1));
        else :
            el = sqrt(lambda2 / lambda1);
    
        return el, alpha
    
    
    for i in tree.iteratorFromPixelsToRoot():
        attrElongation[i], attrOrientation[i] = elongationOrientationFromMoments(attrMoment[i])



def addAttributeDepth(tree, name="depth"):
    attr=tree.addAttribute(name,0)
    if attr==None:
        return
    for j in tree.iteratorFromRootToPixels():
        par=tree[j]
        if par!=-1:
            attr[j] = attr[par]+1

def addAttributeHighest(tree, name="highest"):
    attr=tree.addAttribute(name,None)
    if attr==None:
        return
    addAttributeChildren(tree)
    level = tree.level
    nbLeaves = tree.nbPixels
    children = tree.children
    for i in tree.iteratorFromPixelsToRoot():
        if i<nbLeaves:
            attr[i]=level[i]
        else:
            maxv=attr[children[i][0]]
            for c in children[i]:
                if(attr[c]>maxv):
                    maxv=attr[c]
            attr[i]=maxv
            
def addAttributeLowest(tree, name="lowest"):
    attr=tree.addAttribute(name,None)
    if attr==None:
        return
    addAttributeChildren(tree)
    level = tree.level
    nbLeaves = tree.nbPixels
    children = tree.children
    for i in tree.iteratorFromPixelsToRoot():
        if i<nbLeaves:
            attr[i]=level[i]
        else:
            minv=attr[children[i][0]]
            for c in children[i]:
                if(attr[c]<minv):
                    minv=attr[c]
            attr[i]=minv
            
def addAttributeDynamics(tree,extremaAttributeName="highest",name="dynamics"):
    attr=tree.addAttribute(name,None)
    if attr==None:
        return
    extrema = tree.getAttribute(extremaAttributeName)
    level = tree.level
    for i in tree.iteratorFromRootToPixels():
        parent = tree[i]
        if parent==-1:
            attr[i]=abs(extrema[i])
        else:
            if extrema[i]==extrema[parent]:
                attr[i]=attr[parent]
            else:
                attr[i]=abs(level[parent]-extrema[i])
        
def addAttributePerimeter(tree, name="perimeter", adjacency=None):
    attr=tree.addAttribute(name,0)
    if attr==None:
        return
    addAttributeChildren(tree)
    children = tree.getAttribute("children")
    if adjacency==None:
        adjacency=tree.leavesAdjacency
    nbLeaves = tree.nbPixels
    visited=HiPy.Structures.Image(nbLeaves,False)
    for i in range(nbLeaves):
        attr[i]=4#len(adjacency.getNeighbours(i)) FIXME !!!!
        
    for i in tree.iteratorFromPixelsToRoot(False):
        remove=0;
        for c in children[i]:
            attr[i] = attr[i] + attr[c]
            
            if c<nbLeaves:
                for n in adjacency.getNeighbours(c):
                    if visited[n]:
                        remove = remove + 2
                visited[c]=True
           
        attr[i] = attr[i] - remove
            
def addAttributeCompactness(tree, name="compactness"):
    attr=tree.addAttribute(name,0)
    if attr==None:
        return
    addAttributeArea(tree)
    addAttributePerimeter(tree)
    area=tree.getAttribute("area")
    perimeter=tree.getAttribute("perimeter")
    for i in tree.iteratorFromPixelsToRoot():
        attr[i]=4.0*pi*area[i]/(perimeter[i]*perimeter[i])  
    

        
def addAttributeComplexity(tree, name="complexity"):
    attr=tree.addAttribute(name,0)
    if attr==None:
        return
    addAttributeArea(tree)
    addAttributePerimeter(tree)
    area=tree.getAttribute("area")
    perimeter=tree.getAttribute("perimeter")
    for i in tree.iteratorFromPixelsToRoot():
        attr[i]=perimeter[i]/area[i] 
        
def addAttributeExtrema(tree, name="extrema"):
    ''' 
    true if node is a maxima, false otherwise
    '''
    attr=tree.addAttribute(name,True)
    if attr==None:
        return
    
    for i in tree.iteratorFromPixelsToRoot(False):
        par=tree[i]
        if par!=-1:
            attr[par]=False

def addAttributeIsLeaf(tree, name="isLeaf"):
    ''' 
    True if node is a logical leaves, false otherwise
    
    In a partition hierarchy, logical leaves and structure leaves are the same (pixels are part of the hierarchy)
    
    In a component tree, logical leaves are nodes such that every child is a structure leaf (pixels are not part of the hierarchy)
    '''
    attr=tree.addAttribute(name,True)
    if attr==None:
        return
    if tree.treeType==HiPy.Structures.TreeType.PartitionHierarchy:
        for i in tree.iteratorFromPixelsToRoot(False):
            attr[i]=False
    elif tree.treeType==HiPy.Structures.TreeType.ComponentTree:       
        for i in tree.iteratorFromPixelsToRoot(False):
            par=tree[i]
            if par!=-1:
                attr[par]=False
    else:
        raise Exception("addAttributeIsLeaf: unknown tree type.")