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

from HiPy.Structures import Image, Tree
from HiPy.Hierarchies.ComponentTree import preTreeConstruction, canonizeTree, expandCanonizedParentRelation
from HiPy.Util.Geometry2d import interpolateX2, interpolatePlainMapKhalimsky

from collections import deque


# construction of the TOS
# assumptions: 
#  - (interpolated) pixel values are positive integers
#  - pixel 0 is  p(infiny)
def constructTreeOfShapes(image,interpolation=max,verbose=False):
    
    # interpolation 
    if interpolation!=None:
        image1=interpolateX2(image,interpolation)
    else: 
        image1=image
    #plain map
    if verbose:
        print("Plain Map Creation")
    image2=interpolatePlainMapKhalimsky(image1)
    # sorting
    if verbose:
        print("Sorting")
    (sortedPixels, enqueuedLevels)=sort(image2,max(image)+1)
    # generic tree construction
    if verbose:
        print("Union Find")
    parent=preTreeConstruction(image2,sortedPixels)
    # canonize tree
    if verbose:
        print("Canonize")
    canonizeTree(parent,enqueuedLevels,sortedPixels)
    parent,levels = expandCanonizedParentRelation(parent,enqueuedLevels,sortedPixels)
    # tree is left in interpolated space...
   
    # encapsulated tree for ease of manipulation
    if verbose:
        print("Tree finalization")
    return Tree(parent,levels,image2)

class PriorityQueue:
    def __init__(self, levels):
        self.data = []
        for _ in range(levels):
            self.data.append(deque())
        self.levels = levels
        self.size=0
    
    def push(self,level,element):
        self.data[level].append(element)
        self.size = self.size+1
        
    def isEmpty(self):
        return self.size==0
    
    def isLevelEmpty(self,level):
        return len(self.data[level])==0
    
    def pop(self,level):
        if len(self.data[level])>0:
            self.size=self.size-1
            return self.data[level].popleft()
        return None
    
    def findClosestNonEmpty(self,level):
        if not self.isLevelEmpty(level):
            return level
        lvlb=level-1
        lvlu=level+1
        while lvlb>=0 or lvlu<self.levels:
            if lvlb>=0:
                if not self.isLevelEmpty(lvlb):
                    return lvlb
                lvlb=lvlb-1
            if lvlu<self.levels:
                if not self.isLevelEmpty(lvlu):
                    return lvlu
                lvlu=lvlu+1
        return None




def priorityPush(queue, point, image, currentLevel):
    level = image[point]
    if isinstance(level, (list, tuple)):
        low = level[0]
        up  = level[1]
    else:
        low=level
        up=level
    
    newLevel = min(up,max(low,currentLevel))
    queue.push(newLevel,point)

def priorityPop(queue, currentLevel):
    newLevel = queue.findClosestNonEmpty(currentLevel)
    newPoint = queue.pop(newLevel)
    return (newLevel,newPoint)

def sort(image,levels):
    dejaVu = Image(len(image),False)
    enqueuedLevels = Image(len(image))
    sortedPixels=[]
    queue = PriorityQueue(levels)
    startPoint=0                               #  Uhuhuhuhuuuuuuuuuuuu !!!!!!!
    queue.push(image[startPoint][0],startPoint) 
    currentLevel = image[startPoint][0]
    dejaVu[startPoint]=True
    i=0
    while not queue.isEmpty():
        (currentLevel,currentPoint) = priorityPop(queue,currentLevel)
        enqueuedLevels[currentPoint] = currentLevel 
        i=i+1
        sortedPixels.append(currentPoint)
        for neighbour in image.getNeighbours(currentPoint):
            if not dejaVu[neighbour]:
                priorityPush(queue,neighbour,image,currentLevel)
                dejaVu[neighbour]=True
    return (sortedPixels, enqueuedLevels)