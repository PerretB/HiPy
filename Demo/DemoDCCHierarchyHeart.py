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
from HiPy.Hierarchies.ComponentTree import * #@UnusedWildImport
from HiPy.Hierarchies.DirectedComponentHierarchy import * #@UnusedWildImport


# Create an edge weighted symmetric directed graph.
# The graph nodes are thee image pixels.
# Arcs are symmetric and constructed using a 4 or 8 connectivity.
# Weights are computed using the function "weightFunction" which is not necessarily symmetric.
# Weight function must a be a function that associates a dissimilarity measure to two grey level values 
# (the grey level value of the origin of the arc and the one of its destination)

def createHeartGradientGraph(image, size, weightFunction, connectivity=8):
    graph =  AdjacencyEdgeWeightedGraph(size[0]*size[1])
    width=size[0]
    height=size[1]
    coordLinTo2D = lambda x: (x % width, x // width)
    coord2DToLin= lambda x,y: y*width+x
    dim = width * height
    
    def writeLink(i,j):
        graph.createEdge(i,j,weightFunction(image[i],image[j]))
        
    for i in range(dim):
        x,y = coordLinTo2D(i)
        if(x+1<width):
            writeLink(i,coord2DToLin(x+1,y))
        if(x-1>=0):
            writeLink(i,coord2DToLin(x-1,y))
        if(y+1<height):
            writeLink(i,coord2DToLin(x,y+1))
        if(y-1>=0):
            writeLink(i,coord2DToLin(x,y-1))
        if(connectivity==8):
            if(x+1<width and y-1>=0):
                writeLink(i,coord2DToLin(x+1,y-1))
            if(x-1>=0 and y-1>=0):
                writeLink(i,coord2DToLin(x-1,y-1))
            if(x+1<width and y+1<height):
                writeLink(i,coord2DToLin(x+1,y+1))
            if(x-1>=0 and y+1<height):
                writeLink(i,coord2DToLin(x-1,y+1))        
    
    return graph
   




  


# Illustration on heart segmentation 
# Input: 
#    - file "heart.pgm": the original heart image (Fig. 10 a)
#    - file "heart-ObjectMarker.pgm": the object marker (Fig. 10 b)
#    - file "heart-BackgroundMarker.pgm": the background marker (Fig. 10 c)
# Output:
#    - file "Heart-Segmentation-DCC.png": the segmentation result (Fig. 10 f)
def testHeartSegmentation():
    print("Illustration on heart image segmentation")
    print("Reading images")
    image= readImage("../samples/DCC/heart.pgm")
    size=[image.embedding.width,image.embedding.height]
    markerObj= readImage("../samples/DCC/heart-ObjectMarker.pgm")
    markerFond = readImage("../samples/DCC/heart-BackgroundMarker.pgm")
    
    print("Creating adjacency")
    def weightFunction(v1,v2):
        diff=abs(v1-v2)
        if v2 <= 37 or v2 > 116: # hardcoded classification !  (Fig. 10 e)
            return (1.5)*diff
        else:
            return diff
    graph = createHeartGradientGraph(image, size, weightFunction, 8)
    
    print("Creating hierarchy")
    stack = createGraphStackFromEdgeWeightedGraph(graph)
    parent, completeGraphEdges, Lvls = directedComponentHierarchyStackFast(stack)
    dccTree = buildFinalDCCTree(stack.nbPoints, parent, completeGraphEdges, Lvls,image)
    
    print("Computing attributes")
    addAttributeMarker(dccTree, markerObj, 255, "marker1")  
    addAttributeMarker(dccTree, markerFond, 255, "marker2")
    addAttributeMarkerDirectedComponent(dccTree,"marker2")
    
    print("Filtering")
    def filterRule(tree, n):
        return not(tree.marker1[n] and not tree.marker2_directed[n])
    
    dccTree.filterDirect(filterRule)
    regularizeSelectMax(dccTree)
    regularizeSelectMaxTree(dccTree)

    print("Reconstruction")
    r1=dccTree.reconstructImage()
    saveImage(r1,"Heart-Segmentation-DCC.png")  
    print("=> Result saved in file " + " 'Heart-Segmentation-DCC.png'")
    
    print("Done\n\n")  
    
    
    
def main():
    testHeartSegmentation()   
    
if __name__ == '__main__':
    main()