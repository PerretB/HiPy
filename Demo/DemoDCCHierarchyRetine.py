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
Created on 11 juin 2015

@author: perretb
'''
from HiPy.IO import * #@UnusedWildImport
from HiPy.Hierarchies.ComponentTree import * #@UnusedWildImport
from HiPy.Hierarchies.DirectedComponentHierarchy import * #@UnusedWildImport


# Create a k nearest neighbour graph
#  - image: input image
#  - size: dimensions [width, height] of image
#  - similarityFunc: a function that weights the similarity between two pixel values (the lower the more similar)
#  - k: number of neighbors to select
#  - maxD: nearest neighbors are selected in a window of size (2*maxD+1)x(2*maxD+1)
def createKNNNeighbourhoodGraph(image,size,similarityFunc,k,maxD):
    graph =  DirectedWeightedAdjacency(size[0]*size[1])
    width=size[0]
    height=size[1]
    coordLinTo2D = lambda x: (x % width, x // width)
    coord2DToLin= lambda x,y: y*width+x
    dim = width * height
    for i in range(dim):
        px,py = coordLinTo2D(i)
        measures=[]
        for y in range(py-maxD,py+maxD+1,1):
            if y>=0 and y<height:
                for x in range(px-maxD,px+maxD+1,1):
                    if x>=0 and x<width:
                        if(x!=px or y!=py):
                            pos=coord2DToLin(x,y)
                            op=similarityFunc(image,size,i,pos)
                            if(op!=None):
                                measures.append((pos,op))
        measures=sorted(measures, key = lambda x: x[1])
        for j in range(min(len(measures),k)):
            graph.createEdge(i,measures[j][0])
    return graph
            

# Create a 4 or 8-neighborhood graph for an image of dimensions "size"
# k=4 or 8 indicates the adjacency.
# If "graph" equals None a new graph is created, otherwise the edges are added to the existing graph.
def createSimpleNeighbourhoodGraph(size,k=4, graph=None):
    if graph==None:
        graph =  DirectedWeightedAdjacency(size[0]*size[1])
    width=size[0]
    height=size[1]
    coordLinTo2D = lambda x: (x % width, x // width)
    coord2DToLin= lambda x,y: y*width+x
    dim = width * height

    def writeLink(i,j):
        if not j in graph.getSuccesors(i):
            graph.createEdge(i,j)
         
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
        if(k==8):
            if(x+1<width and y-1>=0):
                writeLink(i,coord2DToLin(x+1,y-1))
            if(x-1>=0 and y-1>=0):
                writeLink(i,coord2DToLin(x-1,y-1))
            if(x+1<width and y+1<height):
                writeLink(i,coord2DToLin(x+1,y+1))
            if(x-1>=0 and y+1<height):
                writeLink(i,coord2DToLin(x-1,y+1))        
    
    return graph

# Pseudo symilarity measure for finding the brighest pixels
#    - im: image
#    - size: dimensions of im
#    - p1: first pixel
#    - p2: second pixel
def simMax(im,size,p1,p2):
        return -im[p2]
    

# Creates the d-component hierarchy induced by the image and the adjacency. 
# Size is the dimensions [width,height] of the image.
def prepareHierarchy(adjacency,image,size):
    print("Creating hierarchy")
    stack = createGraphStackFromVertexWeightedGraph(adjacency,image)
    parent, completeGraphEdges, Lvls = directedComponentHierarchyStackFast(stack)
    dccTree = buildFinalDCCTree(stack.nbPoints, parent, completeGraphEdges, Lvls,image)
    
    print("Computing attributes")
    #computeAttributeLeavesCoordinates(graph, size[0])
    addAttributeAreaDirectedComponent(dccTree)
    addAttributeSimpleMomentsDirectedComponent(dccTree)
    addAttributeElongationOrientationDirectedComponent(dccTree)
    addAttributeInertiaDirectedComponent(dccTree)
    return dccTree

# Retinal image filtering of one image of the DRIVE database
#    imageNumber: number of the image to process
#    method: method number between 1 and 5 (see description of the method testRetineAll)
# Input: 
#    - file "DRIVE/<imageNumber>_test_filtered.png": the preprocessed image of the DRIVE database number <imageNumber>
# Output:
#    - filtering result saved in file "DRIVE/<imageNumber>_test_filtered_result_method_<method>.png"
def filterRetine(imageNumber,method=1):

    if method==1:
        print("Filtering with regularization on D-components on a non-local asymmetric graph.")
    elif method==2:
        print("Filtering on D-components on a non-local asymmetric graph.")
    elif method==3:
        print("Filtering on S-components on a non-local asymmetric graph.")
    elif method==4:
        print("Filtering on connected components on a non-local symmetric graph (symmetrization strategy: max).")
    elif method==5:
        print("Filtering on connected components on a non-local symmetric graph (symmetrization strategy: min).")
    else :
        print("Incorrect method !")
        return
    
    numS="%02d" % (imageNumber)
    fileName=numS+"_test_filtered"
    filePath="../samples/DCC/DRIVE/" +fileName + ".png"
    
    print("Reading image: " + filePath)
    image= readImage(filePath)
    size=[image.embedding.width,image.embedding.height]
    print("Creating adjacency")
    
    adj = createKNNNeighbourhoodGraph(image,size,simMax,3,4)
    adj = adj.getTranspose()
   
    if method==4:
        adj=adj.getSymmetricMax()
    elif method==5:
        adj=adj.getSymmetricMin()
    
    createSimpleNeighbourhoodGraph(size,4, adj)
    
    dccTree=prepareHierarchy(adj,image,size)


    print("Node selection")
    def filterRule_NL_Asym_DCC(dccTree, n):       
        return    dccTree.area[n]<0 or dccTree.area_directed[n]<15 or  (dccTree.inertia_directed[n]<1 and dccTree.area_directed[n]>=1000)  or (  dccTree.elongation_directed[n]>0.2 and dccTree.area_directed[n]<1000)
    
    def filterRule_NL_Asym_SCC(dccTree, n):       
        return    dccTree.area[n]<30 or  (dccTree.inertia[n]<1.3 and dccTree.area[n]>=1000)  or (  dccTree.elongation[n]>0.15 and dccTree.area[n]<1000)
    
    def filterRule_NL_SymMax(dccTree, n):       
        return    dccTree.area[n]<0 or dccTree.area_directed[n]<20 or  (dccTree.inertia_directed[n]<1.3 and dccTree.area_directed[n]>=1000)  or (  dccTree.elongation_directed[n]>0.2 and dccTree.area_directed[n]<1000)
   
    def filterRule_NL_SymMin(dccTree, n):       
        return    dccTree.area[n]<0 or dccTree.area_directed[n]<15 or  (dccTree.inertia_directed[n]<1.35 and dccTree.area_directed[n]>=1000)  or (  dccTree.elongation_directed[n]>0.15 and dccTree.area_directed[n]<1000)
   
    if method==1:
        dccTree.filterDirect(filterRule_NL_Asym_DCC)
        regularizeDiscardMin(dccTree)
    elif method==2:
        dccTree.filterDirect(filterRule_NL_Asym_DCC)
    elif method==3:
        dccTree.filterDirect(filterRule_NL_Asym_SCC)
    elif method==4:
        dccTree.filterDirect(filterRule_NL_SymMax)
    elif method==5:
        dccTree.filterDirect(filterRule_NL_SymMin)
        
    print("Reconstruction")

  
    saveName=fileName+"_result_method_" + str(method)+".png"  
    r1=dccTree.reconstructImage()
    saveImage(r1,"Results/" +saveName)

    
    print("=> Result saved in file " + saveName)
  

# Retinal image filtering of all the images of the DRIVE database
#    method: method number between 1 and 5 (see description in the paper)
#       - method==1: Filtering with regularization on D-components on a non-local asymmetric graph
#       - method==2: Filtering on D-components on a non-local asymmetric graph
#       - method==3: Filtering on S-components on a non-local asymmetric graph
#       - method==4: Filtering on connected components on a non-local symmetric graph (symmetrization strategy: max)
#       - method==5: Filtering on connected components on a non-local symmetric graph (symmetrization strategy: min)
# Input: 
#    - files "DRIVE/XX_test_filtered.png" for XX from 00 to 20: the preprocessed image of the DRIVE database number <imageNumber>
# Output:
#    - filtering result saved in file "DRIVE/XX_test_filtered_result_method_<method>.png"  for XX from 00 to 20
def testRetineAll(method):
    for i in range(1,21):
        filterRetine(i, method)  
        
        
        
    
def main():
    for i in range(1,6):
        testRetineAll(i) 
    
if __name__ == '__main__':
    main()