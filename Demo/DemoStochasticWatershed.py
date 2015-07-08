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
from HiPy.Util.Spatial import meanGray, spatialFilter
from HiPy.Util.Accumulator import BasicAccumulator

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
Created on 26 juin 2015

@author: perretb
'''

from math import * #@UnusedWildImport
from HiPy.IO import * #@UnusedWildImport
from HiPy.Hierarchies.ComponentTree import * #@UnusedWildImport
from HiPy.Hierarchies.PartitionHierarchy import * #@UnusedWildImport
from HiPy.Hierarchies.StochasticWatershed import * #@UnusedWildImport
from HiPy.Processing.Attributes import * #@UnusedWildImport
from HiPy.Util.Histogram import * #@UnusedWildImport
from HiPy.Util.VMath import * #@UnusedWildImport
from HiPy.Util.Color import convertRGBtoLAB
from HiPy.Structures import * #@UnusedWildImport


# determinist random generator
rnd = random.Random() 
rnd.seed(1)


def demoSWS():
    print("reading image...")
    image = readImage('../samples/monsters.png', False)
    image = convertRGBtoLAB(image)
    
    print("constructing gradient graph...")
    adj4 = AdjacencyNdRegular.getAdjacency2d4(image.embedding.size)
    adjacency=image.adjacency = WeightedAdjacency.createAdjacency(adj4, lambda i,j: euclideanDistance(image[i], image[j]))
    
    print("constructing stochastic watershed...")
    sws = constructExactRandomSeedsWatershed(adjacency, verbose=True)
    print("drawing saliency...")
    salSWS = drawSaliencyForVizu(sws,image)
    saveImage(salSWS, "Results/Random seeds SWS Saliency map.png")
    
    print("constructing area watershed...")
    wsh = transformAltitudeBPTtoWatershedHierarchy(constructAltitudeBPT(adjacency))
    addAttributeArea(wsh)
    wsha= transformBPTtoAttributeHierarchy(wsh,"area")
    print("Testing isomorphism...")
    print(Tree.testTreeIsomorphism(sws, wsha))
    
def demoNoiseWS(iteration=11):
    print("reading image...")
    image = readImage('../samples/lenna512.png', False)
    adj4 = AdjacencyNdRegular.getAdjacency2d4(image.embedding.size)
    image = rescaleGray(image,0,1)

    sal0=None
    for i in range(iteration):
        print("-> Iteration " + str(i))
        print("Adding noise...")
        im2 = imageMap(image, lambda x: x+rnd.uniform(-0.001, 0.001), marginal=True)
        im2 = rescaleGray(im2,0,1)
        im2 = convertRGBtoLAB(im2)
        print("Constructing gradient graph...")
        adjacency= WeightedAdjacency.createAdjacency(adj4, lambda i,j: euclideanDistance(im2[i], im2[j]))
        print("Constructing area watershed...")
        wsh = transformAltitudeBPTtoWatershedHierarchy(constructAltitudeBPT(adjacency))
        addAttributeArea(wsh)
        wsha= transformBPTtoAttributeHierarchy(wsh,"area")
        sal=computeSaliencyMap(wsha, adj4)
        
        
        print("Averaging and Drawing saliency...")
        lineGraph=WeightedAdjacency.createLineGraph(sal)
        adjLineGraps=WeightedAdjacency.createReflexiveRelation(lineGraph)
        meanSal=spatialFilter(sal, adjLineGraps, BasicAccumulator.getMeanAccumulator())
        saveImage(drawSaliencyMap(image.embedding.size,meanSal) , "Results/Random noise gradient "+ str(i) +".png")
        if sal0==None:
            sal0=meanSal
            for i in range(len(sal0)):
                sal0[i]=[sal0[i]]
        else:
            for i in range(len(sal0)):
                sal0[i].append(meanSal[i])
    print("Merging results...")

    
    
    for i in range(len(sal0)):
        sal0[i]=sum(sal0[i])/len(sal0[i])
    
        
    saveImage(drawSaliencyMap(image.embedding.size,sal0) , "Results/Random noise combined gradient.png")
    
    print("Ultra metric opening...")
    bpt = transformAltitudeBPTtoWatershedHierarchy(constructAltitudeBPT(sal0))
    saveImage(drawSaliencyForVizu(bpt,image), "Results/Random noise WS bpt.png")
    
    print("Area watershed...")
    addAttributeArea(bpt)
    wsa= transformBPTtoAttributeHierarchy(bpt,"area")
    saveImage(drawSaliencyForVizu(wsa,image), "Results/Random noise WSH.png")
 
def testEdgeFilter():
    image = readImage('../samples/lennaGray256.png', False)
    adj4 = AdjacencyNdRegular.getAdjacency2d4(image.embedding.size)
    im2 = rescaleGray(image,0,1)
    
    print("Constructing gradient graph...")
    adjacency= WeightedAdjacency.createAdjacency(adj4, lambda i,j: abs(im2[i] - im2[j]))
    print("Constructing area watershed...")
    wsh = transformAltitudeBPTtoWatershedHierarchy(constructAltitudeBPT(adjacency))
    addAttributeArea(wsh)
    wsha= transformBPTtoAttributeHierarchy(wsh,"area")
    sal=computeSaliencyMap(wsha, adj4)
    
    lineGraph=WeightedAdjacency.createLineGraph(sal)
    adjLineGraps=WeightedAdjacency.createReflexiveRelation(lineGraph)
    meanSal=spatialFilter(sal, adjLineGraps, BasicAccumulator.getMeanAccumulator())
    
    imSal1=drawSaliencyMap(image.embedding.size,sal) 
    imSal2=drawSaliencyMap(image.embedding.size,meanSal) 
    
    saveImage(imSal1, "./Result/lenna area sal map.png")
    saveImage(imSal2, "./Result/lenna area mean sal map.png")
    

    
    
def main():
    #demoSWS()
    demoNoiseWS()
    #testEdgeFilter()

    
    
if __name__ == '__main__':
    main()