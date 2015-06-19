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
Created on 15 juin 2015

@author: perretb
'''
from math import * #@UnusedWildImport
from HiPy.IO import * #@UnusedWildImport
from HiPy.Hierarchies.ComponentTree import * #@UnusedWildImport
from HiPy.Hierarchies.PartitionHierarchy import * #@UnusedWildImport
from HiPy.Processing.Attributes import * #@UnusedWildImport
from HiPy.Util.Histogram import * #@UnusedWildImport
from HiPy.Util.Color import convertRGBtoLAB


def distance(image, i, j):
    v1 = image[i]
    v2 = image[j]
    return sqrt((v1[0]-v2[0])**2+(v1[1]-v2[1])**2+(v1[2]-v2[2])**2)

def drawSaliencyForVizu(tree,image):
    adj4 = Adjacency2d4(image.embedding.size)
    saliency=computeSaliencyMap(tree,adj4)
    sal=drawSaliencyMap(image.embedding.size,saliency)
    sal = rescaleGray(sal, 0, 1)
    sal = imageMap(sal, lambda x:x**(0.3333))
    sal = normalizeToByte(sal)
    return sal
    
    
def demoBPT():
    print("reading image...")
    image = readImage('../samples/monsters.png', False)
    image = convertRGBtoLAB(image)
    
    print("constructing gradient graph...")
    adj4 = Adjacency2d4(image.embedding.size)
    adjacency=image.adjacency = AdjacencyEdgeWeightedGraph.createAdjacency(adj4, lambda i,j: distance(image,i,j) )
    
    print("constructing BPT...")
    bpt = constructAltitudeBPT(adjacency)
    print("drawing saliency BPT...")
    salBpt = drawSaliencyForVizu(bpt,image)
    saveImage(salBpt, "Results/BPT Saliency map.png")
    
    print("constructing Component Tree...")
    comptTree=transformAltitudeBPTtoComponentTree(bpt)
    print("drawing saliency Comp Tree...")
    salCt = drawSaliencyForVizu(comptTree,image)
    saveImage(salCt, "Results/CompTree Saliency map.png")
    
    if salBpt.equals(salCt):
        print("Good news BPT and Compt Tree have the same saliency !")
    else:
        print("Yearkk! BPT and Compt Tree don't have the same saliency !")
    
    print("constructing watershed hierarchy by altitude...")
    wsh=transformAltitudeBPTtoWatershedHierarchy(bpt)
    print("drawing saliency watershed hierarchy by altitude...")
    salWSh = drawSaliencyForVizu(wsh,image)
    saveImage(salWSh, "Results/Watershed by Altitude Saliency map.png")

    print("constructing watershed hierarchy by area...")
    addAttributeArea(wsh)
    wsha= transformBPTtoAttributeHierarchy(wsh,"area")
    print("drawing saliency watershed hierarchy by area...")
    salWSha = drawSaliencyForVizu(wsha,image)
    saveImage(salWSha, "Results/Watershed by Area Saliency map.png")
   

def main():
    demoBPT()
    
if __name__ == '__main__':
    main()
