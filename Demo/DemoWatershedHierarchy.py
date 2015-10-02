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
from math import *  # @UnusedWildImport
from HiPy.IO import *  # @UnusedWildImport
from HiPy.Hierarchies.ComponentTree import *  # @UnusedWildImport
from HiPy.Hierarchies.WatershedHierarchy import *  # @UnusedWildImport
from HiPy.Processing.Attributes import *  # @UnusedWildImport
from HiPy.Util.Histogram import *  # @UnusedWildImport
from HiPy.Util.VMath import *  # @UnusedWildImport
from HiPy.Util.Color import convertRGBtoLAB
from HiPy.Structures import AdjacencyNdRegular, WeightedAdjacency


def demoWatershedHierarchies():
    print("reading image...")
    image = readImage('../samples/monsters.png')
    image = convertRGBtoLAB(image)

    print("constructing gradient graph...")
    adj4 = AdjacencyNdRegular.getAdjacency2d4(image.embedding.size)
    adjacency = image.adjacency = WeightedAdjacency.createAdjacency(adj4,
                                                                    lambda i, j: euclideanDistance(image[i], image[j]))

    print("constructing BPT...")
    bpt = constructAltitudeBPT(adjacency)

    print("drawing saliency BPT...")
    salBpt = drawSaliencyForVisualisation(bpt, image)
    saveImage(salBpt, "Results/BPT Saliency map.png")

    print("constructing Component Tree...")
    naryTree = bpt.simplifyTreeByAttribute("level")
    print("drawing saliency Comp Tree...")
    salCt = drawSaliencyForVisualisation(naryTree, image)
    saveImage(salCt, "Results/CompTree Saliency map.png")

    if salBpt.equals(salCt):
        print("Good news BPT and simplified tree have the same saliency !")
    else:
        print("Yearkk! BPT and simplified tree don't have the same saliency !")

    print("constructing watershed hierarchy by altitude...")
    wsh = transformAltitudeBPTtoWatershedHierarchy(bpt)
    print("drawing saliency watershed hierarchy by altitude...")
    salWSh = drawSaliencyForVisualisation(wsh, image)
    saveImage(salWSh, "Results/Watershed by Altitude Saliency map.png")

    print("constructing watershed hierarchy by area...")
    addAttributeArea(wsh)
    wsha = transformBPTtoAttributeHierarchy(wsh, "area")
    print("drawing saliency watershed hierarchy by area...")
    salWSha = drawSaliencyForVisualisation(wsha, image)
    saveImage(salWSha, "Results/Watershed by Area Saliency map.png")


def main():
    demoWatershedHierarchies()


if __name__ == '__main__':
    main()
