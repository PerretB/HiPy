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
from HiPy.Hierarchies.WatershedHierarchy import constructAltitudeBPT, transformBPTtoAttributeHierarchy, \
    drawSaliencyForVisualisation
from HiPy.IO import readImage, saveImage
from HiPy.Processing.Attributes import addAttributeArea
from HiPy.Processing.EnergyHierarchy import transformTreeToOptimalMumfordShahEnergyCutHierarchy, \
    addAttributeMumfordShahEnergy, addAttributeEnergyCut, constructMumfordShahEnergyBPT
from HiPy.Structures import AdjacencyNdRegular, WeightedAdjacency
from HiPy.Util import VMath
from HiPy.Util.Color import convertRGBtoLAB

__author__ = 'Benjamin Perret'


def demoEnergyCut():
    print("Construction primal hierarchy by area...")
    image = readImage("../samples/monsters.png")
    image = convertRGBtoLAB(image)
    image.adjacency = AdjacencyNdRegular.getAdjacency2d4(image.embedding.size)
    weightedAdjacency = WeightedAdjacency.createAdjacency(image.adjacency,
                                                          lambda i, j: VMath.euclideanDistance(image[i], image[j]))
    tree = constructAltitudeBPT(weightedAdjacency)
    addAttributeArea(tree)
    tree = transformBPTtoAttributeHierarchy(tree, "area")

    print("Computing Mumford-Shah energy...")
    addAttributeMumfordShahEnergy(tree, 9000, image)

    print("Cutting...")
    addAttributeEnergyCut(tree, "mumfordShahEnergy")

    print("Drawing saliency")
    result = drawSaliencyForVisualisation(tree, image, "energyCut")

    fileName = "Results/Optimal energy cut.png"
    print("Saving result to: " + fileName)
    saveImage(result, fileName)
    print("Done")


def demoOptimalEnergyCutHierarchy():
    print("Construction primal hierarchy by area...")
    image = readImage("../samples/monsters.png")
    image = convertRGBtoLAB(image)
    image.adjacency = AdjacencyNdRegular.getAdjacency2d4(image.embedding.size)
    weightedAdjacency = WeightedAdjacency.createAdjacency(image.adjacency,
                                                          lambda i, j: VMath.euclideanDistance(image[i], image[j]))
    tree = constructAltitudeBPT(weightedAdjacency)
    addAttributeArea(tree)
    tree = transformBPTtoAttributeHierarchy(tree, "area")

    print("Computing optimal energy cut hierarchy w.r.t. Mumford-Shah energy...")
    energyTree = transformTreeToOptimalMumfordShahEnergyCutHierarchy(tree, image)

    print("Drawing saliency")
    result = drawSaliencyForVisualisation(energyTree, image, gammaFactor=1)

    fileName = "Results/Energy cut hierarchy Saliency map.png"
    print("Saving result to: " + fileName)
    saveImage(result, fileName)
    print("Done")


def demoOptimalEnergyHierarchy():
    print("reading image...")
    image = readImage("../samples/monsters.png")
    image = convertRGBtoLAB(image)
    adj4 = AdjacencyNdRegular.getAdjacency2d4(image.embedding.size)

    print("Computing optimal energy hierarchy w.r.t. Mumford-Shah energy...")
    bpt = constructMumfordShahEnergyBPT(image, adj4)

    print("Drawing saliency")
    result = drawSaliencyForVisualisation(bpt, image, gammaFactor=1)

    fileName = "Results/Energy hierarchy Saliency map.png"
    print("Saving result to: " + fileName)
    saveImage(result, fileName)
    print("Done")


def main():
    print("--- Optimal energy cut")
    demoEnergyCut()

    print("--- Optimal energy cut hierarchy")
    demoOptimalEnergyCutHierarchy()

    print("--- Optimal energy hierarchy")
    demoOptimalEnergyHierarchy()

if __name__ == '__main__':
    main()