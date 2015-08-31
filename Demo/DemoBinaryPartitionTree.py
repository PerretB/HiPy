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

from HiPy.Hierarchies.BinaryPartitionTree import constructBPT, singleLinkageFusionWeights, completeLinkageFusionWeights
from HiPy.Hierarchies.WatershedHierarchy import drawSaliencyForVisualisation, constructAltitudeBPT
from HiPy.IO import readImage, saveImage
from HiPy.Structures import AdjacencyNdRegular, WeightedAdjacency
from HiPy.Util.Color import convertRGBtoLAB
from HiPy.Util.VMath import euclideanDistance

__author__ = 'Benjamin Perret'


def testBPTSingleLinkage():
    print("reading image...")
    image = readImage("../samples/remotesensing1.png", False)# "../samples/monsters.png", False)#
    image = convertRGBtoLAB(image)

    print("constructing gradient graph...")
    adj4 = AdjacencyNdRegular.getAdjacency2d4(image.embedding.size)
    adjacency = image.adjacency = WeightedAdjacency.createAdjacency(adj4,
                                                                    lambda i, j: euclideanDistance(image[i], image[j]))

    print("constructing single linkage BPT...")
    bpt = constructBPT(image, adjacency, singleLinkageFusionWeights)

    print("constructing altitude BPT...")
    bpt2 = constructAltitudeBPT(adjacency)

    print("drawing saliencies...")
    salBpt = drawSaliencyForVisualisation(bpt, image)
    salBpt2 = drawSaliencyForVisualisation(bpt2, image)

    if salBpt.equals(salBpt2):
        print("Good news single linkage BPT and altitude BPT have the same saliency !")
    else:
        print("Yearkk!  single linkage BPT and altitude BPT don't have the same saliency !")


def testBPTCompleteLinkage():
    print("reading image...")
    image = readImage("../samples/monsters.png", False)# "../samples/remotesensing1.png", False)#"samples/monsters.png", False)#
    image = convertRGBtoLAB(image)

    print("constructing gradient graph...")
    adj4 = AdjacencyNdRegular.getAdjacency2d4(image.embedding.size)
    adjacency = image.adjacency = WeightedAdjacency.createAdjacency(adj4,
                                                                    lambda i, j: euclideanDistance(image[i], image[j]))

    print("constructing complete linkage BPT...")
    bpt = constructBPT(image, adjacency, completeLinkageFusionWeights)

    print("drawing saliency BPT...")
    salBpt = drawSaliencyForVisualisation(bpt, image)
    saveImage(salBpt, "Results/complete linkage BPT Saliency map.png")


def main():
    testBPTSingleLinkage()
    testBPTCompleteLinkage()

if __name__ == '__main__':
    main()
