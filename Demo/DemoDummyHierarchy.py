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
import random

from HiPy.Hierarchies.Quadtree import constructBinaryQuadTree
from HiPy.Hierarchies.WatershedHierarchy import drawSaliencyForVisualisation, constructAltitudeBPT
from HiPy.IO import readImage, saveImage
from HiPy.Processing.Attributes import addAttributeRank
from HiPy.Structures import AdjacencyNdRegular, WeightedAdjacency

__author__ = 'perretb'


def testQuadTree():
    print("reading image (but who cares for quadtre)...")
    image = readImage('../samples/monsters.png')

    print("constructing QuadTree...")
    bpt = constructBinaryQuadTree(image)

    print("drawing saliency BPT...")
    salBpt = drawSaliencyForVisualisation(bpt, image, gammaFactor=1)
    saveImage(salBpt, "Results/QuadTree Saliency map.png")


def testRandomHierarchy():
    rnd = random.Random()
    rnd.seed(1)
    print("reading image (but who cares for random hierarchy)...")
    image = readImage('../samples/monsters.png')
    adj4 = AdjacencyNdRegular.getAdjacency2d4(image.embedding.size)
    print("constructing random hierarchy...")
    adjacency = image.adjacency = WeightedAdjacency.createAdjacency(adj4, lambda i, j: rnd.random())
    bpt = constructAltitudeBPT(adjacency)
    print("drawing saliency BPT...")
    salBpt = drawSaliencyForVisualisation(bpt, image, gammaFactor=1)
    saveImage(salBpt, "Results/Random Saliency map.png")


def main():
    testQuadTree()
    testRandomHierarchy()

if __name__ == '__main__':
    main()
