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

"""
Created on 10 juin 2015

@author: perretb
"""
from HiPy.IO import *  # @UnusedWildImport
from HiPy.Hierarchies.ComponentTree import *  # @UnusedWildImport
from HiPy.Hierarchies.DirectedComponentHierarchy import *  # @UnusedWildImport


# Create the graph for the neurite application.
# Image is the vesselness image.
def createSemanticGraphNeurite(image, size, tube=1, bulbe=2):
    graph = DirectedWeightedAdjacency(size[0] * size[1])
    width = size[0]
    height = size[1]
    coordLinTo2D = lambda x: (x % width, x // width)
    coord2DToLin = lambda x, y: y * width + x
    dim = width * height

    def writeLink(i, j):
        if image[i] == image[j]:
            graph.createEdge(i, j)
        elif image[i] == tube and image[j] == bulbe:
            graph.createEdge(i, j)
        elif image[i] == 0:
            graph.createEdge(i, j)

    for i in range(dim):
        x, y = coordLinTo2D(i)
        if x + 1 < width:
            writeLink(i, coord2DToLin(x + 1, y))
        if x - 1 >= 0:
            writeLink(i, coord2DToLin(x - 1, y))
        if y + 1 < height:
            writeLink(i, coord2DToLin(x, y + 1))
        if y - 1 >= 0:
            writeLink(i, coord2DToLin(x, y - 1))
        if x + 1 < width and y - 1 >= 0:
            writeLink(i, coord2DToLin(x + 1, y - 1))
        if x - 1 >= 0 and y - 1 >= 0:
            writeLink(i, coord2DToLin(x - 1, y - 1))
        if x + 1 < width and y + 1 < height:
            writeLink(i, coord2DToLin(x + 1, y + 1))
        if x - 1 >= 0 and y + 1 < height:
            writeLink(i, coord2DToLin(x - 1, y + 1))

    return graph


# Illustration on neurite image filtering
# Input: 
#    - file "neurite.png": the original neurite image (Fig. 1 a)
#    - file "neurite-vesselness.png": the classification into background, blobs and vessels (Fig. 9 a)
# Output:
#    - file "Neurite-Filtered-DCC.png": the filtered image (Fig. 9 b)
#    - file "Neurite-Filtered-DCC-Regularized.png": the regularized filtered image (Fig. 1 b)
def testNeuriteFiltering():
    print("\nIllustration on neurite image filtering")
    print("Reading images")
    classeFile = "../samples/DCC/neurite-vesselness.png"
    greyFile = "../samples/DCC/neurite.png"

    print("Creating adjacency")
    LabelV = 70  # gray level of vessels in the classification
    LabelB = 140  # gray level of blobs in the classification

    image = readImage(greyFile)
    size = [image.embedding.width, image.embedding.height]
    classif = readImage(classeFile)
    adj = createSemanticGraphNeurite(classif, size, LabelV, LabelB)

    print("Creating hierarchy")
    stack = createGraphStackFromVertexWeightedGraph(adj, image)
    # here we don't use the optimized version of the algorithm 
    # because it does no allow us to easily compute the attribute nbOut
    parent, DAGsAdj, Lvls, nbOut = directedComponentHierarchyStack(stack)
    ensureTree(parent, Lvls)
    dccTree = buildFinalDCCTree(stack.nbPoints, parent, DAGsAdj, Lvls, image, nbOut)

    print("Filtering")

    def filterRule(tree, n):
        return len(tree.sucs[n]) < 2 or tree.nbOut[n] > 20

    dccTree.filterDirect(filterRule)
    r1 = dccTree.reconstructImage()
    saveImage(r1, "Results/Neurite-Filtered-DCC.png")

    print("=> Result saved in file " + " 'Neurite-Filtered-DCC.png'")

    regularizeSelectMax(dccTree)
    regularizeSelectMaxTree(dccTree)
    r2 = dccTree.reconstructImage()
    saveImage(r2, "Results/Neurite-Filtered-DCC-Regularized.png")

    print("=> Result saved in file " + " 'Neurite-Filtered-DCC-Regularized.png'")

    print("Done\n\n")


def main():
    testNeuriteFiltering()


if __name__ == '__main__':
    main()
