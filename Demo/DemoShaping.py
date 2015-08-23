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
Created on 11 juin 2015

@author: perretb
'''

from HiPy.IO import *  # @UnusedWildImport
from HiPy.Hierarchies.ComponentTree import *  # @UnusedWildImport
from HiPy.Util.Histogram import *  # @UnusedWildImport
from HiPy.Processing.Attributes import *  # @UnusedWildImport
from HiPy.Processing.Shaping import *  # @UnusedWildImport
from HiPy.Structures import *  # @UnusedWildImport


def dummyDemoShaping():
    '''
    Demonstrates that performing a threshold on an increasing attribute value in the original tree space
    is equivalent to thresholding the level of the shape tree.
    '''
    print("Reading image...")
    image = readImage("../samples/blood1.png")
    image.adjacency = AdjacencyNdRegular.getAdjacency2d4(image.embedding.size)

    print("Building tree...")
    tree = constructComponentTree(image, ComponentTreeType.TreeOfShapes)
    addAttributeArea(tree)
    addAttributeChildren(tree)

    print("Building tree of tree on area")
    area = prepareForShaping(tree, tree.area)
    tree2 = constructComponentTree(area)

    print("Filtering tree of tree on area based on level")
    tree2.filterDirect(lambda _, x: tree2.level[x] < 5000)
    print("Shaping first tree based on previous filter result")
    res = tree.reconstructImage("level", shapingCriterion(tree2))

    print("Filtering first tree based on area")
    res2 = tree.reconstructImage("level", lambda x: tree.area[x] < 5000)

    if res.equals(res2):
        print("Hurray: its the same")
    else:
        print("Yeark: does not work")


def demoNonIncreasingFilter():
    '''
    Demonstrate how to select attribute maxima of high value using a shaping
    '''
    print("Reading image...")
    image = readImage("../samples/blood1.png")
    image.adjacency = AdjacencyNdRegular.getAdjacency2d4(image.embedding.size)

    print("Building tree...")
    tree = constructComponentTree(image, ComponentTreeType.TreeOfShapes)
    addAttributeCompactness(tree)

    # prefiltering
    tree.filterDirect(lambda _, x: (tree.area[x] < 200 or tree.area[x] > tree.nbPixels * 0.6))
    updateAttributeAfterFiltering(tree, "compactness")

    print("Building tree of tree on Compactness...")
    compactness = prepareForShaping(tree, tree.compactness)
    tree2 = constructComponentTree(compactness)
    addAttributeExtrema(tree2)

    print("Shaping first tree based on previous filter result...")
    tree2.filterDirect(lambda _, x: tree2.level[x] < 0.6 or tree2.extrema[x] == False)
    res = tree.reconstructImage("level", shapingCriterion(tree2))

    print("Saving result...")
    saveImage(res, "Results/Shaping - high maxima of compactness.png")


def main():
    demoNonIncreasingFilter()
    dummyDemoShaping()


if __name__ == '__main__':
    main()
