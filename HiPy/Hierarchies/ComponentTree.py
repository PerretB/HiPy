# Copyright ESIEE (2015)
#
# benjamin.perret@esiee.fr
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
This module is dedicated to the construction of Min or Max Trees of an image.

This module contains two main definitions:
    * the function 'constructComponentTree'.
    * the enum type ComponentTreeType contains two values 'MinTree' and 'MaxTree'.

Example
-------
The following example construct the Max Tree of an image:

    $ image = readImage('myImage.png')

    $ image.adjacency = AdjacencyNdRegular.getAdjacency2d4(image.embedding.size)

    $ tree = constructComponentTree(image,ComponentTreeType.MaxTree)

Please see :doc:`../Demo/DemoMaxTree.py` for a longer example.

Notes
-----
    The Min/Max Tree structure were proposed in [1]_, [2]_. The algorithm used in this
    implementation was first described in [3]_.

.. [1] Ph. Salembier, A. Oliveras, and L. Garrido, "Anti-extensive connected operators for image \
and sequence processing," IEEE Trans. Image Process., vol. 7, no. 4, pp. 555-570, Apr. 1998.
.. [2] Ro. Jones, "Connected filtering and segmentation using component trees," Comput. Vis. \
Image Understand., vol. 75, no. 3, pp. 215-228, Sep. 1999.
.. [3] Ch. Berger, T. Geraud, R. Levillain, N. Widynski, A. Baillard, and E. Bertin, "Effective \
Component Tree Computation with Application to Pattern Recognition in Astronomical Imaging," \
IEEE ICIP 2007.

Created on 3 juin 2015

@author: Benjamin Perret
'''

from HiPy.Structures import Image, Tree, TreeType, HiPyException, HiPyLogger
from HiPy.Processing.Attributes import addAttributeChildren
from HiPy.Util.UnionFind import findTarjan
from enum import Enum

class ComponentTreeType(Enum):
    '''
    This enum lists the supported component tree type: max and min trees.
    '''
    MinTree = 1
    MaxTree = 2


def constructComponentTree(image, treeType=ComponentTreeType.MaxTree):
    '''
    Construct the Min or Max Tree of the given image.

    Parameters
    ----------
    image : HiPy.Structures.Image
        The image to process.
        The image must be equipped with an adjacency (HiPy.Structures.Adjacency).
        The image values must be totally ordered.
    treeType : HiPy.Hierarchies.ComponentTree.ComponentTreeType (enum)
        The component tree can be either a Max Tree with ComponentTreeType.MaxTree (default)
        or a Min Tree with ComponentTreeType.MinTree

    Returns
    -------
    HiPy.Structures.Tree
        The component tree representing the input image.

    Raises
    ------
    HiPyException
        If treeType is different from ComponentTreeType.MaxTree and ComponentTreeType.MinTree.
    '''

    HiPyLogger.debug("call to: constructComponentTree")

    HiPyLogger.debug("Sorting values")
    if treeType == ComponentTreeType.MaxTree:
        sortedPixels = sorted(range(len(image)), key=lambda x: image[x])
    elif treeType == ComponentTreeType.MinTree:
        sortedPixels = sorted(range(len(image)), key=lambda x: -image[x])
    else:
        errorMsg = "Component Tree construction error : unknown tree type"
        HiPyLogger.error(errorMsg)
        raise HiPyException(errorMsg)

    parent = __preTreeConstruction(image, sortedPixels)

    __canonizeTree(parent, image, sortedPixels)

    parent, levels = __expandCanonizedParentRelation(parent, image, sortedPixels)

    return Tree(TreeType.ComponentTree, parent, levels, image)


def __preTreeConstruction(image, sortedPixels):
    '''
    Generic tree construction using ordered pixels and union find
    '''
    HiPyLogger.debug("call to: __preTreeConstruction")
    parent = image.copy(False)
    ufParent = Image(len(image), None)
    ufRank = Image(len(image))
    reprez = Image(len(image))
    for i in range(len(image)-1, -1, -1):
        currentPoint = sortedPixels[i]
        parent[currentPoint] = currentPoint
        ufParent[currentPoint] = currentPoint
        reprez[currentPoint] = currentPoint
        cpReprez = currentPoint
        for neighbour in image.getNeighbours(currentPoint):
            if ufParent[neighbour] != None:
                cpNeigh = findTarjan(neighbour, ufParent)
                if cpNeigh != cpReprez:
                    parent[reprez[cpNeigh]] = currentPoint
                    if ufRank[cpReprez] < ufRank[cpNeigh]:
                        cpNeigh, cpReprez = cpReprez, cpNeigh
                    ufParent[cpNeigh] = cpReprez
                    reprez[cpReprez] = currentPoint
                    if ufRank[cpReprez] == ufRank[cpNeigh]:
                        ufRank[cpReprez] = ufRank[cpReprez]+1

    return parent


def __canonizeTree(parent, enqueuedLevels, sortedPixels):
    '''
    Parent relation "canonization" (path compression) after preTreeConstruction
    '''
    HiPyLogger.debug("call to: __canonizeTree")
    for pixel in sortedPixels:
        par = parent[pixel]
        if enqueuedLevels[parent[par]] == enqueuedLevels[par]:
            parent[pixel] = parent[par]


def __expandCanonizedParentRelation(canonizedTreeImage, nodeLevels, sortedPixels):
    '''
    Expand a canonized parent relation to a regular parent relation
    (each node is represented individually)
    '''
    HiPyLogger.debug("call to: __expandCanonizedParentRelation")
    levels = nodeLevels.copy(True)
    data = [None]*len(canonizedTreeImage)
    for j in range(len(canonizedTreeImage)-1, -1, -1):
        i = sortedPixels[j]
        if nodeLevels[i] != nodeLevels[canonizedTreeImage[i]]:
            parent = i
        else:
            parent = canonizedTreeImage[i]
        if data[parent] == None:
            data.append(-1)
            data[parent] = len(data)-1
            levels.append(nodeLevels[parent])
        data[i] = data[parent]

    for j in range(len(canonizedTreeImage)-1, -1, -1):
        i = sortedPixels[j]
        if nodeLevels[i] != nodeLevels[canonizedTreeImage[i]]:
            parent = i
            pparent = canonizedTreeImage[parent]
            data[data[parent]] = data[pparent]
    #data[-1]=-1
    return data, levels


def regularizeSelectMaxTree(tree):
    '''
    Regularize the result of a selection criterion using the select max hierarchical strategy.

    A node is mark deleted if and only if it and all his ancestors were marked deleted.

    Parameters
    ----------
    tree : HiPy.Structures.Tree
        The tree containing the selection to be regularized

    Returns
    -------
    void

    '''

    # TODO: refactor location
    deleted = tree.getAttribute("deleted")
    addAttributeChildren(tree)
    children = tree.getAttribute("children")
    def keepBranch(i):
        '''
        Mark the whole branch rooted in node i as selected
        (deleted[n]=False, for all node n in the branch)
        '''
        deleted[i] = False
        for child in children[i]:
            if deleted[child]:
                keepBranch(child)

    for i in tree.iteratorFromPixelsToRoot():
        if not deleted[i]:
            keepBranch(i)

