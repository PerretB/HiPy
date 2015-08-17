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

"""
This module is dedicated to the construction of Component Trees:
i.e. Min Tree, Max Tree, or Tree of Shapes of an image.

This module contains two main definitions:
    * the function 'constructComponentTree'.
    * the enum type ComponentTreeType contains two values 'MinTree' and 'MaxTree'.

Example
-------
The following example construct the Max Tree of an image:

    $ image = readImage('myImage.png')

    $ image.adjacency = AdjacencyNdRegular.getAdjacency2d4(image.embedding.size)

    $ tree = constructComponentTree(image,ComponentTreeType.MaxTree)

More examples in:
    * doc:`../Demo/DemoMaxTree.py`
    * doc:`../Demo/DemoTreeOfShapes.py`
    * doc:`../Demo/DemoShaping.py`

Notes
-----
    The Min/Max Tree structure were proposed in [1]_, [2]_. The algorithm used in this
    implementation was first described in [3]_.

    The Tree of Shapes was described in [4]_. The algorithm used in this
    implementation was first described in [5]_.

.. [1] Ph. Salembier, A. Oliveras, and L. Garrido, "Anti-extensive connected operators for image \
and sequence processing," IEEE Trans. Image Process., vol. 7, no. 4, pp. 555-570, Apr. 1998.
.. [2] Ro. Jones, "Connected filtering and segmentation using component trees," Comput. Vis. \
Image Understand., vol. 75, no. 3, pp. 215-228, Sep. 1999.
.. [3] Ch. Berger, T. Geraud, R. Levillain, N. Widynski, A. Baillard, and E. Bertin, "Effective \
Component Tree Computation with Application to Pattern Recognition in Astronomical Imaging," \
IEEE ICIP 2007.
.. [4] Pa. Monasse, and F. Guichard, "Fast computation of a contrast-invariant image representation," \
Image Processing, IEEE Transactions on , vol.9, no.5, pp.860,872, May 2000
.. [5] Th. Géraud, E. Carlinet, S. Crozet, and L. Najman, "A Quasi-linear Algorithm to Compute the Tree \
of Shapes of n D Images", ISMM 2013.

Created on 3 juin 2015

@author: Benjamin Perret
"""

from HiPy.Structures import Image, Tree, TreeType, HiPyException, HiPyLogger, PriorityQueue
from HiPy.Processing.Attributes import addAttributeChildren
from HiPy.Util.Geometry2d import interpolatePlainMapKhalimsky
from HiPy.Util.UnionFind import findTarjan
from enum import Enum


class ComponentTreeType(Enum):
    """
    This enum lists the supported component tree type: min tree, max tree, and tree of shapes.
    """
    MinTree = 1
    MaxTree = 2
    TreeOfShapes = 3


def constructComponentTree(image: Image, treeType: ComponentTreeType=ComponentTreeType.MaxTree) -> Tree:
    """
    Construct the Min Tree, Max Tree or Tree of shapes of the given image.

    This function is a simple wrapper around the three base function:
        * constructMinTree
        * constructMaxTree
        * constructTreeOfShapes
    Please read the documentation of the specific function, for details on the constraints and the output for each
    type of tree.

    :param image: The image to process.
    :param treeType : The type of tree to construct

    :returns HiPy.Structures.Tree
        The component tree representing the input image.

    Raises
    ------
    HiPyException
        If treeType is not a valid element of the ComponentTreeType enum.
    """

    HiPyLogger.debug("call to: constructComponentTree")

    HiPyLogger.debug("Sorting values")
    if treeType == ComponentTreeType.MaxTree:
        return constructMaxTree(image)
    elif treeType == ComponentTreeType.MinTree:
        return constructMinTree(image)
    elif treeType == ComponentTreeType.TreeOfShapes:
        return constructTreeOfShapes(image)
    else:
        errorMsg = "Component Tree construction error : unknown tree type"
        HiPyLogger.error(errorMsg)
        raise HiPyException(errorMsg)


def constructMaxTree(image: Image) -> Tree:
    """
    Construct the Max Tree of the given image.

    :param image: image: HiPy.Structures.Image
        * The image must be equipped with an adjacency (HiPy.Structures.Adjacency).
        * The image values must be totally ordered.

    :return: the Max Tree of the input image.
    """
    HiPyLogger.debug("call to: constructMaxTree")
    sortedPixels = sorted(range(len(image)), key=lambda x: image[x])
    return __constructTreeFromSortedPixels(image, sortedPixels)


def constructMinTree(image: Image) -> Tree:
    """
    Construct the Min Tree of the given image.

    :param image: image: HiPy.Structures.Image
        * The image must be equipped with an adjacency (HiPy.Structures.Adjacency).
        * The image values must be totally ordered.

    :return: the Min Tree of the input image.
    """
    HiPyLogger.debug("call to: constructMaxTree")
    sortedPixels = sorted(range(len(image)), key=lambda x: -image[x])
    return __constructTreeFromSortedPixels(image, sortedPixels)


def constructTreeOfShapes(image: Image) -> Tree:
    """
    Construct the tree of shapes of the given image.

    :param image: image: HiPy.Structures.Image
        * The image should be equipped with a 2d grid embedding.
        * The image values must be (small) positive integers.
        * The image is interpolated in the Khalimsky grid in order to compute the tree and the functions returns \
    the tree of the interpolated image (resolution doubled).
        * The pixel of index 0 is considered to be the pixel "at infinity" i.e. the exterior of everything.

    :return: the tree of shapes of the input image.

    TODO: Add option to leave interpolated space
    TODO: Add option to change origin pixel
    TODO: Add option to disable Khalimsky interpolation
    """

    HiPyLogger.debug("call to: constructTreeOfShapes")

    image2 = interpolatePlainMapKhalimsky(image)

    HiPyLogger.debug("Sorting")
    (sortedPixels, enqueueLevels) = sortTreeOfShape(image2, max(image) + 1)
    return __constructTreeFromSortedPixels(enqueueLevels, sortedPixels)


def __constructTreeFromSortedPixels(image: Image, sortedPixels: list) -> Tree:
    parent = __preTreeConstruction(image, sortedPixels)
    __canonizeTree(parent, image, sortedPixels)
    parent, levels = __expandCanonizedParentRelation(parent, image, sortedPixels)
    return Tree(TreeType.ComponentTree, parent, levels, image)


def __preTreeConstruction(image: Image, sortedPixels: list) -> Image:
    """
    Generic tree construction using ordered pixels and union find
    :param image: image to process
    :param sortedPixels: the sorted pixel indices

    :returns a pre-parent relation
    """
    HiPyLogger.debug("call to: __preTreeConstruction")
    parent = image.copy(False)
    ufParent = Image(len(image), None)
    ufRank = Image(len(image))
    representing = Image(len(image))
    for i in range(len(image) - 1, -1, -1):
        currentPoint = sortedPixels[i]
        parent[currentPoint] = currentPoint
        ufParent[currentPoint] = currentPoint
        representing[currentPoint] = currentPoint
        cpReprez = currentPoint
        for neighbour in image.getNeighbours(currentPoint):
            if ufParent[neighbour] is not None:
                cpNeigh = findTarjan(neighbour, ufParent)
                if cpNeigh != cpReprez:
                    parent[representing[cpNeigh]] = currentPoint
                    if ufRank[cpReprez] < ufRank[cpNeigh]:
                        cpNeigh, cpReprez = cpReprez, cpNeigh
                    ufParent[cpNeigh] = cpReprez
                    representing[cpReprez] = currentPoint
                    if ufRank[cpReprez] == ufRank[cpNeigh]:
                        ufRank[cpReprez] += 1

    return parent


def __canonizeTree(parent: list, enqueueLevels: list, sortedPixels: list) -> None:
    """
    Parent relation "canonization" (path compression) after preTreeConstruction
    :param parent: a pre-parent relation as constructed by preTreeConstruction
    :param enqueueLevels: the node levels associated to the pre-parent relation
    :param sortedPixels: the sorted pixel indices
    """
    HiPyLogger.debug("call to: __canonizeTree")
    for pixel in sortedPixels:
        par = parent[pixel]
        if enqueueLevels[parent[par]] == enqueueLevels[par]:
            parent[pixel] = parent[par]


def __expandCanonizedParentRelation(canonizedTreeImage: Image, nodeLevels: Image, sortedPixels: list) -> (list, list):
    """
    Expand a canonized parent relation to a regular parent relation
    (each node is represented individually)
    :param canonizedTreeImage: a canonized parent relation
    :param nodeLevels: the node levels associated to the canonized parent relation
    :param sortedPixels: the sorted pixel indices

    :returns
        * parent relation
        * node levels associated to the parent relation
    """
    HiPyLogger.debug("call to: __expandCanonizedParentRelation")
    levels = nodeLevels.copy(True)
    data = [None] * len(canonizedTreeImage)
    for j in range(len(canonizedTreeImage) - 1, -1, -1):
        i = sortedPixels[j]
        if nodeLevels[i] != nodeLevels[canonizedTreeImage[i]]:
            parent = i
        else:
            parent = canonizedTreeImage[i]
        if data[parent] is None:
            data.append(-1)
            data[parent] = len(data) - 1
            levels.append(nodeLevels[parent])
        data[i] = data[parent]

    for j in range(len(canonizedTreeImage) - 1, -1, -1):
        i = sortedPixels[j]
        if nodeLevels[i] != nodeLevels[canonizedTreeImage[i]]:
            parent = i
            pparent = canonizedTreeImage[parent]
            data[data[parent]] = data[pparent]
    # data[-1]=-1
    return data, levels


def priorityPush(queue: PriorityQueue, point: int, image: Image, currentLevel: int) -> None:
    """
    Finds the adequate level to insert a new pixel of the image in the given priority queue and inserts it.
    :param queue: a multi-level priority queue
    :param point: pixel index
    :param image: pixel values
    :param currentLevel: current propagation level
    :return: None
    """
    level = image[point]
    if isinstance(level, (list, tuple)):
        lowLevel = level[0]
        upLevel = level[1]
    else:
        lowLevel = level
        upLevel = level

    newLevel = min(upLevel, max(lowLevel, currentLevel))
    queue.push(newLevel, point)


def priorityPop(queue: PriorityQueue, currentLevel: int) -> (int, int):
    """
    Finds and pop the closest element from the current priority level in the given multi level priority queue
    :param queue: a multi-level priority queue
    :param currentLevel: current priority level
    :return: (level of the closest element, closest element)
    """
    newLevel = queue.findClosestNonEmpty(currentLevel)
    newPoint = queue.pop(newLevel)
    return newLevel, newPoint


def sortTreeOfShape(image: Image, maxLevel: int) -> (list, list):
    """
    Sort image pixels according to propagation through level lines.
    :param image: input image, should be a range image
    :param maxLevel: maximum value in the image
    :return: list of sorted pixel indices, list of enqueued levels
    """
    dejaVu = Image(len(image), False)
    enqueuedLevels = image.copy(False)
    sortedPixels = []
    queue = PriorityQueue(maxLevel)
    startPoint = 0  # Uhuhuhuhuuuuuuuuuuuu !!!!!!!
    queue.push(image[startPoint][0], startPoint)
    currentLevel = image[startPoint][0]
    dejaVu[startPoint] = True
    i = 0
    while not queue.isEmpty():
        (currentLevel, currentPoint) = priorityPop(queue, currentLevel)
        enqueuedLevels[currentPoint] = currentLevel
        i += 1
        sortedPixels.append(currentPoint)
        for neighbour in image.getNeighbours(currentPoint):
            if not dejaVu[neighbour]:
                priorityPush(queue, neighbour, image, currentLevel)
                dejaVu[neighbour] = True
    return sortedPixels, enqueuedLevels


def regularizeSelectMaxTree(tree: Tree) -> None:
    """
    Regularize the result of a selection criterion using the select max hierarchical strategy.

    A node is mark deleted if and only if it and all his ancestors were marked deleted.

    :param tree :  The tree containing the selection to be regularized

    :return void

    """

    # TODO: refactor location
    deleted = tree.getAttribute("deleted")
    addAttributeChildren(tree)
    children = tree.getAttribute("children")

    def keepBranch(j):
        """
        Mark the whole branch rooted in node i as selected
        (deleted[n]=False, for all node n in the branch)

        :param j: root node i of the branch to select
        """
        deleted[j] = False
        for child in children[j]:
            if deleted[child]:
                keepBranch(child)

    for i in tree.iteratorFromPixelsToRoot():
        if not deleted[i]:
            keepBranch(i)
