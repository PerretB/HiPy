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
Created on 15 june 2015

@author: Benjamin Perret
"""

import HiPy.Util.UnionFind as UnionFind
import HiPy.Structures
import HiPy.Processing.Attributes
from HiPy.Util.Histogram import imageMap, rescaleGray, normalizeToByte


def constructAltitudeBPT(adjacency: "AbstractWeightedAdjacency") -> "Tree":
    """
    Construct the binary partition tree by altitude of the given edge image.

    The corresponding MST is stored in the attribute leavesAdjacency of the returned tree

    :param adjacency: a weighted adjacency (an edge valued graph)
    :returns the binary partition tree by altitude of the given adjacency
    """
    HiPy.Structures.HiPyLogger.debug("call to: constructAltitudeBPT")

    nbPoints = adjacency.nbPoints
    edgesI = sorted(range(len(adjacency)), key=lambda x: adjacency[x])
    edges = [[adjacency.source[i], adjacency.target[i], adjacency[i]] for i in edgesI]

    mst, parent = computeMSTBPT(nbPoints, edges)

    levels = [0] * nbPoints
    adjMST = HiPy.Structures.WeightedAdjacency(nbPoints)
    for e in mst:
        levels.append(e[2])
        adjMST.createEdge(*e)

    tree = HiPy.Structures.Tree(HiPy.Structures.TreeType.PartitionHierarchy, parent, levels)
    tree.leavesAdjacency = adjMST
    return tree


def transformAltitudeBPTtoComponentTree(bpt: "Tree") -> "Tree":
    """
    Copy bpt and delete the nodes n such that bpt.level[n]=bpt.level[bpt[n]]
    (and update the parent relation accordingly...

    :param bpt: a binary partition tree
    :returns the partition tree corresponding to the binary partition tree
    """
    nbLeaves = bpt.nbPixels
    nbNodes = len(bpt)
    children = HiPy.Processing.Attributes.addAttributeChildren(bpt)
    level = bpt.level

    count = 0
    deleted = [False] * nbNodes
    deletedMap = [0] * nbNodes

    # from root to leaves, compute the new parent relation,
    # don't care of the holes in the parent tab
    for i in range(nbNodes - 2, nbLeaves - 1, -1):
        par = bpt[i]
        if level[i] == level[par]:
            for c in children[i]:
                bpt[c] = par
            deleted[i] = True
            count += 1
        # inverse of what we want: number of deleted nodes after node i
        deletedMap[i] = count

    # correct the mapping
    for i in bpt.iteratorFromPixelsToRoot(False):
        deletedMap[i] = count - deletedMap[i]

    # new relations with correct size
    newParent = [-1] * (nbNodes - count)
    newLevel = [0] * (nbNodes - count)

    count = 0
    for i in range(0, nbNodes - 1):
        if not deleted[i]:
            par = bpt[i]
            newPar = par - deletedMap[par]
            newParent[count] = newPar
            newLevel[count] = level[i]
            count += 1
    newParent[count] = -1
    newLevel[count] = level[-1]

    newTree = HiPy.Structures.Tree(HiPy.Structures.TreeType.PartitionHierarchy, newParent, newLevel)
    newTree.leavesAdjacency = bpt.leavesAdjacency

    return newTree


def transformBPTtoAttributeHierarchy(bpt: "Tree", attributeName: str) -> "Tree":
    """
    Transforms the given binary partition tree (bpt) according to the given attribute.

    The edges of the minimum spanning tree associated to the binary partition tree are reweighted by the attribute
    values, and the binary partition tree of the reweighted minimum spanning tree is returned.

    :param bpt: a binary partition tree and its associated minimum spanning tree (attribute leavesAdjacency)
    :param attributeName: the attribute used to reweight the minimum spanning tree
    :return: the new binary partition tree
    """
    adj = reweightMSTByAttribute(bpt, attributeName)
    return constructAltitudeBPT(adj)


def filterBPTbyCriterion(bpt: "Tree", filterCriterion: "function int->bool") -> "Tree":
    """
    Transforms the given binary partition tree (bpt) according to the given criterion.

    The edges of the minimum spanning tree associated to the binary partition tree are set to 0 when the criterion is
    not satisfied, and the binary partition tree of the reweighted minimum spanning tree is returned.
    :param bpt: a binary partition tree and its associated minimum spanning tree (attribute leavesAdjacency)
    :param filterCriterion: The criterion is a function that associates True or False to any inner node of \
    the binary partition tree. Given an inner node i of the bpt:
        * if filterCriterion(i)==True: the weight of the corresponding edge in the MST is set to 0
        * otherwise the weight of the corresponding edge in the MST is left unchanged
    :return: the new binary partition tree
    """
    adj = filterMSTByCriterion(bpt, filterCriterion)
    return constructAltitudeBPT(adj)


def transformAltitudeBPTtoWatershedHierarchy(bpt: "Tree") -> "Tree":
    """
    Transforms the given binary partition tree (bpt) according to the given criterion.

    The edges of the minimum spanning tree associated to the binary partition tree are set to 0 when they do not belong
    to a watershed cut, and the binary partition tree of the reweighted minimum spanning tree is returned.
    :param bpt: a binary partition tree and its associated minimum spanning tree (attribute leavesAdjacency)
    :return: the new binary partition tree
    """
    adj = extractWatershedEdges(bpt)
    return constructAltitudeBPT(adj)


def computeMSTBPT(nbPoints: int, sortedEdgeList: list) -> (list, list):
    """
    Compute the minimum spanning tree and the binary partition tree associated to a weighted adjacency given by the list
    of its edges sorted by increasing weight.

    :param nbPoints: number of vertices in the adjacency
    :param sortedEdgeList: list of edges, each edge is in the form (source, destination, weight)
    :return: a list of edge composing the MST, a parent relation representing the BPT
    """
    HiPy.Structures.HiPyLogger.debug("call to: computeMSTBPT")
    nbEdgeMST = nbPoints - 1

    mst = []
    parent = [-1] * nbPoints
    ufParent = [i for i in range(nbPoints)]
    ufRank = [0] * nbPoints
    root = [i for i in range(nbPoints)]
    nbEdgeFound = 0
    i = 0
    while nbEdgeFound < nbEdgeMST:
        e = sortedEdgeList[i]
        comp1 = UnionFind.findTarjan(e[0], ufParent)
        comp2 = UnionFind.findTarjan(e[1], ufParent)
        if comp1 != comp2:
            newParent = len(parent)
            parent.append(-1)
            parent[root[comp1]] = newParent
            parent[root[comp2]] = newParent
            newRoot, _ = UnionFind.unionTarjan(comp1, comp2, ufParent, ufRank)
            root[newRoot] = newParent
            mst.append(e)
            nbEdgeFound += 1
        i += 1

    return mst, parent


def extractWatershedEdges(bpt: "Tree") -> "AbstractWeightedAdjacency":
    """
    Return a copy of the minimum spanning tree associated to the binary partition tree (bpt.leavesAdjacency)
    where the weight of the non-watershed edges are set to 0
    :param bpt: a binary partition tree
    :return: the reweighted minimum spanning tree
    """
    adj = bpt.leavesAdjacency.getCopy()
    nbLeaves = bpt.nbPixels
    nbNodes = len(bpt)
    children = HiPy.Processing.Attributes.addAttributeChildren(bpt)
    level = bpt.level

    minima = [0] * nbNodes
    for i in bpt.iteratorFromPixelsToRoot(False):
        childList = children[i]
        min1 = minima[childList[0]]
        min2 = minima[childList[1]]
        numMin = min1 + min2

        if min1 == 0 or min2 == 0:
            adj[i - nbLeaves] = 0

        if numMin != 0:
            minima[i] = numMin
        else:
            par = bpt[i]
            # if p==-1 there is a unique minima in the image
            if par == -1 or level[i] != level[par]:
                minima[i] = 1
                # else minima[i]=0

    return adj


def correctAttributeValueBPT(bpt: "Tree", attributeName: str) -> None:
    """
    Most attributes computed using functions in module HiPy.Processing.Attributes
    will be incorrect for the BPT. However their values can be corrected a posteriori
    using this function !

    TODO: add precisions!

    :param bpt: a binary partition tree
    :param attributeName: name of the attribute
    :return: nothing, modification is done in place
    """
    children = HiPy.Processing.Attributes.addAttributeChildren(bpt)
    attr = bpt.getAttribute(attributeName)
    level = bpt.level
    for i in bpt.iterateOnLeaves():
        attr[i] = 0

    for i in bpt.iteratorFromPixelsToRoot(False):
        par = bpt[i]
        if par != -1 and level[i] == level[par]:
            vc1 = attr[children[i][0]]
            vc2 = attr[children[i][1]]
            attr[i] = max(vc1, vc2)


def reweightMSTByAttribute(bpt: "Tree", attributeName: str, extinctionValue: bool=True):
    """
    Reweight the MST associated to a binary partition tree according to a given attribute.

    If the attribute represents an extinctionValue, then the new weight of an edge is the
    minimum of the attribute values of the two regions separated by this edge.

    Otherwise the new value of the edge is simply the attribute value of the node associated
    to this edge.

    The function returns a new weighted adjacency containing with the same edges as the MST but
    different weights.

    :param bpt: a binary partition tree
    :param attributeName: name of the attribute
    :param extinctionValue: indicates if the attribute represents an extinction value
    :return: a copy of the MST associated to bpt weighted by the attribute
    """
    nbLeaves = bpt.nbPixels
    adj = bpt.leavesAdjacency.getCopy()
    children = HiPy.Processing.Attributes.addAttributeChildren(bpt)
    attr = bpt.getAttribute(attributeName)
    if extinctionValue:
        for i in bpt.iteratorFromPixelsToRoot(False):
            adj[i - nbLeaves] = min(attr[children[i][0]], attr[children[i][1]])
    else:
        for i in bpt.iteratorFromPixelsToRoot(False):
            adj[i - nbLeaves] = attr[i]
    return adj


def filterMSTByCriterion(bpt: "Tree", filterCriterion: "function int->bool") \
        -> "AbstractWeightedAdjacency":
    """
    Filter the MST associated to a binary partition tree according to a given criterion.

    The criterion is a function that associates True or False to any inner node of
    the binary partition tree

    Given an inner node i of the bpt:
      - if filterCriterion(i)==True: the weight of the corresponding edge in the
      MST is set to 0
      - otherwise the weight of the corresponding edge in the MST is left unchanged

    :param bpt: a binary partition tree
    :param filterCriterion: : a filtering criterion
    :return: the reweighted minimum spanning tree
    """
    nbLeaves = bpt.nbPixels
    adj = bpt.leavesAdjacency.getCopy()
    HiPy.Processing.Attributes.addAttributeChildren(bpt)

    for i in bpt.iteratorFromPixelsToRoot(False):
        if filterCriterion(i):
            adj[i - nbLeaves] = 0

    return adj


def computeSaliencyMapFromAttribute(partitionTree: "Tree",
                                    adjacency: "AbstractWeightedAdjacency",
                                    attribute="level"):
    """
    Compute the saliency values for the given attribute of the edges
    of the given adjacency w.r.t the given partition tree.
    :param partitionTree: a partition hierarchy
    :param adjacency: an adjacency in the pixel domain of the partition hierarchy
    :param attribute: the attribute used to valuate the saliency (should be increasing)
    :return: A copy of the given adjacency weighted by the saliency values
    """
    attr = partitionTree.getAttribute(attribute)
    valFun = lambda i: attr[i]
    return computeSaliencyMap(partitionTree, adjacency, valFun)


def computeSaliencyMap(partitionTree: "Tree",
                       adjacency: "AbstractWeightedAdjacency",
                       valuationFunction):
    """
    Compute the saliency values for the given valuation function of the edges
    of the given adjacency w.r.t the given partition tree.
    :param partitionTree: a partition hierarchy
    :param adjacency: an adjacency in the pixel domain of the partition hierarchy
    :param valuationFunction: a function associating a value to a node index
    :return: A copy of the given adjacency weighted by the saliency values
    """
    HiPy.Processing.Attributes.addAttributeDepth(partitionTree)
    lca = partitionTree.lca
    valFun = lambda i, j: valuationFunction(lca(i, j))
    return HiPy.Structures.WeightedAdjacency.createAdjacency(adjacency, valFun)


def drawSaliencyMap(size: (int, int),
                    saliency: "AbstractAdjacency",
                    interpolationFunction: "function *float->float"=max) -> "Image":
    """
    Represent a saliency map as a contour image.
    :param size: is the size [width, height] of the image.
    :param saliency: must represent a 4 adjacency on the 2d grid, results are unpredictable otherwise
    :param interpolationFunction: is used to interpolate values on the 0-faces of the Khalimsky grid. \
        Its argument is a list containing the value of the neighbours 1-faces
    :returns An contour map image of size [2*width-1,2*height-1]
    """
    width = size[0]
    height = size[1]
    grid = HiPy.Structures.Embedding2dGrid(width, height)
    resWidth = width * 2 - 1
    resHeight = height * 2 - 1
    grid2 = HiPy.Structures.Embedding2dGrid(resWidth, resHeight)

    res = HiPy.Structures.Image(resWidth * resHeight,
                0,
                HiPy.Structures.AdjacencyNdRegular.getAdjacency2d4([resWidth, resHeight]),
                grid2)
    for y in range(height):
        for x in range(width):
            pixLin = grid.getLinearCoordinate(x, y)
            for n in saliency.getOutEdges(pixLin):
                if n[1] > pixLin:
                    neighbourCoord = grid.fromLinearCoordinate(n[1])
                    res.setPixelWCS(n[2], x + neighbourCoord[0], y + neighbourCoord[1])

    for y in range(1, resHeight - 1, 2):
        for x in range(1, resWidth - 1, 2):
            pixLin = grid2.getLinearCoordinate(x, y)
            values = []
            for n in res.getNeighbours(pixLin):
                values.append(res[n])
            res[pixLin] = interpolationFunction(values)

    return res


def drawSaliencyForVisualisation(tree: "Tree", image: "Image", attributeName="level", gammaFactor=0.33333):
    """
    Draw the saliency map associated to the given partition tree in a byte image.

    The tree is assumed to represent the given image.

    Only saliency associated to a 4 adjacency can be drawn consistently.

    :param tree: a partition tree
    :param image: the image associated to the tree
    :param attributeName: name of the attribute used to compute saliency (should be increasing)
    :param gammaFactor: a gamma factor applied to the saliency value
    :return a contour map image with byte values
    """
    adj4 = HiPy.Structures.AdjacencyNdRegular.getAdjacency2d4(image.embedding.size)

    saliency = computeSaliencyMapFromAttribute(tree, adj4, attributeName)
    sal = drawSaliencyMap(image.embedding.size, saliency)

    sal = rescaleGray(sal, 0, 1)
    sal = imageMap(sal, lambda x: x ** gammaFactor)
    sal = normalizeToByte(sal)
    return sal
