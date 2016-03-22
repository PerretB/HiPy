# Copyright ESIEE (2015)
#
# benjamin.perret@esiee.fr
#
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
Utility methods for region adjacency graphs

A RAG is an image with an adjacency relation and a mapping from the original image pixels to the regions of the RAG
"""
import copy

import HiPy

from HiPy.Structures import WeightedAdjacency, Image


def constructRAGFromAdjacency(adjacency, edgeFrontierFunction=None):
    """
    Construct a Region Adjacency Graph from an adjacency and a function indicating which edges are frontiers

    :param adjacency: instance of AbstractAdjacency
    :param edgeFrontierFunction:  a function that associates True to an edge of adjacency if it is a frontier edge
            If this parameter is None then the default function is lambda edgeIndex: adjacency[edgeIndex] > 0
    :return:
        ragAdj
        pixelMapping,
        edgeMapping
    """
    HiPy.Structures.HiPyLogger.debug("call to: constructRAGFromAdjacency")

    nbPix = adjacency.nbPoints

    pixelMapping = Image(nbPix, -1)
    edgeMapping = Image(len(adjacency), -1)
    if not edgeFrontierFunction:
        edgeFrontierFunction = lambda edgeIndex: adjacency[edgeIndex] > 0

    def exploreComponent(node, regionIndex):
        stack = [node]
        pixelMapping[node] = regionIndex
        while stack:
            top = stack.pop()
            for edgeIndex in adjacency.getEdgeIndices(top):
                if not edgeFrontierFunction(edgeIndex):
                    edge = adjacency.getEdgeFromIndex(edgeIndex)
                    neighbour = edge[0] + edge[1] - top
                    if pixelMapping[neighbour] == -1:
                        pixelMapping[neighbour] = regionIndex
                        stack.append(neighbour)

    regionNum = 0
    for i in range(nbPix):
        if pixelMapping[i] == -1:
            exploreComponent(i, regionNum)
            regionNum += 1

    edgeNumber = 0
    ragAdj = WeightedAdjacency(regionNum)
    visited = [False] * nbPix
    flag = [-1] * regionNum
    ragEdgeDict = {}

    def buildEdges(node):
        nonlocal edgeNumber
        regionIndex = pixelMapping[node]

        visited[node] = True
        stack = [node]
        while stack:
            top = stack.pop()
            for edgeIndex in adjacency.getEdgeIndices(top):
                edge = adjacency.getEdgeFromIndex(edgeIndex)
                neighbour = edge[0] + edge[1] - top

                if pixelMapping[neighbour] == regionIndex:
                    if not visited[neighbour]:
                        pixelMapping[neighbour] = regionIndex
                        visited[neighbour] = True
                        stack.append(neighbour)
                else:
                    adjRegion = pixelMapping[neighbour]
                    if adjRegion > regionIndex:
                        edgeRef = edgeNumber

                        if flag[adjRegion] != regionIndex:
                            ragAdj.createEdge(regionIndex, adjRegion)
                            ragEdgeDict[(regionIndex, adjRegion)] = edgeNumber
                            edgeNumber += 1
                            flag[adjRegion] = regionIndex
                        else:
                            edgeRef = ragEdgeDict[(regionIndex, adjRegion)]

                        edgeMapping[edgeIndex] = edgeRef

    for i in range(nbPix):
        if not visited[i]:
            buildEdges(i)
    ragAdj.attributes["pixelRAGMap"] = pixelMapping
    ragAdj.attributes["edgeRAGMap"] = edgeMapping
    return ragAdj, pixelMapping, edgeMapping


def inverseRAGMapEdges(rag, originalAdjacency, inplace=False, unMappedEdgeValue=0):
    HiPy.Structures.HiPyLogger.debug("call to: inverseRAGMapEdges")
    if inplace:
        result = originalAdjacency
    else:
        result = originalAdjacency.getCopy(False)

    edgeMap = rag.attributes["edgeRAGMap"]

    for i in result.iterateOnPixels():
        if edgeMap[i] >= 0:
            result[i] = rag[edgeMap[i]]
        else:
            result[i] = unMappedEdgeValue

    return result


def inverseRAGMapPixels(rag, ragValues):
    pixelMap = rag.attributes["pixelRAGMap"]

    result = pixelMap.getCopy(False)

    for i in result.iterateOnPixels():
        result[i] = ragValues[pixelMap[i]]

    return result


def accumulateOnRAGPixels(rag, values, accumulator):
    HiPy.Structures.HiPyLogger.debug("call to: accumulateOnRAGPixels")
    nbPix = rag.nbPoints
    image = Image(nbPix, 0, adjacency=rag)
    accumulators = [None] * nbPix
    for i in image.iterateOnPixels():
        accumulators[i] = accumulator.copy()
        accumulators[i].reset()
    pixelMap = rag.attributes["pixelRAGMap"]
    for i in pixelMap.iterateOnPixels():
        accumulators[pixelMap[i]].accumulate(values[i])

    for i in image.iterateOnPixels():
        image[i] = accumulators[i].result()

    return image


def accumulateOnRAGEdges(rag, values, accumulator, inplace=False):
    HiPy.Structures.HiPyLogger.debug("call to: accumulateOnRAGEdges")
    nbEdges = len(rag)
    image = Image(nbEdges, 0)
    if inplace:
        result = rag
    else:
        result = rag.getCopy()
    accumulators = [None] * nbEdges
    for i in result.iterateOnPixels():
        accumulators[i] = accumulator.copy()
        accumulators[i].reset()
    edgeMap = rag.attributes["edgeRAGMap"]
    for i in edgeMap.iterateOnPixels():
        accumulators[edgeMap[i]].accumulate(values[i])

    for i in image.iterateOnPixels():
        result[i] = accumulators[i].result()

    return result

