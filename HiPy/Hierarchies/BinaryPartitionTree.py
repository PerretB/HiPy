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

from HiPy.Structures import Tree, TreeType

__author__ = 'Benjamin Perret'

import functools
import heapdict


def constructBPT(image, baseAdjacency, computeFusionWeights):
    adjacency = baseAdjacency.getCopy(optimizedRemove=True, copyData=True)

    parent = [-1] * len(image)
    level = [0] * len(image)
    active = [False] * len(adjacency)
    inEdges = []
    for _ in range(len(image)):
        inEdges.append([])

    heap = heapdict.heapdict()
    __initHeap(image, adjacency, heap, active)

    while len(heap) > 0:
        fusionEdge, weight = heap.popitem()

        if active[fusionEdge]:
            newParent = adjacency.addVertex()
            region1, region2, weight = adjacency.getEdgeFromIndex(fusionEdge)
            parent[region1] = newParent
            parent[region2] = newParent
            parent.append(-1)
            level.append(weight)
            inEdges.append([])

            adjacency.removeEdge(fusionEdge)

            deletedEdges = adjacency.edgeList[region1] | adjacency.edgeList[region2]

            neighbourEdges = __findNeighbours(adjacency, deletedEdges, region1, region2, active, inEdges)

            if len(neighbourEdges) > 0:
                __computeNewEdges(image, adjacency, neighbourEdges, fusionEdge, computeFusionWeights, newParent, heap,
                                  active, inEdges)

    return Tree(TreeType.PartitionHierarchy, parent, level, image)


def __initHeap(image, adjacency, heap, active):
    """
    For each vertex of the adjacency, inserts the out edge of minimal weight in the heap and set its state as active.
    :param image:
    :param adjacency:
    :param heap:
    :param active:
    :return:
    """
    for i in image.iterateOnPixels():
        eightMin = 99999999999
        edgeIndex = -1

        for edge in adjacency.edgeList[i]:
            source, target, weight = adjacency.getEdgeFromIndex(edge)
            if target > source and weight < eightMin:
                eightMin = weight
                edgeMin = edge

        heap[edgeMin] = eightMin
        active[edgeMin] = True


def __findNeighbours(adjacency, deletedEdges, region1, region2, active, inEdges):
    """
    Find all vertex v of adjacency such that there exists an edge (region1,v) or region(2,v).
    The value associated to v in the list inEdges is the list of
    edges coming from {region1, region2} to v.

    All the edges going out of {region1, region2} are set inactive
    :param adjacency:
    :param deletedEdges:
    :param region1:
    :param region2:
    :param active:
    :return:
    """

    neighbourList = []
    adjacency_source = adjacency.source
    adjacency_target = adjacency.target
    for edge in deletedEdges:

        source = adjacency_source[edge]
        target = adjacency_target[edge]

        active[edge] = False

        if source == region1 or source == region2:
            newNeighbour = target
        else:
            newNeighbour = source

        if len(inEdges[newNeighbour]) == 0:
            neighbourList.append(newNeighbour)
        inEdges[newNeighbour].append(edge)

    return neighbourList


def __computeNewEdges(image, adjacency, neighbourList, fusionEdge, computeFusionWeights, newParent, heap, active,
                      inEdges):
    weights = computeFusionWeights(image, adjacency, fusionEdge, neighbourList, inEdges)

    weightMin = 99999999
    edgeMin = -1
    for i in range(len(neighbourList)):
        neighbour = neighbourList[i]
        weight = weights[i]
        fusedEdges = inEdges[neighbour]

        for i in range(1, len(fusedEdges)):
            adjacency.removeEdge(fusedEdges[i])

        edgeIndex = fusedEdges[0]
        adjacency.setEdge(edgeIndex, newParent, neighbour, weight)

        fusedEdges.clear()

        if weight < weightMin:
            weightMin = weight
            edgeMin = edgeIndex

    heap[edgeMin] = weightMin
    active[edgeMin] = True


def computeFusionWeightsLinkage(image, adjacency, fusionEdge, neighbourList, inEdges, fusionFunction):
    weights = []
    for neighbour in neighbourList:
        newWeight = fusionFunction([adjacency[e] for e in inEdges[neighbour]])
        weights.append(newWeight)
    return weights


singleLinkageFusionWeights = functools.partial(computeFusionWeightsLinkage, fusionFunction=min)

completeLinkageFusionWeights = functools.partial(computeFusionWeightsLinkage, fusionFunction=max)
