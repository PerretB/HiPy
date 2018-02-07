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
import heapq
from heapq import heappush
from _heapq import heappop


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
Created on 6 juil. 2015

@author: perretb
'''

import HiPy.Structures
import HiPy.Util.Accumulator
from HiPy.Util.Accumulator import BasicAccumulator


def spatialFilter(image: HiPy.Structures.Image,
                  adjacency: HiPy.Structures.AbstractAdjacency,
                  accumulator: HiPy.Util.Accumulator.AbstractAccumulator) \
        -> HiPy.Structures.Image:
    res = image.getCopy(False)
    for i in image.iterateOnPixels():
        accumulator.reset()
        for e in adjacency.getOutEdges(i):
            accumulator.accumulate(image[e[1]], e[2])
        res[i] = accumulator.result()
    return res


def dilationGray(image, structuringElement):
    adj = constructAdjacencyFromStructuringElement(image, structuringElement)
    return spatialFilter(image, adj, BasicAccumulator.getMaxAccumulator())


def erosionGray(image, structuringElement):
    adj = constructAdjacencyFromStructuringElement(image, structuringElement)
    return spatialFilter(image, adj, BasicAccumulator.getMinAccumulator())


def medianGray(image, structuringElement):
    adj = constructAdjacencyFromStructuringElement(image, structuringElement)
    return spatialFilter(image, adj, BasicAccumulator.getMedianAccumulator())


def meanGray(image, structuringElement):
    adj = constructAdjacencyFromStructuringElement(image, structuringElement)
    return spatialFilter(image, adj, BasicAccumulator.getMeanAccumulator())


def constructAdjacencyFromStructuringElement(image, structuringElement):
    return HiPy.Structures.AdjacencyNdRegular(image.embedding, structuringElement)


def getCrossStructuringElement():
    return getDiamondStructuringElement()


def getSquareStructuringElement(maxHalfWidth=1):
    """
    Return a square shape structuring element
    :param maxHalfWidth: maxHalfWidth: the width (and height) of the square will be 2*maxHalfWidth + 1
    :return:
    """
    points = []
    maxWidth = 2 * maxHalfWidth + 1
    for y in range(-maxHalfWidth, maxHalfWidth+1, 1):
        for x in range(-maxHalfWidth, maxHalfWidth + 1, 1):
            points.append((x, y))
    return points


def getDiamondStructuringElement(maxHalfWidth=1):
    """
    Return a diamond shape structuring element
    :param maxHalfWidth: the width (and height) of the diamond will be 2*maxHalfWidth + 1
    :return:
    """
    points = []
    maxWidth = 2 * maxHalfWidth + 1
    for i in range(maxWidth):
        height = i - maxHalfWidth
        halfWidth = maxHalfWidth - abs(i - maxHalfWidth)
        for x in range(-halfWidth, halfWidth+1, 1):
            points.append((x, height))
    return points


def simpleXGradient(image):
    adj = HiPy.Structures.AdjacencyNdRegular(image.embedding, neighbourList=[(-1, 0), (1, 0)], weights=[-1, 1])
    return spatialFilter(image, adj, BasicAccumulator.getWeightedSumAccumulator())


def simpleYGradient(image):
    adj = HiPy.Structures.AdjacencyNdRegular(image.embedding, neighbourList=[(0, -1), (0, 1)], weights=[-1, 1])
    return spatialFilter(image, adj, BasicAccumulator.getWeightedSumAccumulator())


if __name__ == "__main__":
    print(getSquareStructuringElement(1))
    print(getSquareStructuringElement(2))
    print(getSquareStructuringElement(3))