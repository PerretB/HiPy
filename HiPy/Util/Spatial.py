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
from HiPy.IO import readImage, saveImage

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
from HiPy.Util.Histogram import normalizeToByte

'''
Created on 6 juil. 2015

@author: perretb
'''

from HiPy.Structures import AdjacencyNdRegular, WeightedAdjacency
from HiPy.Util.Accumulator import BasicAccumulator


def spatialFilter(image, adjacency, accumulator):
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
    return AdjacencyNdRegular(image.embedding, structuringElement)


def getCrossStructuringElement():
    return [(0, -1), (-1, 0), (0, 0), (1, 0), (0, 1)]


def getSquareStructuringElement():
    return [(-1, -1), (0, -1), (1, -1), (-1, 0), (0, 0), (1, 0), (-1, 1), (0, 1), (1, 1)]


def simpleXGradient(image):
    adj = AdjacencyNdRegular(image.embedding, neighbourList=[(-1, 0), (1, 0)], weights=[-1, 1])
    return spatialFilter(image, adj, BasicAccumulator.getWeightedSumAccumulator())


def simpleYGradient(image):
    adj = AdjacencyNdRegular(image.embedding, neighbourList=[(0, -1), (0, 1)], weights=[-1, 1])
    return spatialFilter(image, adj, BasicAccumulator.getWeightedSumAccumulator())

    # image = readImage("../../samples/lennaGray256.png")
    # image.adjacency = AdjacencyNdRegular.getAdjacency2d4(image.embedding.size)
    # bigSquareAdj = WeightedAdjacency.createKHopAjacency(image.adjacency, 6, lambda *_:1)
    # res = spatialFilter(image,bigSquareAdj,BasicAccumulator.getMaxAccumulator())
#     # saveImage(res,"bigdil.png")
# image = readImage("../../samples/lennaGray256.png", grayScale=True)
# xg = simpleXGradient(image)
# yg = simpleYGradient(image)
#
# xg = normalizeToByte(xg)
# yg = normalizeToByte(yg)
# saveImage(xg, "f:/xg.png")
# saveImage(yg, "f:/yg.png")