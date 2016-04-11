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

from HiPy.Hierarchies.BinaryPartitionTree import constructBPT, singleLinkageFusionWeights, completeLinkageFusionWeights
from HiPy.Hierarchies.Quadtree import constructBinaryQuadTree
from HiPy.Hierarchies.WatershedHierarchy import drawSaliencyForVisualisation, constructAltitudeBPT
from HiPy.IO import readImage, saveImage
from HiPy.Processing.Attributes import addAttributeHOG, addAttributeSimpleMoments2d, addAttributeArea
from HiPy.Structures import AdjacencyNdRegular, WeightedAdjacency
from HiPy.Util.Color import convertRGBtoLAB
from HiPy.Util.Histogram import normalizeToByte
from HiPy.Util.VMath import euclideanDistance
import math
__author__ = 'Benjamin Perret'


# Breisen... WHAT !
def drawLine(image, val, x1, y1, x2, y2):

    x1 = int(x1)
    y1 = int(y1)
    x2 = int(x2)
    y2 = int(y2)

    if abs(x1 - x2) <= abs(y1 - y2):
        if y1 > y2:
            x1, x2 = x2, x1
            y1, y2 = y2, y1

        slope = (x2 - x1) / (y2 - y1)
        for i in range(0, y2-y1+1):
            y = int(y1 + i)
            x = int(i * slope + x1)
            image.setPixelWCS(max(val, image.getPixelWCS(x, y)), x, y)

    else:
        if x1 > x2:
            x1, x2 = x2, x1
            y1, y2 = y2, y1

        slope = (y2 - y1) / (x2 - x1)
        for i in range(0, x2-x1+1):
            x = int(x1 + i)
            y = int(i * slope + y1)
            image.setPixelWCS(max(val, image.getPixelWCS(x, y)), x, y)

def demoHOG():
    print("reading image...")
    image = readImage("../samples/elefant.bmp", grayScale=True)

    bpt = constructBinaryQuadTree(image)

    hog = addAttributeHOG(bpt, grayImage=image)

    moments = addAttributeSimpleMoments2d(bpt)

    result = image.getCopy(False)

    area = addAttributeArea(bpt)
    scale1 = 8
    scale = scale1*scale1+1
    binwidth = math.pi/len(hog[0])
    k = 4

    for i in bpt.iteratorFromPixelsToRoot(includeRoot=False):
        if area[i] >= scale < area[bpt[i]]:
            h = hog[i]
            xmean = moments[i][1]/moments[i][0]
            ymean = moments[i][2]/moments[i][0]

            for j in range(len(h)):
                angle = (j+0.5)*binwidth + math.pi/2 #to have something tangent to border
                x1 = math.cos(angle)*k
                y1 = math.sin(angle)*k
                x2 = -x1
                y2 = -y1
                x1 += xmean
                x2 += xmean
                y1 += ymean
                y2 += ymean

                drawLine(result, h[j], x1, y1, x2, y2)

    result = normalizeToByte(result)

    saveImage(result, "Results/demoHog.png")




def main():
    demoHOG()


if __name__ == '__main__':
    main()
