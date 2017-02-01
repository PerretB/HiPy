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
from HiPy.Hierarchies.BinaryPartitionTree import constructBPT

import HiPy.Processing.Attributes as Attributes
from collections import deque
import math
import HiPy.Util.VMath as VMath
from HiPy.Structures import WeightedAdjacency

__author__ = 'Benjamin Perret'


def constructMumfordShahEnergyBPT(image, baseAdjacency):
    adjacency = WeightedAdjacency.createAdjacency(baseAdjacency)

    vectorial = type(image[0]) is list
    # attributes on regions
    regionStats = []  # (area, perimeterLength, mean, mean�)
    optimalEnergies = []

    # attributes on edges
    edgeLength = [1] * len(adjacency)

    def computeLinearPieceForStatVectorial(regionStat):
        return LinearPiece(0,  regionStat[0]*sum(VMath.diffV(regionStat[3], VMath.multV(regionStat[2],regionStat[2]))),
                           regionStat[1])

    def computeLinearPieceForStatScalar(regionStat):
        return LinearPiece(0, regionStat[0] * (regionStat[3] - regionStat[2] * regionStat[2]),
                           regionStat[1])

    def computeFusionStatsScalar(region1, region2, fusionEdge):
        statR1 = regionStats[region1]
        statR2 = regionStats[region2]
        area = statR1[0] + statR2[0]
        perimeterLength = statR1[1] + statR2[1] - 2 * edgeLength[fusionEdge]
        mean = (statR1[2] * statR1[0] + statR2[2] * statR2[0]) / area
        mean2 = (statR1[3] * statR1[0] + statR2[3] * statR2[0]) / area
        return [area, perimeterLength, mean, mean2]

    def computeFusionStatsVectorial(region1, region2, fusionEdge):
        statR1 = regionStats[region1]
        statR2 = regionStats[region2]
        area = statR1[0] + statR2[0]
        perimeterLength = statR1[1] + statR2[1] - 2 * edgeLength[fusionEdge]
        mean = VMath.divS(VMath.addV(VMath.multS(statR1[2], statR1[0]), VMath.multS(statR2[2], statR2[0])), area)
        mean2 = VMath.divS(VMath.addV(VMath.multS(statR1[3], statR1[0]), VMath.multS(statR2[3], statR2[0])), area)
        return [area, perimeterLength, mean, mean2]

    if vectorial:
        computeLinearPieceForStat = computeLinearPieceForStatVectorial
        computeFusionStats = computeFusionStatsVectorial
    else:
        computeLinearPieceForStat = computeLinearPieceForStatScalar
        computeFusionStats = computeFusionStatsScalar

    def computeFusionEnergy(region1, region2, fusionEdge):
        statFusion = computeFusionStats(region1, region2, fusionEdge)
        energyFun = optimalEnergies[region1].sum(optimalEnergies[region2])
        apparitionScale = energyFun.infimum(computeLinearPieceForStat(statFusion))
        return apparitionScale, energyFun

    for i in image.iterateOnPixels():
        v = image[i]
        if vectorial:
            regionStats.append([1, baseAdjacency.countOutEdges(i, includeExternal=True), v, VMath.multV(v, v)])
        else:
            regionStats.append([1, baseAdjacency.countOutEdges(i, includeExternal=True), v, v * v])
        optimalEnergies.append(PiecewiseLinearEnergyFunction(
            computeLinearPieceForStat(regionStats[i])))

    edgeLength = [1] * len(adjacency)
    for i in adjacency.iterateOnPixels():
        adjacency[i] = computeFusionEnergy(adjacency.source[i], adjacency.target[i], i)[0]

    def computeFusionWeights(image, adjacency, fusionEdge, neighbourList, inEdges):

        # attribute merged regions
        newRegion = len(regionStats)
        region1 = adjacency.source[fusionEdge]
        region2 = adjacency.target[fusionEdge]
        statNewRegion = computeFusionStats(region1, region2, fusionEdge)
        regionStats.append(statNewRegion)
        optimalEnergies.append(computeFusionEnergy(region1, region2, fusionEdge)[1])

        weights = []
        # attribute edges and reweighting
        for n in neighbourList:
            edgeLength[inEdges[n][0]] = sum([edgeLength[j] for j in inEdges[n]])

            weight = computeFusionEnergy(newRegion, n, inEdges[n][0])[0]

            weights.append(weight)
        return weights

    return constructBPT(image, adjacency, computeFusionWeights)


def trace(mat):
    l = len(mat)
    c = 0
    for i in range(l):
        c += mat[i][i]
    return c


@Attributes.autoCreateAttribute("mumfordShahEnergy", False)
def addAttributeMumfordShahEnergy(tree, attribute, lambdaValue, levelImage):
    """
    For partition hierarchy
    TODO: check adaptation for component tree
    :param tree:
    :param attribute:
    :param lambdaValue:
    :param levelImage:
    :return:
    """
    area = Attributes.addAttributeArea(tree)
    perimeter = Attributes.addAttributePerimeterPartitionHierarchy(tree, levelImage.adjacency)
    levelStatistics = Attributes.addAttributeLevelStatistics(tree, levelImage)
    if type(levelStatistics[0][1]) is list:
        for i in tree.iteratorFromPixelsToRoot():
            attribute[i] = trace(levelStatistics[i][1]) * area[i] + lambdaValue * perimeter[i]
    else:
        for i in tree.iteratorFromPixelsToRoot():
            attribute[i] = levelStatistics[i][1] * area[i] + lambdaValue * perimeter[i]


@Attributes.autoCreateAttribute("energyCut", 0)
def addAttributeEnergyCut(tree: "Tree", attribute, energyAttributeName="energy"):
    """
    Computes the cut that maximises the given energy attribute.

    See Guigues Thesis
    :param tree:
    :param attribute:
    :param energyAttributeName:
    :return:
    """
    optimalNodes = attribute.getCopy(False)
    optimalEnergy = attribute.getCopy(False)
    children = Attributes.addAttributeChildren(tree)
    energyAttribute = tree.getAttribute(energyAttributeName)
    for i in tree.iterateOnPixels():
        optimalNodes[i] = True
        optimalEnergy[i] = energyAttribute[i]

    for i in tree.iteratorFromPixelsToRoot(includePixels=False):
        energyChildren = 0
        for child in children[i]:
            energyChildren += optimalEnergy[child]
        if energyAttribute[i] <= energyChildren:
            optimalNodes[i] = True
            optimalEnergy[i] = energyAttribute[i]
        else:
            optimalNodes[i] = False
            optimalEnergy[i] = energyChildren

    stack = [len(tree) - 1]
    while stack:
        element = stack.pop()
        if not optimalNodes[element]:
            attribute[element] = 1
            stack.extend(children[element])


def transformTreeToOptimalMumfordShahEnergyCutHierarchy(tree, levelImage):
    """
    Transform the given partition tree into an optimal energy cut hierarchy using the piecewise constant
    Mumford-Shah energy.

    levelImage must be quipped with an appropriate adjacency in order to compute perimeter length of regions.

    See transformTreeToOptimalEnergyCutHierarchy

    :param tree: Tree to transform
    :param levelImage: Image used to compute region statistics
    :return: See transformTreeToOptimalEnergyCutHierarchy
    """
    area = Attributes.addAttributeArea(tree)
    stats = Attributes.addAttributeLevelStatistics(tree, levelImage)
    perimeter = Attributes.addAttributePerimeterPartitionHierarchy(tree, levelImage.adjacency, False)
    variance, created = tree.addAttribute("variance", 0)
    for i in tree.iteratorFromPixelsToRoot():
        variance[i] = trace(stats[i][1]) * area[i]
    return transformTreeToOptimalEnergyCutHierarchy(tree, variance, perimeter)


def transformTreeToOptimalEnergyCutHierarchy(tree, dataFidelityAttribute, regularizationAttribute):
    """
    Transform the given hierarchy into its optimal energy cut hierarchy for the given energy terms.

    In the optimal energy cut hierarchy, any horizontal cut corresponds to an optimal energy cut in the original
    hierarchy.

    See Guigues thesis

    :param tree: a partition hierarchy
    :param dataFidelityAttribute:
    :param regularizationAttribute:
    :return:
    """
    children = Attributes.addAttributeChildren(tree)
    optimalEnergies = [None] * len(tree)
    apparitionScales, created = tree.addAttribute("apparitionScales", 0)
    for i in tree.iteratorOnPixels():
        optimalEnergies[i] = PiecewiseLinearEnergyFunction(
            LinearPiece(0, dataFidelityAttribute[i], regularizationAttribute[i]))
        apparitionScales[i] = -dataFidelityAttribute[i] / regularizationAttribute[i]

    for i in tree.iteratorFromPixelsToRoot(includePixels=False):
        energyI = LinearPiece(0, dataFidelityAttribute[i], regularizationAttribute[i])
        energyChildren = optimalEnergies[children[i][0]]
        for childIndex in range(1, len(children[i])):
            child = children[i][childIndex]
            energyChildren = energyChildren.sum(optimalEnergies[child])
        apparitionScales[i] = energyChildren.infimum(energyI)
        optimalEnergies[i] = energyChildren

    __filterNonPersistentNodes(tree, apparitionScales)

    return tree.simplifyTreeByAttribute("apparitionScales", "apparitionScales")


def __filterNonPersistentNodes(tree, apparitionScales):
    for i in tree.iteratorFromRootToPixels(includeRoot=False):
        apparitionScales[i] = max(0, min(apparitionScales[i], apparitionScales[tree[i]]))


class LinearPiece:
    """
    One piece of a PiecewiseLinearEnergyFunction
    """

    def __init__(self, originX, originY, slope):
        self.originX = originX
        self.originY = originY
        self.slope = slope

    def __call__(self, x):
        return self.originY + self.slope * (x - self.originX)

    def __str__(self):
        return "(" + str(self.originX) + ", " + str(self.originY) + ", " + str(self.slope) + ")"

    def __repr__(self):
        return self.__str__()


class PiecewiseLinearEnergyFunction(deque):
    """
    Piecewise linear energy function as modelled in Laurent Guigues thesis.

    An energy function is a concave non decreasing piecewise linear positive function.
    """

    def __init__(self, *args):
        if args:
            self.extend(args)

    def sum(self, plFunction, maxPieces=10):
        """
        Compute the sum between two PiecewiseLinearEnergyFunction.
        The computation is by default limited to the maxPieces largest pieces (right most)
        :param plFunction:
        :param maxPieces:
        :return:
        """
        result = PiecewiseLinearEnergyFunction()
        count = 0
        i1 = len(self) - 1
        i2 = len(plFunction) - 1
        while i1 >= 0 and i2 >= 0 and count < maxPieces:
            piece1 = self[i1]
            piece2 = plFunction[i2]
            newSlope = piece1.slope + piece2.slope
            if piece1.originX > piece2.originX:
                newOriginX = piece1.originX
                newOriginY = piece1.originY + piece2(piece1.originX)
                if piece1.originX == piece2.originX:
                    i2 -= 1
                i1 -= 1
            else:
                newOriginX = piece2.originX
                newOriginY = piece2.originY + piece1(piece2.originX)
                i2 -= 1
            result.appendleft(LinearPiece(newOriginX, newOriginY, newSlope))
            count += 1

        firstPiece = result[0]
        if firstPiece.originX > 0:
            firstPiece.originY -= firstPiece.slope * firstPiece.originX
            firstPiece.originX = 0
        return result

    def infimum(self, linearPiece):
        """
        Infimum between the current piecewise linear function and the given linear piece

        Modification is done in place
        :param linearPiece:
        :return:
        """
        i = len(self) - 1

        lastPiece = self[i]
        if linearPiece.slope == lastPiece.slope:
            y = linearPiece(lastPiece.originX)
            if y > lastPiece.originY:
                return float("inf")
            elif y == lastPiece.originY:
                return lastPiece.originX
            else:
                i -= 1

        flag = True
        while i >= 0 and flag:
            piece = self[i]
            xi = (linearPiece.originX * linearPiece.slope - piece.originX * piece.slope - (
                linearPiece.originY - piece.originY)) / (linearPiece.slope - piece.slope)
            if xi > piece.originX:
                flag = False
            else:
                self.pop()
            i -= 1

        newPiece = LinearPiece(xi, linearPiece(xi), linearPiece.slope)
        self.append(newPiece)

        return xi
