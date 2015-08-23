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

import HiPy.Processing.Attributes as Attributes
from collections import deque
import math


__author__ = 'Benjamin Perret'


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
            attribute[i] = sum(levelStatistics[i][1]) * area[i] + lambdaValue * perimeter[i]
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
    optimalNodes = attribute.copy(False)
    optimalEnergy = attribute.copy(False)
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
        variance[i] = sum(stats[i][1]) * area[i]
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
    apparitionScales, created = tree.addAttribute("apparitionScales",0)
    for i in tree.iteratorOnPixels():
        optimalEnergies[i] = PiecewiseLinearEnergyFunction(
            LinearPiece(0, dataFidelityAttribute[i], regularizationAttribute[i]))
        apparitionScales[i] = -dataFidelityAttribute[i]/regularizationAttribute[i]

    for i in tree.iteratorFromPixelsToRoot(includePixels=False):
        energyI = LinearPiece(0, dataFidelityAttribute[i], regularizationAttribute[i])
        energyChildren = optimalEnergies[children[i][0]]
        for childIndex in range(1,len(children[i])):
            child = children[i][childIndex]
            energyChildren = energyChildren.sum(optimalEnergies[child])
        apparitionScales[i] = energyChildren.infimum(energyI)
        optimalEnergies[i] = energyChildren

    __filterNonPersistentNodes(tree, apparitionScales)

    return tree.simplifyTreeByAttribute("apparitionScales","apparitionScales")


def __filterNonPersistentNodes(tree, apparitionScales):
    for i in tree.iteratorFromRootToPixels(includeRoot=False):
        apparitionScales[i] = max(0,min(apparitionScales[i], apparitionScales[tree[i]]))


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
