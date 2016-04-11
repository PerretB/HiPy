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
Created on 19 juin 2015

@author: perretb
"""

import random
from math import *  # @UnusedWildImport
import HiPy.Hierarchies.WatershedHierarchy
import HiPy.Structures
import HiPy.Util.VMath as VMath
from functools import wraps
import inspect

from HiPy.Util.Spatial import simpleXGradient, simpleYGradient


def autoCreateAttribute(defaultName, defaultValue=0):
    def decorator(fun):
        # argsSpecFun = inspect.getfullargspec(fun)
        # argsWithDefaultVal = dict(zip(argsSpecFun.args[-len(argsSpecFun.defaults):],argsSpecFun.defaults))
        # if "attribute" not in argsWithDefaultVal:
        #    raise Exception("HiPy Library error: method fun should give a default attribute name value")
        # defaultAttribute = argsWithDefaultVal["attribute"]

        @wraps(fun)
        def wrapper(tree, *args, **kwargs):
            attributeName = kwargs.pop("attributeName", defaultName)
            attr, created = tree.addAttribute(attributeName, defaultValue)
            if created:
                fun(tree, attr, *args, **kwargs)
            return attr

        wrapper.original = fun
        return wrapper

    return decorator


def updateAttributeAfterFiltering(tree, attributeName, defaultValue=0):
    """
    Update attribute values according to a filtering.
    
    The attribute value of each node marked as deleted is set to the attribute value of its closest non deleted ancestor.
    
    If no such ancestor exist, its attribute value is set to "defaultValue".
    """
    deleted = tree.deleted
    attr = tree.getAttribute(attributeName)
    for i in tree.iteratorFromRootToPixels():
        if deleted[i]:
            par = tree[i]
            if par != -1:
                attr[i] = attr[par]
            else:
                attr[i] = defaultValue


@autoCreateAttribute("randomColor", None)
def addAttributeRandomColor(tree: "HiPy.Structures.Tree", attribute):
    for i in tree.iteratorFromPixelsToRoot(False):
        attribute[i] = [random.randint(100, 255), random.randint(100, 255), random.randint(100, 255)]


@autoCreateAttribute("area", 0)
def addAttributeArea(tree: "HiPy.Structures.Tree", attribute):
    for i in tree.iteratorFromPixelsToRoot(True):
        if i < tree.nbPixels:
            attribute[i] = 1
        par = tree[i]
        if par != -1:
            attribute[par] += attribute[i]


@autoCreateAttribute("children", [])
def addAttributeChildren(tree, attribute):
    for i in tree.iteratorFromPixelsToRoot(True):
        par = tree[i]
        if par != -1:
            attribute[par].append(i)


@autoCreateAttribute("childrenLogical", [])
def addAttributeChildrenLogical(tree, attribute):
    if tree.treeType == HiPy.Structures.TreeType.ComponentTree:
        for i in tree.iteratorFromLeavesToRoot():
            par = tree[i]
            if par != -1:
                attribute[par].append(i)

    else:
        addAttributeChildren(tree, attribute)


@autoCreateAttribute("volume", 0)
def addAttributeVolume(tree, attribute):
    area = addAttributeArea(tree)

    for i in tree.iteratorFromPixelsToRoot(True):
        par = tree[i]
        if par != -1:
            lvl = tree.level[i]
            lvlPar = tree.level[par]
            v = abs(lvl - lvlPar) * area[i]
            attribute[i] = attribute[i] + v
            attribute[par] = attribute[par] + v


@autoCreateAttribute("marker", False)
def addAttributeMarker(tree, attribute, markerImage, markerValue):
    """
    Compute the attribute "attributeName" on the graph.
    For each node, the attribute is set to true if the  connected 
    component of the node contains a point of the image "markerImage" 
    having the value "markerValue"
    """
    children = addAttributeChildren(tree)

    for i in range(tree.nbPixels):
        attribute[i] = (markerImage[i] == markerValue)

    for i in tree.iteratorFromPixelsToRoot(False):
        flag = False
        for j in children[i]:
            if attribute[j]:
                flag = True
                break
        attribute[i] = flag


@autoCreateAttribute("moments", None)
def addAttributeSimpleMoments2d(tree, attribute):
    """
    Compute the "raw" moments [M00 M10 M01 M11 M20 M02 M21 M12 M30 M03]  of each component
    """
    nbLeaves = tree.nbPixels
    embedding = tree.leavesEmbedding

    for i in tree.iteratorFromPixelsToRoot():
        if i < nbLeaves:
            c = embedding.fromLinearCoordinate(i)
            x = c[0]
            y = c[1]
            xx = x * x
            yy = y * y
            attribute[i] = [1, x, y, x * y, xx, yy, xx * y, x * yy, x * xx, y * yy]
        else:
            attribute[i] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    for i in tree.iteratorFromPixelsToRoot():
        par = tree[i]
        if par != -1:
            m = attribute[i]
            mp = attribute[par]
            for j in range(len(m)):
                mp[j] += m[j]


def computeCentralMoments2d(rawMoment):
    """                           0   1   2   3   4   5   6   7   8   9
    Compute the central moments [u00 u10 u01 u11 u20 u02 u21 u12 u30 u03]  from the raw moments
    """
    m = rawMoment
    mx = m[1] / m[0]
    my = m[2] / m[0]
    u00 = m[0]
    u01 = u10 = 0
    u11 = m[3] - mx * m[2]
    u20 = m[4] - mx * m[1]
    u02 = m[5] - my * m[2]
    u21 = m[6] - 2 * mx * m[3] - my * m[4] + 2 * mx * mx * m[2]
    u12 = m[7] - 2 * my * m[3] - mx * m[5] + 2 * my * my * m[1]
    u30 = m[8] - 3 * mx * m[4] + 2 * mx * mx * m[1]
    u03 = m[9] - 3 * my * m[5] + 2 * my * my * m[2]
    return [u00, u10, u01, u11, u20, u02, u21, u12, u30, u03]


def computeScaleInvariantMoments2d(centralMoment):
    """                             3+    0   1   2   3   4   5   6   
    Compute the scale invariant moments [n11 n20 n02 n21 n12 n30 n03]  from the central moments
    """
    u00 = centralMoment[0]

    def sc(u, i, j):
        return u / u00 ** (1 + (i + j) / 2)

    n11 = sc(centralMoment[3], 1, 1)
    n20 = sc(centralMoment[4], 2, 0)
    n02 = sc(centralMoment[5], 0, 2)
    n21 = sc(centralMoment[6], 2, 1)
    n12 = sc(centralMoment[7], 1, 2)
    n30 = sc(centralMoment[8], 3, 0)
    n03 = sc(centralMoment[9], 0, 3)
    return [n11, n20, n02, n21, n12, n30, n03]


def computeHuInvariant(scaleInvariantMoments):
    n = scaleInvariantMoments
    #  0   1   2   3   4   5   6
    # [n11 n20 n02 n21 n12 n30 n03]
    I1 = n[1] + n[2]
    I2 = I1 ** 2 + 4 * n[0] ** 2
    I3 = (n[5] - 3 * n[4]) ** 2 + (3 * n[3] - n[6]) ** 2
    I4 = (n[5] + n[4]) ** 2 + (n[3] + n[6]) ** 2
    I5 = (n[5] - 3 * n[4]) * (n[5] + n[4]) * ((n[5] + n[4]) ** 2 - 3 * (n[3] + n[6])) + (3 * n[3] - n[6]) * (
        n[3] + n[6]) * (3 * (n[5] + n[4]) ** 2 - (n[3] + n[6]) ** 2)
    I6 = (n[1] - n[2]) * ((n[5] + n[4]) ** 2 - (n[3] + n[6]) ** 2) + 4 * n[0] * (n[5] + n[4]) * (n[3] + n[6])
    I7 = (3 * n[3] - n[6]) * (n[5] + n[4]) * ((n[5] + n[4]) ** 2 - 3 * (n[3] + n[6])) - (n[5] - 3 * n[4]) * (
        n[3] + n[6]) * (3 * (n[5] + n[4]) ** 2 - (n[3] + n[6]) ** 2)
    return [I1, I2, I3, I4, I5, I6, I7]


@autoCreateAttribute("huInvariant", 0)
def addAttributeHuInvariant(tree, attribute, momentAttributeName="moments"):
    """
    Compute the Hu invariants (1 to 7) of each component
    """
    attrMoment = tree.getAttribute(momentAttributeName)

    for i in tree.iteratorFromPixelsToRoot():
        attribute[i] = computeHuInvariant(computeScaleInvariantMoments2d(computeCentralMoments2d(attrMoment[i])))

@autoCreateAttribute("inertia", 0)
def addAttributeInertia2d(tree, attribute, momentAttributeName="moments"):
    """
    Compute the moment of inertia of each component 
    """
    attrMoment = tree.getAttribute(momentAttributeName)

    def computeInertia(m):
        xmean = m[1] / m[0]
        ymean = m[2] / m[0]
        u20 = m[4] - xmean * m[1]
        u02 = m[5] - ymean * m[2]
        return (u20 + u02) / (m[0] * m[0])

    for i in tree.iteratorFromPixelsToRoot():
        attribute[i] = computeInertia(attrMoment[i])


def addAttributeElongationOrientation2d(tree, nameElongation="elongation", nameOrientation="orientation",
                                        momentAttributeName="moments"):
    """
    Compute the elongation and orientation of each component
    """
    attrElongation, created1 = tree.addAttribute(nameElongation)
    attrOrientation, created2 = tree.addAttribute(nameOrientation)
    if not created1 and not created2:
        return attrElongation, attrOrientation
    attrMoment = tree.getAttribute(momentAttributeName)



    # Compute the elongation and the orientation from an array of moments [M00 M10 M01 M11 M20 M02]    
    def elongationOrientationFromMoments(m):
        sign = lambda x: copysign(1, x)
        xmean = m[1] / m[0]
        ymean = m[2] / m[0]
        xvar = m[4] / m[0] - xmean * xmean
        yvar = m[5] / m[0] - ymean * ymean
        xycovar = m[3] / m[0] - xmean * ymean
        if xvar - yvar != 0:
            a = 0.5 * atan(2 * xycovar / (xvar - yvar))
            if (sign(a) * sign(xycovar)) < 0.0:
                a += pi / 2.0
            alpha = a
        else:
            if xycovar == 0:
                alpha = 0
            else:
                alpha = pi / 4.0
        lambda1 = max(0, 0.5 * (xvar + yvar + sqrt(4 * xycovar * xycovar + (xvar - yvar) * (xvar - yvar))))
        lambda2 = max(0, 0.5 * (xvar + yvar - sqrt(4 * xycovar * xycovar + (xvar - yvar) * (xvar - yvar))))

        if lambda1 == 0:  # ill posed case
            el = 1;
        elif lambda2 == 0:  # ill posed case
            el = sqrt((lambda2 + 1) / (lambda1 + 1));
        else:
            el = sqrt(lambda2 / lambda1);

        return el, alpha

    for i in tree.iteratorFromPixelsToRoot():
        attrElongation[i], attrOrientation[i] = elongationOrientationFromMoments(attrMoment[i])

    return attrElongation, attrOrientation


@autoCreateAttribute("levelStatistics", None)
def addAttributeLevelStatistics(tree, attribute, levelImage=None):
    """
    [Mean Variance]
    :param tree:
    :param attribute:
    :param levelImage:
    :return:
    """
    if not levelImage:
        levelImage = tree.level
    children = addAttributeChildren(tree)
    area = addAttributeArea(tree)
    meanSquare = attribute.getCopy(False)

    if type(levelImage[0]) is list:
        dim = len(levelImage[0])
        for i in tree.iteratorOnLeaves():
            attribute[i] = [levelImage[i], [0] * dim]
            meanSquare[i] = VMath.multV(levelImage[i], levelImage[i])
        for i in tree.iteratorFromPixelsToRoot(includePixels=False):
            mean = [0] * dim
            mean2 = [0] * dim
            for child in children[i]:
                mean = VMath.addV(mean, VMath.multS(attribute[child][0], area[child]))
                mean2 = VMath.addV(mean2, VMath.multS(meanSquare[child], area[child]))

            mean = VMath.divS(mean, area[i])
            mean2 = VMath.divS(mean2, area[i])
            meanSquare[i] = mean2
            attribute[i] = [mean, VMath.diffV(mean2, VMath.multV(mean, mean))]
    else:
        for i in tree.iteratorOnLeaves():
            attribute[i] = [levelImage[i], 0]
            meanSquare[i] = levelImage[i] ** 2
        for i in tree.iteratorFromPixelsToRoot(includePixels=False):
            mean = 0
            mean2 = 0
            for child in children[i]:
                mean += attribute[child][0] * area[child]
                mean2 += meanSquare[child] * area[child]
            mean /= area[i]
            mean2 /= area[i]

            meanSquare[i] = mean2
            attribute[i] = [mean, mean2 - mean * mean]


@autoCreateAttribute("depth", 0)
def addAttributeDepth(tree, attribute):
    for j in tree.iteratorFromRootToPixels():
        par = tree[j]
        if par != -1:
            attribute[j] = attribute[par] + 1


@autoCreateAttribute("highest", None)
def addAttributeHighest(tree, attribute):
    children = addAttributeChildren(tree)
    level = tree.level
    nbLeaves = tree.nbPixels
    for i in tree.iteratorFromPixelsToRoot():
        if i < nbLeaves:
            attribute[i] = level[i]
        else:
            maxValue = attribute[children[i][0]]
            for c in children[i]:
                if attribute[c] > maxValue:
                    maxValue = attribute[c]
            attribute[i] = maxValue


@autoCreateAttribute("height", 0)
def addAttributeHeight(tree, attribute):
    highest = addAttributeHighest(tree)
    level = tree.level
    attribute[-1] = highest[-1]
    for i in tree.iteratorFromPixelsToRoot(includeRoot=False):
        attribute[i] = highest[i] - level[tree[i]]


@autoCreateAttribute("lowest", None)
def addAttributeLowest(tree, attribute):
    children = addAttributeChildren(tree)
    level = tree.level
    nbLeaves = tree.nbPixels
    for i in tree.iteratorFromPixelsToRoot():
        if i < nbLeaves:
            attribute[i] = level[i]
        else:
            minValue = attribute[children[i][0]]
            for c in children[i]:
                if attribute[c] < minValue:
                    minValue = attribute[c]
            attribute[i] = minValue


@autoCreateAttribute("dynamics", None)
def addAttributeDynamics(tree, attribute, extremaAttributeName="highest"):
    extrema = tree.getAttribute(extremaAttributeName)
    if not extrema:
        raise HiPy.Structures.HiPyException("Attribute '" + extremaAttributeName +
                                            "' cannot be found in order to compute dynamics.")
    level = tree.level
    for i in tree.iteratorFromRootToPixels():
        parent = tree[i]
        if parent == -1:
            attribute[i] = abs(extrema[i])
        else:
            if extrema[i] == extrema[parent]:
                attribute[i] = attribute[parent]
            else:
                attribute[i] = abs(level[parent] - extrema[i])


@autoCreateAttribute("perimeter", 0)
def addAttributePerimeterComponentTree(tree, attribute, baseAdjacency=None):
    children = addAttributeChildren(tree)
    if baseAdjacency is None:
        baseAdjacency = tree.leavesAdjacency
    nbLeaves = tree.nbPixels
    visited = HiPy.Structures.Image(nbLeaves, False)
    for i in range(nbLeaves):
        attribute[i] = baseAdjacency.countOutEdges(i, includeExternal=True)
    for i in tree.iteratorFromPixelsToRoot(False):
        remove = 0
        for c in children[i]:
            attribute[i] += attribute[c]
            if c < nbLeaves:
                for n in baseAdjacency.getNeighbours(c):
                    if visited[n]:
                        remove += 2
                visited[c] = True
        attribute[i] -= remove


@autoCreateAttribute("perimeter", 0)
def addAttributePerimeterPartitionHierarchy(tree, attribute, baseAdjacency, includeExternal=True):
    frontierLength = addAttributeFrontierLengthPartitionHierarchy(tree, baseAdjacency)
    for i in tree.iteratorOnPixels():
        attribute[i] = baseAdjacency.countOutEdges(i, includeExternal=includeExternal)
    for i in tree.iteratorFromPixelsToRoot(includeRoot=False):
        attribute[i] -= 2 * frontierLength[i]
        attribute[tree[i]] += attribute[i]
    attribute[len(tree) - 1] -= 2 * frontierLength[len(tree) - 1]


@autoCreateAttribute("frontierLength", 0)
def addAttributeFrontierLengthPartitionHierarchy(tree, attribute, adjacency):
    nodeMapping = HiPy.Hierarchies.WatershedHierarchy.computeSaliencyMap(tree, adjacency, lambda i: i)
    for i in range(nodeMapping.nbPoints):
        for edge in nodeMapping.getOutEdges(i):
            if edge[1] > edge[0]:
                attribute[edge[2]] += 1


@autoCreateAttribute("compactness", 0)
def addAttributeCompactness(tree, attribute):
    area = addAttributeArea(tree)
    perimeter = addAttributePerimeterComponentTree(tree)
    for i in tree.iteratorFromPixelsToRoot():
        attribute[i] = 4.0 * pi * area[i] / (perimeter[i] * perimeter[i])


@autoCreateAttribute("complexity", 0)
def addAttributeComplexity(tree, attribute):
    area = addAttributeArea(tree)
    perimeter = addAttributePerimeterComponentTree(tree)
    for i in tree.iteratorFromPixelsToRoot():
        attribute[i] = perimeter[i] / area[i]


@autoCreateAttribute("extrema", True)
def addAttributeExtrema(tree, attribute):
    """ 
    true if node is a maxima, false otherwise
    """
    for i in tree.iteratorFromPixelsToRoot(False):
        par = tree[i]
        if par != -1:
            attribute[par] = False


@autoCreateAttribute("isLeaf", True)
def addAttributeIsLeaf(tree, attribute):
    """ 
    True if node is a logical leaves, false otherwise
    
    In a partition hierarchy, logical leaves and structure leaves are the same (pixels are part of the hierarchy)
    
    In a component tree, logical leaves are nodes such that every child is a structure leaf (pixels are not part of the hierarchy)
    """
    if tree.treeType == HiPy.Structures.TreeType.PartitionHierarchy:
        for i in tree.iteratorFromPixelsToRoot(False):
            attribute[i] = False
    elif tree.treeType == HiPy.Structures.TreeType.ComponentTree:
        for i in range(tree.nbPixels):
            attribute[i] = False
        for i in range(tree.nbPixels, len(tree)):
            par = tree[i]
            if par != -1:
                attribute[par] = False
    else:
        raise Exception("addAttributeIsLeaf: unknown tree type.")


@autoCreateAttribute("rank", 0)
def addAttributeRank(tree, attribute):
    """
    Indicates the merging order. Leaves have rank 0.

    Target: Binary Partition Tree 
    """
    nbPixels = tree.nbPixels
    level = tree.level

    for i in tree.iteratorFromPixelsToRoot(False):
        attribute[i] = i - nbPixels + 1

    for i in tree.iteratorFromPixelsToRoot(False, False):
        par = tree[i]
        if level[i] == level[par]:
            attribute[par] = attribute[i]


@autoCreateAttribute("hog", None)
def addAttributeHOG(tree, attribute, grayImage, orientationBins=8):
    xGradient = simpleXGradient(grayImage)
    yGradient = simpleYGradient(grayImage)
    children = addAttributeChildren(tree)
    magnitude = grayImage.getCopy(False)
    orientation = grayImage.getCopy(False)
    binWidth = pi / orientationBins

    def getOrientationBin(angleOri):
        if angleOri < 0: #unsigned
            angle = pi + angleOri
        else:
            angle = angleOri
        binNum = int(floor(angle/binWidth))
        if binNum == orientationBins:
            binNum -= 1

        return binNum

    for i in grayImage.iterateOnPixels():
        magnitude[i] = sqrt(xGradient[i]*xGradient[i] + yGradient[i]*yGradient[i])
        orientation[i] = atan2(yGradient[i], xGradient[i])

    for i in attribute.iterateOnPixels():
        attribute[i] = [0]*orientationBins

    for i in tree.iteratorOnPixels():
        attribute[i][getOrientationBin(orientation[i])] = magnitude[i]

    for i in tree.iteratorFromPixelsToRoot(includePixels=False):
        for c in children[i]:
            attribute[i] = VMath.addV(attribute[i], attribute[c])

    for i in attribute.iterateOnPixels():
        attribute[i] = VMath.divS(attribute[i], VMath.normEpsilon(attribute[i]))
