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
The "Structure" module contains the definitions of the fundamental data structures
used in the HiPy library.

Image is the base class to store data.
An Image is basically an integer indexed list with utility functions.

An image may be associated with an Adjacency that describes the neighbouring relations among its
pixels.
The classes of the Adjacency hierarchy are:
    * AbstractAdjacency: the base class of all adjacencies.
    * AbstractWeightedAdjacency: mixin between AbstactAdjacency and Image, as a valued object it \
    can be processed as an image!
    * WeightedAdjacency: extends AbstractWeightedAdjacency, describes a non-directed weighted \
    adjacency (list based representation, suitable only for sparse graphs).
    * DirectedWeightedAdjacency: extends AbstractWeightedAdjacency, describes a directed weighted \
    adjacency (list based representation, suitable only for sparse graphs).
    * AdjacencyNdRegular: extends AbstractAdjacency, implicit representation of a translation \
    invariant adjacency
    * AdjacencyTree: extends AbstractAdjacency, wraps a tree into an Adjacency

An image may be associated with an Embedding that maps its linear integer coordinate
to another space.
The classes of the Embedding hierarchy are:
    * AbstractEmbedding: the base class of all Embeddings.
    * Embedding2dGrid: maps the integers to 2d pixel coordinates.

Created on 3 juin 2015

@author: perretb
"""
import logging
import copy
import HiPy.Processing.Attributes
from HiPy.Util import VMath
from heapq import heappop, heappush
from enum import Enum
from collections import deque


class HiPyException(Exception):
    """Base class for exceptions in this module."""
    pass


def initLogger():
    """
    Creates a simple logger with a console handler for the module
    """
    logger = logging.getLogger("HiPyLogger")
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.removeHandler(handler)
    logger.addHandler(handler)
    return logger


HiPyLogger = initLogger()


class TreeType(Enum):
    """
    This enum lists the supported tree types: component tree and Partition hierarchy.
    """
    ComponentTree = 1
    PartitionHierarchy = 2


def isImmutable(var):
    return type(var) in (str, int, bool, float, tuple)


class Image(list):
    """
    A basic image represented as a linear array of values (can be anything).

    An image may be associated with a (valued) adjacency relation (sounds like a graph)

    An image may be associated with an embedding, mapping its linear space to another target space
    """

    def __init__(self, length, initValue=0, adjacency=None, embedding=None):
        """
        Creates a new image of "length" elements initialized to "initValue" (deep copied).

        It is strongly recommended to provide an adjacency and an embedding as many other
        functions will expect to find them.

        If initValue is identified as a possibly mutable type, it is deep copied in every location
        """
        super().__init__()
        if isImmutable(initValue):
            for _ in range(length):
                self.append(initValue)
        else:
            for _ in range(length):
                self.append(copy.deepcopy(initValue))
        self.embedding = embedding
        self.adjacency = adjacency

    def setAll(self, image):
        """
        Copy values of image into current image.

        If image[0] is identified as possibly mutable, then every value is deep copied.
        """
        if isImmutable(image[0]):
            for i in range(len(image)):
                self[i] = image[i]
        else:
            for i in range(len(image)):
                self[i] = copy.deepcopy(image[i])

    def getCopy(self, copyData=False):
        """
        Warning: this is not a true deep copy!
        Image data are not shared but adjacency and embedding are shared

        If copyData = True and self[0] is identified as a possible mutable type, then values are deep copied
        """
        nIm = Image(len(self), 0, self.adjacency, self.embedding)

        if copyData:
            if isImmutable(self[0]):
                for i in range(len(self)):
                    nIm[i] = self[i]
            else:
                for i in range(len(self)):
                    nIm[i] = copy.deepcopy(self[i])
        return nIm

    def getPixel(self, i):
        """
        Return the value of pixel i
        """
        return self[i]

    def setPixel(self, i, val):
        """
        Set the value of pixel i to val
        """
        self[i] = val

    def setPixelWCS(self, val, *args):
        """
        Set the value of the pixel of coordinate args in the image embedding to val
        """
        self[self.embedding.getLinearCoordinate(*args)] = val

    def getPixelWCS(self, *args):
        """
        Return the value of the pixel of coordinate args in the image embedding
        """
        return self[self.embedding.getLinearCoordinate(*args)]

    def getNeighbours(self, i):
        """
        Return the neighbours (direct succesors and predecors) of pixel i in the image adjacency
        """
        return self.adjacency.getNeighbours(i)

    def equals(self, image, equalFunction=lambda x, y: x == y):
        """
        Test if the current image is equal to the given image.

        Two images are equal if they have the same length and their pixel values are equal
        w.r.t. the given equalFunction (== by default)
        """
        if len(self) != len(image):
            return False
        for i in range(len(image)):
            if not equalFunction(self[i], image[i]):
                return False
        return True

    def iterateOnPixels(self):
        """
        Iterator on all pixels
        """
        return range(len(self))

    def addAttribute(self, name: str, defaultValue: "void*"=None, resetIfExist: bool=False) -> ("Image", bool):
        """
        Creates a new attribute image called with the given name and initialized with the given default value.
        The attribute image has the same size as the current image. The attribute is automatically added to
        the members of the current instance.

        If an attribute with the same name already exists nothing is done unless resetIfExist is True.

        In all cases the attribute with the given name is returned. The method also returns if a new attribute image
        has been created.
        :param name: name of the attribute to create
        :param defaultValue: the initialisation value of the attribute image
        :param resetIfExist: indicates is a new attribute image must be created is an attribute with the same name \
            already exists.
        :return: (attribute image, new attribute image created?)
        """
        if name not in self.__dict__ or resetIfExist:
            image = Image(len(self), defaultValue, self.adjacency, self.embedding)
            self.__dict__[name] = image
            return image, True

        return self.__dict__[name], False

    def deleteAttribute(self, name):
        """
        Remove the attribute "name"
        """
        if name in self.__dict__:
            del self.__dict__[name]

    def getAttribute(self, name):
        """
        Return the attribute "name"
        """
        if name in self.__dict__:
            return self.__dict__[name]
        else:
            return None


class AbstractEmbedding(object):
    """
    Abstract embedding: ie a bijective mapping from linear indices
    (pixels in the abstract Image space) to any target space
    """

    def getLinearCoordinate(self, *args):
        """
        Convert coordinates from the target space to the abstract linear Image space
        """
        raise NotImplementedError("Unsuported method" + " getLinearCoordinate")

    def fromLinearCoordinate(self, i):
        """
        Convert coordinates from the abstract linear Image space to the target space
        """
        raise NotImplementedError("Unsuported method" + " fromLinearCoordinate")

    def isInBoundsWCS(self, *args):
        """
        Test if the given point is in the bounds of the embedding
        """
        raise NotImplementedError("Unsuported method" + " isInBoundsWCS")

    def isInBoundsLinear(self, i):
        """
        Test if the given point (in abstract linear image space) is in the bounds of the embedding
        """
        return self.isInBoundsWCS(self.fromLinearCoordinate(i))


class Embedding2dGrid(AbstractEmbedding):
    """
    This class represents an embedding from discrete linear coordinates to discrete 2d coordinates
    organized on a regular grid given by its width and height (typically a pixel grid).
    """

    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.size = [width, height]

    def getLinearCoordinate(self, *args):
        return args[1] * self.width + args[0]

    def fromLinearCoordinate(self, i):
        return i % self.width, i // self.width

    def isInBoundsWCS(self, *args):
        return 0 <= args[0] < self.width and 0 <= args[1] < self.height


class AbstractAdjacency(object):
    """
    Abstract adjacency relation
    """

    def __init__(self, nbPoints):
        self.nbPoints = nbPoints

    def countEdges(self):
        """
        Count the number of edges in the adjacency relation.

        Not that in a non-directed adjacency the edge (i,j) and
        the edge (j,i) are equivalent and count for one edge.
        """
        count = 0
        for i in range(self.nbPoints):
            count += len(self.getSuccessors(i))
        return count

    def countOutEdges(self, i, includeExternal=False):
        """
        Count the number of edges of source vertex i.

        If includeExternal is True, it also count edges that are going out of the domain
        :param i: vertex index
        :param includeExternal: include edges going out of the definition domain?
        :return: number of out edges for vertex i
        """
        raise NotImplementedError("Unsupported method" + " countOutEdges")

    def getSuccessors(self, i):
        """
        :param i: vertex index
        :return list of points j such that i->j
        """
        raise NotImplementedError("Unsupported method" + " getSuccessors")

    def getPredecessors(self, i):
        """
        :param i: vertex index
        :return list of points j such that j->i
        """
        raise NotImplementedError("Unsupported method" + " getPredecessors")

    def getNeighbours(self, i):
        """
        :param i: vertex index
        :return list of points j such that j->i or i->j (Successors U Predecessors)
        """
        raise NotImplementedError("Unsupported method" + " getNeighbours")

    def getEdges(self, i):
        """
        :param i: vertex index
        :return list of in  edges of the form [j,i,*] or [i,j,*] such that j->i or i->j
        (* can be any auxiliary data, usually weights)

        """
        return self.getInEdges(i) + self.getOutEdges(i)

    def getOutEdges(self, i):
        """
        :param i: vertex index
        :return list of outHead edges of the form [i,j,*] such that i->j
        (* can be any auxiliary data, usually weights)

        """
        raise NotImplementedError("Unsupported method" + " getOutEdges")

    def getInEdges(self, i):
        """
        :param i: vertex index
        :return list of in  edges of the form [j,i,*] such that j->i
        (* can be any auxiliary data, usually weights)
        """
        raise NotImplementedError("Unsupported method" + " getInEdges")


class AdjacencyNdRegular(AbstractAdjacency):
    """
    Implicit representation of a shift-invariant adjacency relation
    in the n dimensional regular grid
    """

    def __init__(self, embedding, neighbourList, weights=None):
        """
        :param embedding: Defines the adjacency on the domain given by the embedding.

        :param neighbourList: The neighbouring relation is given by the neighbour list of the point 0.
        EG. a 4-adjacency in 2d is given by the neighbour list [ (0,-1), (-1,0), (1,0), (0,1)]

        :param weights: Weights should be a list of same length as neighbour list giving the weight
        of the corresponding edge.
        If no weights are provided, all edges are assigned a wight of 1.
        """
        super().__init__(VMath.mult(embedding.size))
        self.embedding = embedding
        self.neighbourList = neighbourList
        self.nbNeighbours = len(neighbourList)
        self.weights = weights if weights is not None else [1] * self.nbNeighbours

    def countEdges(self, includeExternal=False):
        count = 0
        isInBounds = self.embedding.isInBoundsWCS
        linear = self.embedding.getLinearCoordinate
        for i in range(self.nbPoints):
            coordi = self.embedding.fromLinearCoordinate(i)
            for neighbour in self.neighbourList:
                coordn = VMath.addV(coordi, neighbour)
                linearn = linear(*coordn)
                if linearn >= i and isInBounds(*coordn):
                    count += 1
        return count

    def countOutEdges(self, i, includeExternal=False):
        if includeExternal:
            return self.nbNeighbours
        else:
            count = 0
            coordi = self.embedding.fromLinearCoordinate(i)
            isInBounds = self.embedding.isInBoundsWCS
            for neighbour in self.neighbourList:
                coordNeighbour = VMath.addV(coordi, neighbour)
                if isInBounds(*coordNeighbour):
                    count += 1
            return count

    def getNeighbours(self, i):
        coordi = self.embedding.fromLinearCoordinate(i)
        neighbourList = []
        isInBounds = self.embedding.isInBoundsWCS
        linear = self.embedding.getLinearCoordinate
        for neighbour in self.neighbourList:
            coordNeighbour = VMath.addV(coordi, neighbour)
            if isInBounds(*coordNeighbour):
                neighbourList.append(linear(*coordNeighbour))
        return neighbourList

    def getPredecessors(self, i):
        return self.getNeighbours(i)

    def getSuccessors(self, i):
        return self.getNeighbours(i)

    def getOutEdges(self, i):
        coordi = self.embedding.fromLinearCoordinate(i)
        neighbourList = []
        isInBounds = self.embedding.isInBoundsWCS
        linear = self.embedding.getLinearCoordinate
        weights = self.weights
        for i in range(self.nbNeighbours):
            neighbour = self.neighbourList[i]
            coordNeighbour = VMath.addV(coordi, neighbour)
            if isInBounds(*coordNeighbour):
                neighbourList.append([neighbour, linear(*coordNeighbour), weights[i]])
        return neighbourList

    def getInEdges(self, i):
        coordi = self.embedding.fromLinearCoordinate(i)
        neighbourList = []
        isInBounds = self.embedding.isInBoundsWCS
        linear = self.embedding.getLinearCoordinate
        weights = self.weights
        for i in range(self.nbNeighbours):
            neighbour = self.neighbourList[i]
            coordNeighbour = VMath.addV(coordi, neighbour)
            if isInBounds(*coordNeighbour):
                neighbourList.append([linear(*coordNeighbour), neighbour, weights[i]])
        return neighbourList

    @staticmethod
    def getAdjacency2d4(size):
        """
        Creates an instance of AdjacencyNdRegular representing a classical 4 adjacency
        :param size: (width, height) of the underlying 2d grid
        :returns AdjacencyNdRegular representing 4 adjacency
        """
        return AdjacencyNdRegular(Embedding2dGrid(size[0], size[1]),
                                  [(0, -1), (-1, 0), (1, 0), (0, 1)])

    @staticmethod
    def getAdjacency2d8(size):
        """
        Creates an instance of AdjacencyNdRegular representing a classical 8 adjacency
        :param size: (width, height) of the underlying 2d grid
        :returns AdjacencyNdRegular representing 8 adjacency
        """
        return AdjacencyNdRegular(Embedding2dGrid(size[0], size[1]),
                                  [(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)])


class AdjacencyTree(AbstractAdjacency):
    """
    Represent the adjacency relation in a tree.

    The adjacency is oriented from the root toward the leaves
    """
    # TODO: mixin in Tree

    def __init__(self, tree):
        super(AdjacencyTree, self).__init__(len(tree))
        self.tree = tree

    def getNeighbours(self, i):
        neighbourList = []
        if self.tree[i] != -1:
            neighbourList.append(self.tree[i])
        neighbourList.extend(self.tree.children[i])
        return neighbourList

    def getPredecessors(self, i):
        if self.tree[i] != -1:
            return [self.tree[i]]
        return []

    def getSuccessors(self, i):
        return self.tree.children[i]

    def getOutEdges(self, i):
        neighbourList = []
        for child in self.tree.children[i]:
            neighbourList.append([i, child])
        return neighbourList

    def getInEdges(self, i):
        if self.tree[i] != -1:
            return [self.tree[i], i]
        return []

    # @todo: refactor
    def __getNeigbourPointsExtended__(self, i, maxDist, dist=(lambda x, y: 1)):
        neighbourList = [i]
        tree = self.tree
        children = self.tree.children

        cur = i
        par = tree[cur]
        distr = 0
        while par != -1:
            distr = distr + dist(cur, par)
            if distr <= maxDist:
                neighbourList.append(par)
                cur = par
                par = tree[par]
            else:
                par = -1


        # search towards leaves
        stack = [[i, 0]]
        while len(stack) != 0:
            elem = stack.pop()
            cur = elem[0]
            diste = elem[1]
            for child in children[cur]:
                distd = diste + dist(cur, child)
                if distd <= maxDist:
                    neighbourList.append(child)
                    stack.append([child, distd])

        return neighbourList


class AbstractWeightedAdjacency(AbstractAdjacency, Image):
    """
    An abstract edge weighted adjacency.

    A weighted adjacency is an image of weights with two attributes: source and target.
    Each pixel i of the image I represent an edge from source[i] to target[i] weighted by I[i].

     Known concrete implementations are
      - WeightedAdjacency
      - DirectedWeightedAdjacency
    """

    def __init__(self, size: int) -> None:
        """
            Creates a new empty adjacency on a set of size elements.
            :param size: number of elements
            """
        AbstractAdjacency.__init__(self, size)
        Image.__init__(self, 0)
        self.addAttribute('source')
        self.addAttribute('target')


class WeightedAdjacency(AbstractWeightedAdjacency):
    """
    Non directed adjacency represented with list of edges.

    Note:
      - methods getSuccessors, getPredecessors, and getNeighbours are equivalent
      - methods getEdges, getOutEdges, and getInEdges are equivalent
    """

    def __init__(self, size, optimizedRemove=False):
        """
        Create a new empty adjacency on a set of size elements
        """
        AbstractWeightedAdjacency.__init__(self, size)
        self.optimizedRemove = optimizedRemove
        self.edgeList = []
        if optimizedRemove:
            self.appendFunct = WeightedAdjacency.addSet#lambda collection, item: collection.add(item)
            for _ in range(size):
                self.edgeList.append(set())
        else:
            self.appendFunct = WeightedAdjacency.addList#lambda collection, item: collection.append(item)
            for _ in range(size):
                self.edgeList.append([])

    @staticmethod
    def addList(collection, item):
        collection.append(item)

    @staticmethod
    def addSet(collection, item):
        collection.add(item)

    def addVertex(self):
        newVertex = len(self.edgeList)
        if self.optimizedRemove:
            self.edgeList.append(set())
        else:
            self.edgeList.append([])
        return newVertex

    def countEdges(self):
        return len(self)

    def countOutEdges(self, i, includeExternal=False):
        return len(self.edgeList[i])

    def createEdge(self, source, target, weight=1):
        """
        Create a new edge in the graph, linking node source to node dest.

        As the graph is undirected, by convention, the method will ensure that source < target

        Warning: does not verify is the edge already exists

        Returns the index of the new edge
        """
        i = len(self)
        if source > target:
            source, target = target, source
        self.append(weight)
        self.source.append(source)
        self.target.append(target)

        self.appendFunct(self.edgeList[source],i)
        if source != target:
            self.appendFunct(self.edgeList[target],i)
        return i

    def setEdge(self, index, source, target, weight):
        if source > target:
            source, target = target, source
        prevSource = self.source[index]
        prevTarget = self.target[index]

        cond1 = prevSource != source
        cond2 = prevSource != target
        cond3 = prevTarget != target
        cond4 = prevTarget != source

        if cond1:
            if cond2:
                self.edgeList[prevSource].remove(index)
            if cond4:
                self.appendFunct(self.edgeList[source],index)
        if cond3:
            if cond4:
                self.edgeList[prevTarget].remove(index)
            if cond2:
                self.appendFunct(self.edgeList[target],index)

        self.source[index] = source
        self.target[index] = target
        self[index] = weight

    def removeEdge(self, i):
        source = self.source[i]
        target = self.target[i]
        self.edgeList[source].remove(i)
        if source != target:
            self.edgeList[target].remove(i)

    def getEdgeFromIndex(self, i):
        return self.source[i], self.target[i], self[i]

    def getCopy(self, copyData=True, optimizedRemove=False):
        """
        Returns a copy of the current adjacency.
        Guaranties that edges indices are consistent between the copy and the current object
        (ie copy[i]==original[i] for all i)

        Parameter copyData is silently ignored
        """
        adj = WeightedAdjacency(self.nbPoints, optimizedRemove=optimizedRemove)
        for i in range(len(self)):
            adj.createEdge(self.source[i], self.target[i], self[i])
        return adj

    def getSuccessors(self, i):
        # ugly hack to symmetries the adjacency on the fly
        return (self.source[e] + self.target[e] - i for e in self.edgeList[i])

    def getPredecessors(self, i):
        return self.getSuccessors(i)

    def getNeighbours(self, i):
        return self.getSuccessors(i)

    def getEdges(self, i):
        return ((self.source[e], self.target[e], self[e]) for e in self.edgeList[i])

    def getOutEdges(self, i):
        return ((i, self.source[e] + self.target[e] - i, self[e]) for e in self.edgeList[i])

    def getInEdges(self, i):
        return ((self.source[e] + self.target[e] - i, i, self[e]) for e in self.edgeList[i])

    @staticmethod
    def createAdjacency(baseAdjacency, weightingFunction=None):
        """
        Create a new adjacency equivalent to the given adjacency but
        with a different weighting function.

        Warning: the base adjacency is assumed to be symmetric!

        Typical use is to transform an implicit k-adjacency into an explicit weighted adjacency.
        """
        adj = WeightedAdjacency(baseAdjacency.nbPoints)
        if weightingFunction is not None:
            for i in range(adj.nbPoints):
                for j in baseAdjacency.getSuccessors(i):
                    if j > i:
                        adj.createEdge(i, j, weightingFunction(i, j))
        else:
            for i in range(adj.nbPoints):
                for j in baseAdjacency.getSuccessors(i):
                    if j > i:
                        adj.createEdge(i, j)

        return adj

    @staticmethod
    def createKHopAjacency(baseAdjacency, kValue, weightingFunction=lambda i, j, n: n):
        """
        Construct an explicit non directed K-hop adjacency from a base (non directed adjacency).

        A weighting function f may be provided (given two nodes i, j , f(i,j,k) will be
        the weight of the edge linking i to j in k hops).

        """
        nbPoints = baseAdjacency.nbPoints
        newAdj = WeightedAdjacency(nbPoints)
        flags = [None] * nbPoints

        for i in range(nbPoints):
            heap = []
            heappush(heap, (0, i))
            flags[i] = i
            while heap:
                (k, n) = heappop(heap)
                if n > i:
                    newAdj.createEdge(i, n, weightingFunction(i, n, k))

                if k < kValue:
                    k += 1
                    for j in baseAdjacency.getSuccessors(n):
                        if flags[j] != i:
                            heappush(heap, (k, j))
                            flags[j] = i

        return newAdj

    @staticmethod
    def createLineGraph(baseAdjacency, weightingFunction=lambda i, j: 1):
        """
        Create the line graph of the given non directed base adjacency
        ie. a new adjacency relation were the points are the edges of the base adjacency
        and two points are adjacent is their corresponding edges in the base adjacency
        share an extremity.

        The base adjacency must be  a WeightedAdjacency (Use WeightedAdjacency.createAdjacency
        in order to convert any non-directed adjacency to a WeightedAdjaceny).
        The method insure the consistency between the edge indices in the base adjacency and
        the points in the new adjacency.

        The given weighting function computes the weight of the edge linking the edge i to
        the edge j (of the base adjacency).
        """
        nbPoints = baseAdjacency.countEdges()
        newAdj = WeightedAdjacency(nbPoints)
        source = baseAdjacency.source
        target = baseAdjacency.target
        baseEdgeList = baseAdjacency.edgeList
        for i in range(nbPoints):
            for vertex in [source[i], target[i]]:
                for edge in baseEdgeList[vertex]:
                    if edge > i:
                        newAdj.createEdge(i, edge, weightingFunction(i, edge))
        return newAdj

    @staticmethod
    def createReflexiveRelation(baseAdjacency, weightingFunction=lambda i: 1):
        """
        Create a new Weighted Adjacency equals to the reflexive closure of
        the given base adjacency.

        The given weighting function computes the weight of the reflexive edges and
        takes a single parameter: the element i where the edge (i,i) is created.
        """
        newAdj = WeightedAdjacency(baseAdjacency.nbPoints)
        for i in range(newAdj.nbPoints):
            for e in baseAdjacency.getOutEdges(i):
                if e[1] > i:
                    newAdj.createEdge(*e)

            newAdj.createEdge(i, i, weightingFunction(i))

        return newAdj


class DirectedWeightedAdjacency(AbstractWeightedAdjacency):
    """
    Directed adjacency represented using doubly linked lists of out edges.

    The following operations are done in O(1) constant time:
      - Edge insertion
      - Edge deletion
      - Fusion of edge lists
      - Getting the list of out edges for a given source node

    The following operations are not supported (their implementation would be very innefficient):
      - getPredecessors
      - getNeighbours
    However, the method getTranspose which will transpose the current graph gives
    an easy solution to get those functions if needed.
    """

    def __init__(self, size) -> None:
        """
             Create a new empty adjacency on a set of size elements
            """
        AbstractWeightedAdjacency.__init__(self, size)
        self.outHead = [-1] * size
        self.outTail = [-1] * size
        self.prevEdge = []
        self.nextEdge = []

    @staticmethod
    def createAdjacency(baseAdjacency, weightingFunction=None):
        """
        Create a new adjacency equivalent to the given adjacency but with
        a different weighting function.

        Typical use is to transform an implicit k-adjacency into an explicit weighted adjacency.
        """
        adj = DirectedWeightedAdjacency(baseAdjacency.nbPoints)
        if weightingFunction is not None:
            for i in range(adj.nbPoints):
                for j in baseAdjacency.getSuccessors(i):
                    adj.createEdge(i, j, weightingFunction(i, j))
        else:
            for i in range(adj.nbPoints):
                for j in baseAdjacency.getSuccessors(i):
                    adj.createEdge(i, j)

        return adj

    def countOutEdges(self, i, includeExternal=False):
        count = 0
        e = self.outHead[i]
        while e != -1:
            # dest is the adjacent vertex
            count += 1
            # next edge in the edge list
            e = self.nextEdge[e]
        return count

    def getSuccessors(self, i):
        nodes = []
        e = self.outHead[i]
        while e != -1:
            # dest is the adjacent vertex
            nodes.append(self.target[e])
            # next edge in the edge list
            e = self.nextEdge[e]
        return nodes

    def getPredecessors(self, i):
        raise Exception("Unsupported method" + " getPredecessors")

    def getNeighbours(self, i):
        raise Exception("Unsupported method" + " getNeighbours")

    def getOutEdges(self, i):
        edges = []
        e = self.outHead[i]
        while e != -1:
            # dest is the adjacent vertex
            edges.append([self.source[e], self.target[e], self[e]])
            # next edge in the edge list
            e = self.nextEdge[e]
        return edges

    def getInEdges(self, i):
        raise Exception("Unsupported method" + " getInEdges")

    def createEdge(self, source, target, weight=1):
        """
        Create a new edge in the graph, linking node source to node dest.
        """
        i = len(self)
        self.append(weight)
        self.source.append(source)
        self.target.append(target)
        self.prevEdge.append(-1)
        self.nextEdge.append(-1)
        self.setEdge(source, target, i)

    def setEdge(self, source, _, i):
        """
        Set the edge i of the graph, the edge is put at the head of the outHead linked list
        of the source vertex
        """
        if self.outHead[source] == -1:
            self.outTail[source] = i
        else:
            head = self.outHead[source]
            self.prevEdge[head] = i
            self.nextEdge[i] = head
        self.outHead[source] = i

    def concatEdgeOut(self, vertex1, vertex2):
        """
        Concatenate edge lists, outHead edges of vertex n2 are added to vertex n1
        """
        if self.outTail[vertex1] != -1:
            self.nextEdge[self.outTail[vertex1]] = self.outHead[vertex2]
        if self.outHead[vertex2] != -1:
            self.prevEdge[self.outHead[vertex2]] = self.outTail[vertex1]
        if self.outTail[vertex2] != -1:
            self.outTail[vertex1] = self.outTail[vertex2]
        self.outHead[vertex2] = self.outTail[vertex2] = -1

    def removeEdgeOut(self, n, e):
        """
        Remove the edge e of the list of outHead edges of vertex n
        """
        if self.outHead[n] == e:
            self.outHead[n] = self.nextEdge[e]
        if self.outTail[n] == e:
            self.outTail[n] = self.prevEdge[e]
        if self.nextEdge[e] != -1:
            self.prevEdge[self.nextEdge[e]] = self.prevEdge[e]
        if self.prevEdge[e] != -1:
            self.nextEdge[self.prevEdge[e]] = self.nextEdge[e]
        self.nextEdge[e] = self.prevEdge[e] = -1

    def getCopy(self, copyData=True):
        """
        Creates a copy of the current adjacency

        Parameter copyData is silently ignored
        """
        adj = DirectedWeightedAdjacency(self.nbPoints)
        for i in range(len(self)):
            adj.createEdge(self.source[i], self.target[i], self[i])
        return adj

    def getSymmetricMax(self):
        """
        Symmetries the graph with the max strategy: whenever an edge (p,q) is found
        the edge (q,p) is added to the graph
        """
        graph = DirectedWeightedAdjacency(self.nbPoints)
        for i in range(self.nbPoints):
            for j in self.getSuccessors(i):
                if not j in graph.getSuccessors(i):
                    graph.createEdge(i, j)
                if not i in self.getSuccessors(j):
                    graph.createEdge(j, i)
        return graph

    def getSymmetricMin(self):
        """
        Symmetrize the graph with the min strategy: an edge (p,q) is preserved only if
        the edge (q,p) is also in the graph
        """
        graph = DirectedWeightedAdjacency(self.nbPoints)
        for i in range(self.nbPoints):
            for j in self.getSuccessors(i):
                if j > i and i in self.getSuccessors(j):
                    graph.createEdge(i, j)
                    graph.createEdge(j, i)
        return graph

    def getTranspose(self):
        """
        Transpose the graph: each edge (p,q) is transformed into the edge (q,p)
        """
        graph = DirectedWeightedAdjacency(self.nbPoints)
        for i in range(self.nbPoints):
            for j in self.getSuccessors(i):
                graph.createEdge(j, i)
        return graph


class Tree(Image):
    """
    A tree is an image whose pixel values encode the parent relation.

    A tree always own at least two attributes:
      - "level" storing the level of each node
      - "reconstructedValue" storing the result of a filtering

    A tree has a treeType: component trees or partition hierarchies (see class TreeType).
    Also those two types of trees have a lot in common there is a difference at the leave level.

    The leaves of the data structures always represent the pixels or the original
    image. Those pixels are usualy not part of the component tree (contrarily to the
    partition hierarchy). However those pixels are kept in the data structure as they are useful in
    many processings.

    In consequence all the methods follow the following naming conventions:
        * "pixels" stands for the leaves of the data structure
        * "leaves" stands for the leaves of the represented tree
    If the tree type is "TreeType.PartitionHierachy" then those two notions are the same, if it is
    "TreeType.ComponentTree" then they are different. More precisely, in case of component tree,
    the leaves of the logical tree are the nodes whose children are all leaves in the
    data structure (ie. pixels).
    """

    def __init__(self, treeType, parent, levels, image=None):
        super(Tree, self).__init__(len(parent), None)
        HiPyLogger.debug("constructor: Tree")
        self.treeType = treeType
        self.adjacency = AdjacencyTree(self)
        self.setAll(parent)
        self.leavesAdjacency = image.adjacency if image is not None else None
        self.leavesEmbedding = image.embedding if image is not None else None
        self.nbPixels = len(image) if image is not None else Tree._countLeaves(parent)
        self.addAttribute("level")
        self.level.setAll(levels)
        self.addAttribute('reconstructedValue')
        HiPy.Processing.Attributes.addAttributeChildren(self)

    def getParent(self, i):
        """
        Return the parent of the node i, -1 if i is a root
        :return parent of the node i, -1 if i is a root
        """
        return self[i]

    def getRoot(self):
        """
        Return the index of the root node
        :return: index of the root node
        """
        return len(self) - 1

    def filterDirect(self, filteringCriterion):
        """
        Apply a given filtering criterion f on the tree.
        A filtering criterion is a function that associates the value true or false to a tree
        and a node.

        For each node n the attribute deleted is set to the value of filteringCriterion(self,n)
        """
        self.addAttribute("deleted", False)
        for i in self.iteratorFromPixelsToRoot():
            self.deleted[i] = filteringCriterion(self, i)

    def extractBranch(self, node, attributeName=None):
        """
        Extract values along a branch of the tree, from the given node toward the root.

        If attributeName==False, the result is the list of nodes from the given node to the root.

        Otherwise the result is the list of node attribute values from the given node to the root
        for the attribute "attributeName"
        """
        attr = None
        if attributeName is not None:
            attr = self.getAttribute(attributeName)

        res = []
        res.append(node if attr == None else attr[node])
        node = self[node]
        while node != -1:
            res.append(node if attr == None else attr[node])
            node = self[node]

        return res

    def reconstructImage(self, attributeName="level", criterion=None):
        """
        Reconstruct an image using the value of the attribute "attributeName".
        The final pixel value is given by the value of its closest selected (see below) parent.

        If a criterion is provided (a function that associates True or False to each node),
        then a node x is selected if criterion(x)==False

        If no criterion is provided and the attribute "deleted" is defined, a node x is selected
        if deleted[x]==False

        If no criterion is provided and the attribute "deleted" is not defined then every nodes
        are selected.
        """
        if not criterion:
            if "deleted" in self.__dict__:
                criterion = (lambda x: self.deleted[x])
            else:
                criterion = (lambda _: False)
        root = len(self) - 1
        for i in self.iteratorFromRootToPixels():
            if i >= self.nbPixels and (i == root or not criterion(i)):
                self.reconstructedValue[i] = self.__dict__[attributeName][i]
            else:
                self.reconstructedValue[i] = self.reconstructedValue[self[i]]

        image = Image(self.nbPixels, 0, self.leavesAdjacency, self.leavesEmbedding)
        for i in range(len(image)):
            image[i] = self.reconstructedValue[i]

        return image

    def iteratorFromLeavesToRoot(self, includeLeaves=True, includeRoot=True):
        """
        Provides an iterator on the nodes of the tree going from the leaves to the root.

        """
        return TreeIterator(self, True, False, includeLeaves, includeRoot)

    def iteratorFromRootToLeaves(self, includeLeaves=True, includeRoot=True):
        """
        Provides an iterator on the nodes of the tree going from the root to the leaves.
        """
        return TreeIterator(self, True, True, includeLeaves, includeRoot)

    def iteratorFromPixelsToRoot(self, includePixels=True, includeRoot=True):
        """
        Provides an iterator on the nodes of the tree going from the pixels to the root.

        """
        return TreeIterator(self, False, False, includePixels, includeRoot)

    def iteratorFromRootToPixels(self, includePixels=True, includeRoot=True):
        """
        Provides an iterator on the nodes of the tree going from the root to the pixels.
        """
        return TreeIterator(self, False, True, includePixels, includeRoot)

    def iteratorOnPixels(self):
        """
        Creates an iterator on the pixels of the tree
        """
        return range(self.nbPixels)

    def iteratorOnLeaves(self):
        """
        Creates an iterator on the leaves of the tree.
        """

        def leavesIterator():
            HiPy.Processing.Attributes.addAttributeIsLeaf(self)
            isLeaf = self.isLeaf
            i = self.nbPixels
            size = len(self)
            while i < size:
                while i < size and not isLeaf[i]:
                    i += 1
                if i < size:
                    yield i
                i += 1

        if self.treeType == TreeType.ComponentTree:
            return leavesIterator()
        else:
            return range(self.nbPixels)

    def lca(self, i, j):
        """
        Lowest common ancestor between nodes i and j
        Warning: Assume that the attribute depth is present !
        """
        depth = self.depth
        while i != j:
            if depth[i] > depth[j]:
                i = self[i]
            else:
                j = self[j]

        return i

    def nbNodes(self):
        """
        Count the number of logical nodes.

        For component trees, pixels that are not part of the logical structure must be substracted
        leading to len(tree)-tree.nbPixels
        For partition hierarchies or others the result is equal to len(tree)
        """
        if self.treeType == TreeType.ComponentTree:
            return len(self) - self.nbPixels
        else:
            return len(self)

    @staticmethod
    def _countLeaves(parent):
        """
        Count Leaves in a parent relation
        """
        leaves = [True] * len(parent)
        for i in range(len(parent)):
            par = parent[i]
            if par != -1:
                leaves[par] = False

        count = 0
        for i in range(len(parent)):
            if leaves[i]:
                count += 1

        return count

    @staticmethod
    def testTreeIsomorphism(tree1, tree2):
        """
        Test if tree1 and tree2 are isomorph assuming that leaves are ordered
        """
        if len(tree1) != len(tree2) or tree1.nbPixels != tree2.nbPixels:
            return False

        # both tree have same size so we need to find an injection m from the nodes of t1 to
        # the nodes of t2 such that for any node n of t1 m(parent(n))=parent(m(n))
        mapT1T2 = [-1] * len(tree1)

        for i in range(len(mapT1T2) - 1):  # root is root !
            # pixel mapping is constant
            if i < tree1.nbPixels:
                mapT1T2[i] = i
            # parent(n)
            pT1 = tree1[i]

            # parent(m(n))
            pT2 = tree2[mapT1T2[i]]
            if mapT1T2[pT1] == -1:
                mapT1T2[pT1] = pT2
            elif mapT1T2[pT1] != pT2:
                return False

        return True

    def simplifyTreeByCriterion(self, criterion, levelFunction) -> "Tree":
        """
        Creates a copy of the current Tree and deletes the nodes such that the criterion function is true.

        The level of the nodes in the new tree is given by the level function.

        The criterion is a function that takes two argument: a node index and its parent index and returns
        True (this node must be deleted), or False (do not delete node i). The parent index is the updated index (after
        node removal) and may differ from its parent index in the current tree.

        The level function is function that associates a value to a node index.

        :param criterion:
        :param levelFunction:
        :return:
        """
        nbLeaves = self.nbPixels
        nbNodes = len(self)
        children = HiPy.Processing.Attributes.addAttributeChildren(self)
        copyParent = self.getCopy(True)

        count = 0
        deleted = [False] * nbNodes
        deletedMap = [0] * nbNodes

        # from root to leaves, compute the new parent relation,
        # don't care of the holes in the parent tab
        for i in range(nbNodes - 2, nbLeaves - 1, -1):
            par = copyParent[i]
            if criterion(i, par):
                for c in children[i]:
                    copyParent[c] = par
                deleted[i] = True
                count += 1
            # inverse of what we want: number of deleted nodes after node i
            deletedMap[i] = count

        # correct the mapping
        for i in self.iteratorFromPixelsToRoot(False):
            deletedMap[i] = count - deletedMap[i]

        # new relations with correct size
        newParent = [-1] * (nbNodes - count)
        newLevel = [0] * (nbNodes - count)

        count = 0
        for i in range(0, nbNodes - 1):
            if not deleted[i]:
                par = copyParent[i]
                newPar = par - deletedMap[par]
                newParent[count] = newPar
                newLevel[count] = levelFunction(i)
                count += 1
        newParent[count] = -1
        newLevel[count] = levelFunction(self.getRoot())

        newTree = Tree(self.treeType, newParent, newLevel)
        newTree.leavesAdjacency = self.leavesAdjacency

        return newTree

    def simplifyTreeByAttribute(self, attributeName, levelName="level"):
        """
        Copy the current tree and delete the nodes n such that attribute[n]=attribute[self[n]]
        (and update the parent relation accordingly...

        The new node levels are taken from the attribute levelName
        """
        attribute = self.getAttribute(attributeName)
        level = self.getAttribute(levelName)
        criterion = lambda i, par: attribute[i] == attribute[par]
        levelFunction = lambda i: level[i]
        return self.simplifyTreeByCriterion(criterion, levelFunction)


class TreeIterator(object):
    """
    The TreeIterator class is intended to perform pre or post order traversal of trees.

    It manages several common cases as excluding or including leaves or the root.

    It also manages the logical difference between Component trees
    (where leaves represent pixels of the original image but not nodes of the logical tree)
    and partition hierarchies (where leaves are at the same time pixels of the original image
    and nodes of the logical tree).
    """

    def __init__(self, tree, logical=True, reverseOrder=False, includeLeaves=True,
                 includeRoot=True):
        self.tree = tree

        specialLogic = tree.treeType == TreeType.ComponentTree and logical
        if specialLogic and not includeLeaves:
            HiPy.Processing.Attributes.addAttributeIsLeaf(tree)
            self.nextMethod = self.nextLogical
        else:
            self.nextMethod = self.nextStructural

        if specialLogic or not includeLeaves:
            bmin = tree.nbPixels
        else:
            bmin = 0

        bmax = len(tree)
        if not includeRoot:
            bmax = bmax - 1

        if not reverseOrder:
            self.curVal = bmin
            self.limit = bmax
            self.step = 1
        else:
            self.curVal = bmax - 1
            self.limit = bmin - 1
            self.step = -1

    def __iter__(self):
        return self

    def nextLogical(self):
        """
        "next" iterator method when performing a logical browsing (only logical nodes are traversed)
        """
        tree = self.tree
        isLeaf = tree.isLeaf
        while self.curVal != self.limit and isLeaf[self.curVal]:
            self.curVal += self.step

        if self.curVal != self.limit:
            i = self.curVal
            self.curVal += self.step
            return i
        else:
            raise StopIteration()

    def nextStructural(self):
        """
        "next" iterator method when performing a structural browsing (every node is traversed)
        """
        if self.curVal != self.limit:
            i = self.curVal
            self.curVal += self.step
            return i
        else:
            raise StopIteration()

    def __next__(self):
        return self.nextMethod()


class PriorityQueue:
    """
    A simple multi-level priority queue.
    """

    def __init__(self, levels):
        """
        Create a new multi-level priority queue.
        :param levels: number of levels
        :return: void
        """
        self.data = []
        for _ in range(levels):
            self.data.append(deque())
        self.levels = levels
        self.size = 0

    def push(self, level, element):
        """
        Insert a new element at given priority level
        :param level: priority level
        :param element: new element
        :return: void
        """
        self.data[level].append(element)
        self.size += 1

    def isEmpty(self):
        """
        Test if the queue is empty
        :return: True is queue is empty, False otherwise
        """
        return self.size == 0

    def isLevelEmpty(self, level):
        """
        Test is particular priority level of the queue is empty
        :param level: priority level
        :return: True is the priority level of the queue is empty, False otherwise
        """
        return len(self.data[level]) == 0

    def pop(self, level):
        """
        Get the first element of the queue at given priority level
        :param level: priority level
        :return: first element at given priority level, None if this level is empty
        """
        if len(self.data[level]) > 0:
            self.size -= 1
            return self.data[level].popleft()
        return None

    def findClosestNonEmpty(self, level):
        """
        Finds the closest non empty level to the given priority level.
        :param level: priority level.
        :return: the closest non empty level
        """
        if not self.isLevelEmpty(level):
            return level
        lvlb = level - 1
        lvlu = level + 1
        while lvlb >= 0 or lvlu < self.levels:
            if lvlb >= 0:
                if not self.isLevelEmpty(lvlb):
                    return lvlb
                lvlb -= 1
            if lvlu < self.levels:
                if not self.isLevelEmpty(lvlu):
                    return lvlu
                lvlu += 1
        return None
