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


'''
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
'''
import logging
import copy
import HiPy.Processing.Attributes as Attributes
from HiPy.Util import VMath
from heapq import heappop, heappush
from enum import Enum


class HiPyException(Exception):
    """Base class for exceptions in this module."""
    pass


def initLogger():
    '''
    Creates a simple logger with a console handler for the module
    '''
    logger = logging.getLogger("HiPyLogger")
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.removeHandler(handler)
    logger.addHandler(handler)
    return logger

HiPyLogger = initLogger()

class TreeType(Enum):
    '''
    This enum lists the supported tree types: component tree and Partition hierarchy.
    '''
    ComponentTree = 1
    PartitionHierarchy = 2



class Image(list):
    '''
    A basic image represented as a linear array of values (can be anything).

    An image may be associated with a (valued) adjacency relation (sounds like a graph)

    An image may be associated with an embedding, mapping its linear space to another target space
    '''

    def __init__(self, length, initValue=0, adjacency=None, embedding=None):
        '''
        Creates a new image of "length" elements initialized to "initValue" (deep copied).

        It is strongly recommended to provide an adjacency and an embedding as many other
        functions will expect to find them.
        '''
        super().__init__()
        for _ in range(length):
            self.append(copy.deepcopy(initValue))
        self.embedding = embedding
        self.adjacency = adjacency

    def setAll(self, image):
        '''
        Copy (deep copy) values of image into current image
        '''
        for i in range(len(image)):
            self[i] = copy.deepcopy(image[i])

    def copy(self, copyData=False):
        '''
        Warning: this is not a true deep copy!
        Image data are not shared but adjacency and embedding are shared
        '''
        nIm = Image(len(self), 0, self.adjacency, self.embedding)

        if copyData:
            for i in range(len(self)):
                nIm[i] = self[i]
        return nIm

    def getPixel(self, i):
        '''
        Return the value of pixel i
        '''
        return self[i]

    def setPixel(self, i, val):
        '''
        Set the value of pixel i to val
        '''
        self[i] = val

    def setPixelWCS(self, val, *args):
        '''
        Set the value of the pixel of coordinate args in the image embedding to val
        '''
        self[self.embedding.getLinearCoordinate(*args)] = val

    def getPixelWCS(self, *args):
        '''
        Return the value of the pixel of coordinate args in the image embedding
        '''
        return self[self.embedding.getLinearCoordinate(*args)]

    def getNeighbours(self, i):
        '''
        Return the neighbours (direct succesors and predecors) of pixel i in the image adjacency
        '''
        return self.adjacency.getNeighbours(i)

    def equals(self, image, equalFunction=lambda x, y: x == y):
        '''
        Test if the current image is equal to the given image.

        Two images are equal if they have the same length and their pixel values are equal
        w.r.t. the given equalFunction (== by default)
        '''
        if len(self) != len(image):
            return False
        for i in range(len(image)):
            if not equalFunction(self[i], image[i]):
                return False
        return True

    def iterateOnPixels(self):
        '''
        Iterator on all pixels
        '''
        return range(len(self))

    def addAttribute(self, name, defaultValue=None, resetIfExist=False):
        '''
        Return the new attribute image if it did not exist or was reseted and None otherwise
        '''
        if not name in self.__dict__ or resetIfExist:
            image = Image(len(self), defaultValue, self.adjacency, self.embedding)
            self.__dict__[name] = image
            return image

        return None

    def deleteAttribute(self, name):
        '''
        Remove the attribute "name"
        '''
        if name in self.__dict__:
            del self.__dict__[name]

    def getAttribute(self, name):
        '''
        Return the attribute "name"
        '''
        if name in self.__dict__:
            return self.__dict__[name]
        else:
            return None


class AbstractEmbedding(object):
    '''
    Abstract embedding: ie a bijective mapping from linear indices
    (pixels in the abstract Image space) to any target space
    '''

    def getLinearCoordinate(self, *args):
        '''
        Convert coordinates from the target space to the abstract linear Image space
        '''
        raise NotImplementedError("Unsuported method" + " getLinearCoordinate")

    def fromLinearCoordinate(self, i):
        '''
        Convert coordinates from the abstract linear Image space to the target space
        '''
        raise NotImplementedError("Unsuported method" + " fromLinearCoordinate")

    def isInBoundsWCS(self, *args):
        '''
        Test if the given point is in the bounds of the embedding
        '''
        raise NotImplementedError("Unsuported method" + " isInBoundsWCS")

    def isInBoundsLinear(self, i):
        '''
        Test if the given point (in abstract linear image space) is in the bounds of the embedding
        '''
        return self.isInBoundsWCS(self.fromLinearCoordinate(i))

class Embedding2dGrid(AbstractEmbedding):
    '''
    This class represents an embedding from discrete linear coordinates to discrete 2d coordinates
    organized on a regular grid given by its width and height (typically a pixel grid).
    '''
    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.size = [width, height]

    def getLinearCoordinate(self, *args):
        return args[1]*self.width+args[0]

    def fromLinearCoordinate(self, i):
        return (i%self.width, i//self.width)

    def isInBoundsWCS(self, *args):
        return args[0] >= 0 and args[0] < self.width and args[1] >= 0 and args[1] < self.height


class AbstractAdjacency(object):
    '''
    Abstract adjacency relation
    '''

    def __init__(self, nbPoints):
        self.nbPoints = nbPoints

    def countEdges(self):
        '''
        Count the number of edges in the adjacency relation.

        Not that in a non-directed adjacency the edge (i,j) and
        the edge (j,i) are equivalent and count for one edge.
        '''
        count = 0
        for i in range(self.nbPoints):
            count += len(self.getSuccesors(i))
        return count

    def getSuccesors(self, i):
        '''
        return list of points j such that i->j
        '''
        raise NotImplementedError("Unsuported method" + " getSuccesors")

    def getPredecessors(self, i):
        '''
        return list of points j such that j->i
        '''
        raise NotImplementedError("Unsuported method" + " getPredecessors")

    def getNeighbours(self, i):
        '''
        return list of points j such that j->i or i->j (Succesors U Predecessors)
        '''
        raise NotImplementedError("Unsuported method" + " getNeighbours")

    def getEdges(self, i):
        '''
        return list of in  edges of the form [j,i,*] or [i,j,*] such that j->i or i->j
        (* can be any auxilary data, usually weights)
        '''
        return self.getInEdges(i)+self.getOutEdges(i)

    def getOutEdges(self, i):
        '''
        return list of outHead edges of the form [i,j,*] such that i->j
        (* can be any auxilary data, usually weights)
        '''
        raise NotImplementedError("Unsuported method" + " getOutEdges")

    def getInEdges(self, i):
        '''
        return list of in  edges of the form [j,i,*] such that j->i
        (* can be any auxilary data, usually weights)
        '''
        raise NotImplementedError("Unsuported method" + " getInEdges")

class AdjacencyNdRegular(AbstractAdjacency):
    '''
    Implicit representation of a shift-invariant adjacency relation
    in the n dimensional regular grid
    '''

    def __init__(self, embedding, neighbourList, weights=None):
        '''
        Defines the adjacency on the domain given by the embedding.

        The neighbouring relation is given by the neighbour list of the point 0.
        EG. a 4-adjacency in 2d is given by the neighbour list [ (0,-1), (-1,0), (1,0), (0,1)]

        Weights should be a list of same length as neighbour list giving the weight
        of the corresponding edge.
        If no weights are provided, all edges are assigned a wight of 1.
        '''
        super().__init__(VMath.mult(embedding.size))
        self.embedding = embedding
        self.neighbourList = neighbourList
        self.nbNeighbours = len(neighbourList)
        self.weights = weights if weights != None else [1]*self.nbNeighbours

    def countEdges(self):
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

    def getSuccesors(self, i):
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
        '''
        Creates an instance of AdjacencyNdRegular representing a classical 4 adjacency
        '''
        return AdjacencyNdRegular(Embedding2dGrid(size[0], size[1]), \
        [(0, -1), (-1, 0), (1, 0), (0, 1)])

    @staticmethod
    def getAdjacency2d8(size):
        '''
        Creates an instance of AdjacencyNdRegular representing a classical 8 adjacency
        '''
        return AdjacencyNdRegular(Embedding2dGrid(size[0], size[1]), \
        [(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)])



class AdjacencyTree(AbstractAdjacency):
    '''
    Represent the adjacency relation in a tree.

    The adjacency is oriented from the root toward the leaves
    '''
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

    def getSuccesors(self, i):
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
            distr = distr+dist(cur, par)
            if distr <= maxDist:
                neighbourList.append(par)
                cur = par
                par = tree[par]
            else:
                par = -1


        #search towards leaves
        stack = [[i, 0]]
        while len(stack) != 0:
            elem = stack.pop()
            cur = elem[0]
            diste = elem[1]
            for child in children[cur]:
                distd = diste+dist(cur, child)
                if distd <= maxDist:
                    neighbourList.append(child)
                    stack.append([child, distd])

        return neighbourList

class AbstractWeightedAdjacency(AbstractAdjacency, Image):
    '''
    An abstract edge weighted adjacency.

    A weighted adjacency is an image of weights with two attributes: source and target.
    Each pixel i of the image I represent an edge from source[i] to target[i] weighted by I[i].

     Known concrete implementations are
      - WeightedAdjacency
      - DirectedWeightedAdjacency
    '''
    def __init__(self, size):
        '''
         Create a new empty adjacency on a set of size elements
        '''
        AbstractAdjacency.__init__(self, size)
        Image.__init__(self, 0)
        self.addAttribute('source')
        self.addAttribute('target')

    def getSuccesors(self, i):
        '''
        return list of points j such that i->j
        '''
        raise NotImplementedError("Unsuported method" + " getSuccesors")





class WeightedAdjacency(AbstractWeightedAdjacency):
    '''
    Non directed adjacency represented with list of edges.

    Note:
      - methods getSuccesors, getPredecessors, and getNeighbours are equivalent
      - methods getEdges, getOutEdges, and getInEdges are equivalent
    '''

    def __init__(self, size):
        '''
         Create a new empty adjacency on a set of size elements
        '''
        AbstractWeightedAdjacency.__init__(self, size)
        self.edgeList = []
        for _ in range(size):
            self.edgeList.append([])

    def countEdges(self):
        return len(self)

    def createEdge(self, source, target, weight=1):
        '''
        Create a new edge in the graph, linking node source to node dest.

        As teh graph is undirected, by convention, the method will ensure that source < target

        Warning: does not verify is the edge already exists
        '''
        i = len(self)
        if source > target:
            source, target = target, source
        self.append(weight)
        self.source.append(source)
        self.target.append(target)

        self.edgeList[source].append(i)
        if source != target:
            self.edgeList[target].append(i)

    def getCopy(self):
        '''
        Returns a copy of the current adjacency.
        Garanties that edges indices are consistant between the copy and the current object
        (ie copy[i]==original[i] for all i)
        '''
        adj = WeightedAdjacency(self.nbPoints)
        for i in range(len(self)):
            adj.createEdge(self.source[i], self.target[i], self[i])
        return adj


    def getSuccesors(self, i):
        '''
        return list of points j such that i->j
        '''
        # ugly hack to symmetrise the adjacency on the fly
        return [self.source[e] + self.target[e] - i for e in self.edgeList[i]]

    def getPredecessors(self, i):
        '''
        return list of points j such that j->i
        '''
        return self.getSuccesors(i)

    def getNeighbours(self, i):
        '''
        return list of points j such that j->i or i->j (Succesors U Predecessors)
        '''
        return self.getSuccesors(i)

    def getEdges(self, i):
        '''
        return list of in  edges of the form [j,i,*] or [i,j,*] such that j->i or i->j
        (* can be any auxilary data, usually weights)
        '''
        return [[self.source[e], self.target[e], self[e]] for e in self.edgeList[i]]

    def getOutEdges(self, i):
        '''
        return list of outHead edges of the form [i,j,*] such that i->j
        (* can be any auxilary data, usually weights)
        '''
        return [[i, self.source[e]+self.target[e]-i, self[e]] for e in self.edgeList[i]]

    def getInEdges(self, i):
        '''
        return list of in  edges of the form [j,i,*] such that j->i
        (* can be any auxilary data, usually weights)
        '''
        return [[self.source[e]+self.target[e]-i, i, self[e]] for e in self.edgeList[i]]

    @staticmethod
    def createAdjacency(baseAdjacency, weightingFunction=None):
        '''
        Create a new adjacency equivalent to the given adjacency but
        with a different weighting function.

        Warning: the base adjacency is assumed to be symmetric!

        Typical use is to transform an implicit k-adjacency into an explicit weighted adjacency.
        '''
        adj = WeightedAdjacency(baseAdjacency.nbPoints)
        if weightingFunction != None:
            for i in range(adj.nbPoints):
                for j in baseAdjacency.getSuccesors(i):
                    if j > i:
                        adj.createEdge(i, j, weightingFunction(i, j))
        else:
            for i in range(adj.nbPoints):
                for j in baseAdjacency.getSuccesors(i):
                    if j > i:
                        adj.createEdge(i, j)

        return adj


    @staticmethod
    def createKHopAjacency(baseAdjacency, kValue, weightingFunction=lambda i, j, n: n):
        '''
        Construct an explicit non directed K-hop adjacency from a base (non directed adjacency).

        A weighting function f may be provided (given two nodes i, j , f(i,j,k) will be
        the weight of the edge linking i to j in k hops).

        '''
        nbPoints = baseAdjacency.nbPoints
        newAdj = WeightedAdjacency(nbPoints)
        flags = [None]*nbPoints

        for i in range(nbPoints):
            heap = []
            heappush(heap, (0, i))
            flags[i] = i
            while heap != []:
                (k, n) = heappop(heap)
                if n > i:
                    newAdj.createEdge(i, n, weightingFunction(i, n, k))

                if k < kValue:
                    k += 1
                    for j in baseAdjacency.getSuccesors(n):
                        if flags[j] != i:
                            heappush(heap, (k, j))
                            flags[j] = i

        return newAdj

    @staticmethod
    def createLineGraph(baseAdjacency, weightingFunction=lambda i, j: 1):
        '''
        Create the line graph of the given non directed base adjacency
        ie. a new adjacency relation were the points are the edges of the base adjacency
        and two points are adjacent is their correspind edges in the base adjacency
        share an extremity.

        The base adjacency must be  a WeightedAdjacency (Use WeightedAdjacency.createAdjacency
        in order to convert any non-directed adjacency to a WeightedAdjaceny).
        The medthod insure the consistency between the edges indice in the base adjacency and
        the points in the new adjacency.

        The given weighting function computes the weight of the edge linking the edge i to
        the edge j (of the base adjacency).
        '''
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
        '''
        Create a new Weighted Adjacency equals to the reflexive closure of
        the given base adjacency.

        The given weighting function computes the weight of the reflexive edges and
        takes a single parameter: the element i where the edge (i,i) is created.
        '''
        newAdj = WeightedAdjacency(baseAdjacency.nbPoints)
        for i in range(newAdj.nbPoints):
            for e in baseAdjacency.getOutEdges(i):
                if e[1] > i:
                    newAdj.createEdge(*e)

            newAdj.createEdge(i, i, weightingFunction(i))

        return newAdj

class DirectedWeightedAdjacency(AbstractWeightedAdjacency):
    '''
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
    '''

    def __init__(self, size):
        '''
         Create a new empty adjacency on a set of size elements
        '''
        AbstractWeightedAdjacency.__init__(self, size)
        self.outHead = [-1]*size
        self.outTail = [-1]*size
        self.prevEdge = []
        self.nextEdge = []

    @staticmethod
    def createAdjacency(baseAdjacency, weightingFunction=None):
        '''
        Create a new adjacency equivalent to the given adjacency but with
        a different weighting function.

        Typical use is to transform an implicit k-adjacency into an explicit weighted adjacency.
        '''
        adj = DirectedWeightedAdjacency(baseAdjacency.nbPoints)
        if weightingFunction != None:
            for i in range(adj.nbPoints):
                for j in baseAdjacency.getSuccesors(i):
                    adj.createEdge(i, j, weightingFunction(i, j))
        else:
            for i in range(adj.nbPoints):
                for j in baseAdjacency.getSuccesors(i):
                    adj.createEdge(i, j)

        return adj



    def getSuccesors(self, i):
        nodes = []
        e = self.outHead[i]
        while e != -1:
            # dest is the adjancent vertex
            nodes.append(self.target[e])
            # next edge in the edge list
            e = self.nextEdge[e]
        return nodes

    def getPredecessors(self, i):
        raise Exception("Unsuported method" + " getPredecessors")

    def getNeighbours(self, i):
        raise Exception("Unsuported method" + " getNeighbours")

    def getOutEdges(self, i):
        edges = []
        e = self.outHead[i]
        while e != -1:
            # dest is the adjancent vertex
            edges.append([self.source[e], self.target[e], self[e]])
            # next edge in the edge list
            e = self.nextEdge[e]
        return edges

    def getInEdges(self, i):
        raise Exception("Unsuported method" + " getInEdges")

    def createEdge(self, source, target, weight=1):
        '''
        Create a new edge in the graph, linking node source to node dest.
        '''
        i = len(self)
        self.append(weight)
        self.source.append(source)
        self.target.append(target)
        self.prevEdge.append(-1)
        self.nextEdge.append(-1)
        self.setEdge(source, target, i)

    def setEdge(self, source, _, i):
        '''
        Set the edge i of the graph, the edge is put at the head of the outHead linked list
        of the source vertex
        '''
        if self.outHead[source] == -1:
            self.outTail[source] = i
        else:
            head = self.outHead[source]
            self.prevEdge[head] = i
            self.nextEdge[i] = head
        self.outHead[source] = i

    def concatEdgeOut(self, vertex1, vertex2):
        '''
        Concatenate edge lists, outHead edges of vertex n2 are added to vertex n1
        '''
        if self.outTail[vertex1] != -1:
            self.nextEdge[self.outTail[vertex1]] = self.outHead[vertex2]
        if self.outHead[vertex2] != -1:
            self.prevEdge[self.outHead[vertex2]] = self.outTail[vertex1]
        if self.outTail[vertex2] != -1:
            self.outTail[vertex1] = self.outTail[vertex2]
        self.outHead[vertex2] = self.outTail[vertex2] = -1

    def removeEdgeOut(self, n, e):
        '''
        Remove the edge e of the list of outHead edges of vertex n
        '''
        if self.outHead[n] == e:
            self.outHead[n] = self.nextEdge[e]
        if self.outTail[n] == e:
            self.outTail[n] = self.prevEdge[e]
        if self.nextEdge[e] != -1:
            self.prevEdge[self.nextEdge[e]] = self.prevEdge[e]
        if self.prevEdge[e] != -1:
            self.nextEdge[self.prevEdge[e]] = self.nextEdge[e]
        self.nextEdge[e] = self.prevEdge[e] = -1


    def getCopy(self):
        '''
        Creates a copy of the current adjacency
        '''
        adj = DirectedWeightedAdjacency(self.nbPoints)
        for i in range(len(self)):
            adj.createEdge(self.source[i], self.target[i], self[i])
        return adj

    def getSymmetricMax(self):
        '''
        Symmetrize the graph with the max strategy: whenever an edge (p,q) is found
        the edge (q,p) is added to the graph
        '''
        graph = DirectedWeightedAdjacency(self.nbPoints)
        for i in range(self.nbPoints):
            for j in self.getSuccesors(i):
                if not j in graph.getSuccesors(i):
                    graph.createEdge(i, j)
                if not i in self.getSuccesors(j):
                    graph.createEdge(j, i)
        return graph

    def getSymmetricMin(self):
        '''
        Symmetrize the graph with the min strategy: an edge (p,q) is preserved only if
        the edge (q,p) is also in the graph
        '''
        graph = DirectedWeightedAdjacency(self.nbPoints)
        for i in range(self.nbPoints):
            for j in self.getSuccesors(i):
                if j > i and i in self.getSuccesors(j):
                    graph.createEdge(i, j)
                    graph.createEdge(j, i)
        return graph


    def getTranspose(self):
        '''
        Transpose the graph: each edge (p,q) is transformed into the edge (q,p)
        '''
        graph = DirectedWeightedAdjacency(self.nbPoints)
        for i in range(self.nbPoints):
            for j in self.getSuccesors(i):
                graph.createEdge(j, i)
        return graph




class Tree(Image):
    '''
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
    '''
    def __init__(self, treeType, parent, levels, image=None):
        super(Tree, self).__init__(len(parent), None)
        HiPyLogger.debug("constructor: Tree")
        self.treeType = treeType
        self.adjacency = AdjacencyTree(self)
        self.setAll(parent)
        self.leavesAdjacency = image.adjacency if image != None else None
        self.leavesEmbedding = image.embedding if image != None else None
        self.nbPixels = len(image) if image != None else Tree._countLeaves(parent)
        self.addAttribute("level")
        self.level.setAll(levels)
        self.addAttribute('reconstructedValue')
        Attributes.addAttributeChildren(self)


    def getParent(self, i):
        '''
        Return the parent of the node i, -1 if i is a root
        '''
        return self[i]



    def filterDirect(self, filteringCriterion):
        '''
        Apply a given filtering criterion f on the tree.
        A filtering criterion is a function that associates the value true or false to a tree
        and a node.

        For each node n the attribute deleted is set to the value of filteringCriterion(self,n)
        '''
        self.addAttribute("deleted", False)
        for i in self.iteratorFromPixelsToRoot():
            self.deleted[i] = filteringCriterion(self, i)

    def extractBranch(self, node, attributeName=None):
        '''
        Extract values along a branch of the tree, from the given node toward the root.

        If attributeName==False, the result is the list of nodes from the given node to the root.

        Otherwise the result is the list of node attribute values from the given node to the root
        for the attribute "attributeName"
        '''
        attr = None
        if  attributeName != None:
            attr = self.getAttribute(attributeName)

        res = []
        res.append(node if attr == None else attr[node])
        node = self[node]
        while node != -1:
            res.append(node if attr == None else attr[node])
            node = self[node]

        return res


    def reconstructImage(self, attributeName="level", criterion=None):
        '''
        Reconstruct an image using the value of the attribute "attributeName".
        The final pixel value is given by the value of its closest selected (see below) parent.

        If a criterion is provided (a function that associates True or False to each node),
        then a node x is selected if criterion(x)==False

        If no criterion is provided and the attribute "deleted" is defined, a node x is selected
        if deleted[x]==False

        If no criterion is provided and the attribute "deleted" is not defined then every nodes
        are selected.
        '''
        if criterion == None:
            if "deleted" in self.__dict__:
                criterion = (lambda x: self.deleted[x])
            else:
                criterion = (lambda _: False)
        root = len(self)-1
        for i in self.iteratorFromRootToPixels():
            if i >= self.nbPixels and (i == root or  not criterion(i)):
                self.reconstructedValue[i] = self.__dict__[attributeName][i]
            else:
                self.reconstructedValue[i] = self.reconstructedValue[self[i]]

        image = Image(self.nbPixels, 0, self.leavesAdjacency, self.leavesEmbedding)
        for i in range(len(image)):
            image[i] = self.reconstructedValue[i]

        return image

    def iteratorFromLeavesToRoot(self, includeLeaves=True, includeRoot=True):
        '''
        Provides an iterator on the nodes of the tree going from the leaves to the root.

        '''
        return TreeIterator(self, True, False, includeLeaves, includeRoot)

    def iteratorFromRootToLeaves(self, includeLeaves=True, includeRoot=True):
        '''
        Provides an iterator on the nodes of the tree going from the root to the leaves.
        '''
        return TreeIterator(self, True, True, includeLeaves, includeRoot)

    def iteratorFromPixelsToRoot(self, includePixels=True, includeRoot=True):
        '''
        Provides an iterator on the nodes of the tree going from the pixels to the root.

        '''
        return TreeIterator(self, False, False, includePixels, includeRoot)

    def iteratorFromRootToPixels(self, includePixels=True, includeRoot=True):
        '''
        Provides an iterator on the nodes of the tree going from the root to the pixels.
        '''
        return TreeIterator(self, False, True, includePixels, includeRoot)

    def iteratorOnPixels(self):
        '''
        Creates an iterator on the pixels of the tree
        '''
        return range(self.nbPixels)

    def iteratorOnLeaves(self):
        '''
        Creates an iterator on the leaves of the tree.
        '''
        if self.treeType == TreeType.ComponentTree:
            Attributes.addAttributeIsLeaf(self)
            isLeaf = self.isLeaf
            i = self.nbPixels
            size = len(self)
            while i < size:
                while i < size and not isLeaf[i]:
                    i += 1
                if i < size:
                    yield i
                i += 1
        else:
            return range(self.nbPixels)

    def lca(self, i, j):
        '''
        Lowest common ancestor between nodes i and j
        Warning: Assume that the attribute depth is present !
        '''
        depth = self.depth
        while i != j:
            if depth[i] > depth[j]:
                i = self[i]
            else:
                j = self[j]

        return i

    def nbNodes(self):
        '''
        Count the number of logical nodes.

        For component trees, pixels that are not part of the logical structure must be substracted
        leading to len(tree)-tree.nbPixels
        For partition hierarchies or others the result is equal to len(tree)
        '''
        if self.treeType == TreeType.ComponentTree:
            return len(self)-self.nbPixels
        else:
            return len(self)



    @staticmethod
    def _countLeaves(parent):
        '''
        Count Leaves in a parent relation
        '''
        leaves = [True]*len(parent)
        for i in range(len(parent)):
            par = parent[i]
            if par != -1:
                leaves[par] = False

        count = 0
        for i in range(len(parent)):
            if leaves[i]:
                count = count+1

        return count

    @staticmethod
    def testTreeIsomorphism(tree1, tree2):
        '''
        Test if tree1 and tree2 are isomorph assuming that leaves are ordered
        '''
        if len(tree1) != len(tree2) or tree1.nbPixels != tree2.nbPixels:
            return False

        #both tree have same size so we need to find an injection m from the nodes of t1 to
        # the nodes of t2 such that for any node n of t1 m(parent(n))=parent(m(n))
        mapT1T2 = [-1]*len(tree1)

        for i in range(len(mapT1T2)-1): #root is root !
            #pixel mapping is constant
            if i < tree1.nbPixels:
                mapT1T2[i] = i
            #parent(n)
            pT1 = tree1[i]

            #parent(m(n))
            pT2 = tree2[mapT1T2[i]]
            if mapT1T2[pT1] == -1:
                mapT1T2[pT1] = pT2
            elif mapT1T2[pT1] != pT2:
                return False

        return True


class TreeIterator(object):
    '''
    The TreeIterator class is intended to perform pre or post order traversal of trees.

    It manages several common cases as exluding or including leaves or the root.

    It also manages the logical difference between Component trees
    (where leaves represent pixels of the original image but not nodes of the logical tree)
    and partition hierarchies (where leaves are at the same time pixels of the original image
    and nodes of the logical tree).
    '''
    def __init__(self, tree, logical=True, reverseOrder=False, includeLeaves=True, \
    includeRoot=True):
        self.tree = tree

        specialLogic = tree.treeType == TreeType.ComponentTree and logical
        if specialLogic and not includeLeaves:
            Attributes.addAttributeIsLeaf(tree)
            self.nextMethod = self.nextLogical
        else:
            self.nextMethod = self.nextStructural

        if specialLogic or not includeLeaves:
            bmin = tree.nbPixels
        else:
            bmin = 0

        bmax = len(tree)
        if not includeRoot:
            bmax = bmax-1

        if not reverseOrder:
            self.curVal = bmin
            self.limit = bmax
            self.step = 1
        else:
            self.curVal = bmax-1
            self.limit = bmin-1
            self.step = -1




    def __iter__(self):
        return self

    def nextLogical(self):
        '''
        "next" iterator method when perfoming a logical browsing (only logical nodes are traversed)
        '''
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
        '''
        "next" iterator method when perfoming a structural browsing (every node is traversed)
        '''
        if self.curVal != self.limit:
            i = self.curVal
            self.curVal += self.step
            return i
        else:
            raise StopIteration()

    def __next__(self):
        return self.nextMethod()

