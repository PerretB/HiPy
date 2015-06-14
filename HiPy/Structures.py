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
Created on 3 juin 2015

@author: perretb
'''
import copy

def enum(**enums):
    return type('Enum', (), enums)

class Image(list):
    '''
    A basic image represented as a linear array of values (can be anything).
    
    An image may be associated with a (valued) adjacency relation (sounds like a graph)
    
    An image may be associated with an embedding, mapping its linear space to another target space
    '''
    
    def __init__(self, length, initValue=0,adjacency=None, embedding=None):
        '''
        Creates a new image of "length" elements initialized to "initValue" (deep copied).
        
        It is strongly recommended to provide an adjacency and an embedding as many other functions will expect to find them.
        '''
        for _ in range(length):
            self.append(copy.deepcopy(initValue))
        self.embedding=embedding
        self.adjacency=adjacency   
    
    def setAll(self,image):
        '''
        Copy values of image into current image
        '''
        for i in range(len(image)):
            self[i]=image[i]
            
    def copy(self,copyData=False):
        '''
        Warning: this is not a true deep copy! 
        Image data are not shared but adjacency and embedding are shared
        '''
        nIm = Image(len(self),0,self.adjacency,self.embedding)
       
        if copyData:
            for i in range(len(self)):
                nIm[i]=self[i]
        return nIm

    def getPixel(self,i):
        '''
        Return the value of pixel i
        '''
        return self[i]

    def setPixel(self,i,val):
        '''
        Set the value of pixel i to val
        '''
        self[i]=val
        
    def setPixelWCS(self,val,*args):
        '''
        Set the value of the pixel of coordinate args in the image embedding to val
        '''
        self[self.embedding.getLinearCoordinate(*args)]=val
        
    def getPixelWCS(self,*args):
        '''
        Return the value of the pixel of coordinate args in the image embedding
        '''
        return self[self.embedding.getLinearCoordinate(*args)]

    def getNeighbours(self, i):
        '''
        Return the neighbours (direct succesors and predecors) of pixel i in the image adjacency
        '''
        return self.adjacency.getNeighbours(i)

    def equals(self, image, equalFunction=lambda x,y:x==y):
        '''
        Test if the current image is equal to the given image.
        
        Two images are equal if they have the same length and their pixel values are equal w.r.t. the given equalFunction (== by default)
        '''
        if len(self)!=len(image):
            return False
        for i in range(len(image)):
            if not equalFunction(self[i],image[i]):
                return False
        return True

class Embedding(object):
    '''
    Abstract embedding: ie a bijective mapping from linear indices 
    (pixels in the abstract Image space) to any target space
    '''
    
    def getLinearCoordinate(self, *args):
        '''
        Convert coordinates from the target space to the abstract linear Image space
        '''
        raise Exception("Unsuported method" + " getLinearCoordinate")
    
    def fromLinearCoordinate(self, i):
        '''
        Convert coordinates from the abstract linear Image space to the target sapce
        '''
        raise Exception("Unsuported method" + " fromLinearCoordinate")

class Embedding2dGrid(Embedding):
    
    def __init__(self,width,height):
        self.width=width
        self.height=height
        self.size=[width,height]
    
    def getLinearCoordinate(self, *args):
        return args[1]*self.width+args[0]
    
    def fromLinearCoordinate(self, i):
        return (i%self.width,i//self.width)

class Adjacency(object):
    '''
    Abstract directed adjacency relation
    '''
    
    def getSuccesors(self, i):
        '''
        return list of points j such that i->j
        '''
        raise Exception("Unsuported method" + " getSuccesors")
    
    def getPredecessors(self, i):
        '''
        return list of points j such that j->i
        '''
        raise Exception("Unsuported method" + " getPredecessors")
    
    def getNeighbours(self, i):
        ''' 
        return list of points j such that j->i or i->j (Succesors U Predecessors)
        '''
        raise Exception("Unsuported method" + " getNeighbours")
    
    def getEdges(self,i):
        ''' 
        return list of in  edges of the form [j,i,*] or [i,j,*] such that j->i or i->j(* can be any auxilary data, usually weights)
        '''
        return self.getInEdges(i)+self.getOutEdges(i)
    
    def getOutEdges(self,i):
        ''' 
        return list of out edges of the form [i,j,*] such that i->j (* can be any auxilary data, usually weights)
        '''
        raise Exception("Unsuported method" + " getOutEdges")
        
    def getInEdges(self,i):
        ''' 
        return list of in  edges of the form [j,i,*] such that j->i (* can be any auxilary data, usually weights)
        '''
        raise Exception("Unsuported method" + " getInEdges")
    

class Adjacency2d4(Adjacency):
    '''
    Represent a symmetric 4-neighbourhood in a regular 2D rectangular grid
    '''
    
    def __init__(self,size):
        self.size=size
    
    def __getCoordsLin(self,x,y):
        return y*self.size[0]+x
    
    def __getCoords2D(self,i):
        return (i%self.size[0],i//self.size[0])
    
    def getNeighbours(self, i):
        (x,y)=self.__getCoords2D(i)
        nl=[]
        if x-1>=0:
            nl.append(self.__getCoordsLin(x-1,y))
        if y-1>=0:
            nl.append(self.__getCoordsLin(x,y-1))
        if x+1<self.size[0]:
            nl.append(self.__getCoordsLin(x+1,y))
        if y+1<self.size[1]:
            nl.append(self.__getCoordsLin(x,y+1))
        return nl
    
    def getPredecessors(self, i):
        return self.getNeighbours(i)
    
    def getSuccesors(self, i):
        return self.getNeighbours(i)

    def getOutEdges(self,i):
        (x,y)=self.getCoords2D(i)
        nl=[]
        if x-1>=0:
            nl.append([i,self.__getCoordsLin(x-1,y)])
        if y-1>=0:
            nl.append([i,self.__getCoordsLin(x,y-1)])
        if x+1<self.size[0]:
            nl.append([i,self.__getCoordsLin(x+1,y)])
        if y+1<self.size[1]:
            nl.append([i,self.__getCoordsLin(x,y+1)])
        return nl
    
    def getInEdges(self,i):
        (x,y)=self.__getCoords2D(i)
        nl=[]
        if x-1>=0:
            nl.append([self.__getCoordsLin(x-1,y),i])
        if y-1>=0:
            nl.append([self.__getCoordsLin(x,y-1),i])
        if x+1<self.size[0]:
            nl.append([self.__getCoordsLin(x+1,y),i])
        if y+1<self.size[1]:
            nl.append([self.__getCoordsLin(x,y+1),i])
        return nl
    

class AdjacencyTree(Adjacency):
    '''
    Represent the adjacency relation in a tree.
    
    The adjacency is oriented from the root toward the leaves
    '''

    def __init__(self,tree):
        self.tree=tree
        
    def getNeighbours(self, i):
        nl=[]
        if self.tree[i]!=-1:
            nl.append(self.tree[i])
        nl.extend(self.tree.children[i])
        return nl

    def getPredecessors(self, i):
        if self.tree[i]!=-1:
            return [self.tree[i]]
        return [];
            
    def getSuccesors(self, i):
        return self.tree.children[i]

    def getOutEdges(self,i):
        nl=[]
        for c in self.tree.children[i]:
            nl.append([i,c])
        return nl
    
    def getInEdges(self,i):
        if self.tree[i]!=-1:
            return [self.tree[i],i]
        return [];

    # @todo: refactor
    def getNeigbourPointsExtended(self,i,maxDist,dist=(lambda x,y:1)):
        nl=[i];
        tree=self.tree
        children=self.tree.children
        
        cur=i
        par=tree[cur]
        dr=0    
        while par!=-1:
            dr=dr+dist(cur,par)
            if(dr<=maxDist):
                nl.append(par)
                cur=par
                par=tree[par]
            else:
                par=-1
           
                    
        #search towards leaves
        stack=[[i,0]]
        while len(stack)!=0:
            e=stack.pop()
            cur=e[0]
            d=e[1]
            for c in children[cur]:
                dd = d+dist(cur,c)
                if(dd<=maxDist):
                    nl.append(c)
                    stack.append([c,dd])
                
        return nl


class AdjacencyEdgeWeightedGraph(Adjacency):
    '''
    Generic adjacency represented using link lists of out edges
    '''
    
    def __init__(self,size):
        '''
         Create a new adjacency on a set of size elements
        '''
        self.nbPoints=size
        self.out=[-1]*size
        self.outTail=[-1]*size
        # an edge is an array of two elements [source,destination] so this is kind of a hybrid graph representation
        self.edges=[]
        self.prevEdge=[]
        self.nextEdge=[]

 
    
    def getSuccesors(self, i):
        nodes = []
        e = self.out[i]
        while(e != -1):
            # dest is the adjancent vertex
            nodes.append(self.edges[e][1])
            # next edge in the edge list
            e = self.nextEdge[e]
        return nodes
    
    def getPredecessors(self, i):
        raise Exception("Unsuported method" + " getPredecessors")
    
    def getNeighbours(self, i):
        raise Exception("Unsuported method" + " getNeighbours")
    
    def getOutEdges(self,i):
        edges = []
        e = self.out[i]
        while(e != -1):
            # dest is the adjancent vertex
            edges.append(self.edges[e])
            # next edge in the edge list
            e = self.nextEdge[e]
        return edges
        
    def getInEdges(self,i):
        raise Exception("Unsuported method" + " getInEdges")
    
    def createEdge(self, source, dest, weight=None):
        '''
        Create a new edge in the graph, linking node source to node dest.
        '''
        i=len(self.prevEdge)
        if weight!=None:
            self.edges.append([source,dest,weight])
        else:
            self.edges.append([source,dest])
        self.prevEdge.append(-1)
        self.nextEdge.append(-1)
        self.setEdge(source, dest, i)

    def setEdge(self, source, dest, i):
        '''
        Set the edge i of the graph, the edge is put at the head of the out linked list of the source vertex
        '''
        n1 = source
        if self.out[n1] == -1:
            self.outTail[n1] = i
        else:
            e1 = self.out[n1]
            self.prevEdge[e1] = i
            self.nextEdge[i] = e1
        self.out[n1] = i

    def concatEdgeOut(self, n1, n2):
        '''
        Concatenate edge lists, out edges of vertex n2 are added to vertex n1
        '''
        if self.outTail[n1] != -1:
            self.nextEdge[self.outTail[n1]] = self.out[n2]
        if self.out[n2] != -1:    
            self.prevEdge[self.out[n2]] = self.outTail[n1]
        if self.outTail[n2] != -1:
            self.outTail[n1] = self.outTail[n2]
        self.out[n2] = self.outTail[n2] = -1

    def removeEdgeOut(self, n, e):
        '''
        Remove the edge e of the list of out edges of vertex n
        '''
        if self.out[n] == e:
            self.out[n] = self.nextEdge[e]
        if self.outTail[n] == e:
            self.outTail[n] = self.prevEdge[e]
        if self.nextEdge[e] != -1:
            self.prevEdge[self.nextEdge[e]] = self.prevEdge[e]
        if self.prevEdge[e] != -1:
            self.nextEdge[self.prevEdge[e]] = self.nextEdge[e]
        self.nextEdge[e] = self.prevEdge[e] = -1
   

    def getSymmetricMax(self):
        '''
        Symmetrize the graph with the max strategy: whenever an edge (p,q) is found the edge (q,p) is added to the graph
        '''
        g = AdjacencyEdgeWeightedGraph(self.nbPoints)
        for i in range(self.nbPoints):
            for j in self.getSuccesors(i):
                if not j in g.getSuccesors(i):
                    g.createEdge(i,j)
                if not i in self.getSuccesors(j):
                    g.createEdge(j,i)
        return g

    def getSymmetricMin(self):
        '''
        Symmetrize the graph with the min strategy: an edge (p,q) is preserved only if the edge (q,p) is also in the graph
        '''
        g = AdjacencyEdgeWeightedGraph(self.nbPoints)
        for i in range(self.nbPoints):
            for j in self.getSuccesors(i):
                if j>i and i in self.getSuccesors(j):
                        g.createEdge(i,j)
                        g.createEdge(j,i)
        return g

 
    def getTranspose(self):
        '''
        Transpose the graph: each edge (p,q) is transformed into the edge (q,p)
        '''
        g = AdjacencyEdgeWeightedGraph(self.nbPoints)
        for i in range(self.nbPoints):
            for j in self.getSuccesors(i):
                g.createEdge(j,i)
        return g




class Tree(Image):
    '''
    A tree is an image whose pixel values encode the parent relation
    '''
    def __init__(self, parent,levels,image):
        super(Tree,self).__init__(len(parent),None)
        self.setAll(parent)
        self.leavesAdjacency = image.adjacency
        self.leavesEmbedding = image.embedding
        self.nbLeaves=len(image)
        self.level=Image(len(levels))
        self.level.setAll(levels)
        self.level.adjacency=AdjacencyTree(self)
        self.addAttribute('reconstructedValue')
            
    def getParent(self,i):
        return self[i]
        
    def addAttribute(self,name,defaultValue=None):
        '''
        Return the new attribute image if it did not exist and None otherwise
        '''
        if(not name in self.__dict__):
            im=Image(len(self),defaultValue,AdjacencyTree(self))
            self.__dict__[name]=im
            return im
        return None
            
    def getAttribute(self,name):
        return self.__dict__[name]

    # apply a given selection criterion f on the tree
    # the selection criterion is a function that associates the value true or false to a node
    def filterDirect(self, f):
        attr=self.addAttribute("deleted",False)
        if attr==None:
            for i in self.iterateFromLeavesToRoot():
                self.deleted[i] = False
    
        for i in self.iterateFromLeavesToRoot():
            self.deleted[i] = f(self,i)

    def extractBranch(self,node,attributeName=None):
        '''
        Extract values along a branch of the tree, from the given node toward the root.
        
        If attributeName==False, the result is the list of nodes from the given node to the root.
        
        Otherwise the result is the list of node attribute values from the given node to the root for the attribute "attributeName"
        '''
        attr=None
        if  attributeName!=None:
            attr=self.getAttribute(attributeName)
        
        res=[]
        res.append(node if attr==None else attr[node])   
        node=self[node]
        while node !=-1:
            res.append(node if attr==None else attr[node])   
            node=self[node]
            
        return res


    def reconstructImage(self,attributeName="level",criterion=None):
        '''
        Reconstruct an image using the value of the attribute "attributeName".
        The final pixel value is given by the value of its closest selected (see below) parent.
        
        If a criterion is provided (a function that associates True or False to each node), 
        then a node x is selected if criterion(x)==True
        
        If no criterion is provided and the attribute "deleted" is defined, a node x is selected if deleted[x]==False 
        
        If no criterion is provided and the attribute "deleted" is not defined then every nodes are selected.
        '''
        if criterion==None:
            if "deleted" in self.__dict__:
                criterion = (lambda x:not self.deleted[x])
            else:
                criterion = (lambda _:True)
        root= len(self)-1
        for i in self.iterateFromRootToLeaves(True):
            if i>=self.nbLeaves and (i==root or  criterion(i)) :
                self.reconstructedValue[i]=self.__dict__[attributeName][i]
            else:
                self.reconstructedValue[i]=self.reconstructedValue[self[i]]
                    
        im = Image(self.nbLeaves,0,self.leavesAdjacency,self.leavesEmbedding)
        for i in range(len(im)):
            im[i]=self.reconstructedValue[i]
        return im
    
    def iterateFromLeavesToRoot(self,includePixels=True):
        '''
        Provides an iterator on the nodes of the tree going from the leaves to the root.
        
        Leaves (pixels) can be included or not
        '''
        if includePixels:
            return range(len(self))
        else:
            return range(self.nbLeaves,len(self),1)
        
    def iterateFromRootToLeaves(self,includePixels=True):
        '''
        Provides an iterator on the nodes of the tree going from the root to the leaves.
        
        Leaves (pixels) can be included or not
        '''
        if includePixels:
            return range(len(self)-1,-1,-1)
        else:
            return range(len(self)-1,self.nbLeaves-1,-1)