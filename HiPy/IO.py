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


from HiPy.Structures import Image, Embedding2dGrid, AdjacencyEdgeWeightedGraph
import os

PILAvailable=False
# for image I/O
try:
    from PIL import Image as PILimage
    PILAvailable=True
except ImportError:
    print("Error: PIL library not available !")
    
    
def ensureDirectoryExists(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)    

def readImage(filename, grayScale=True): 
    if not PILAvailable:
        raise Exception("PIL needed to load image")
    image= PILimage.open(filename)
    size = image.size
    coord = lambda x: (x % size[0], x // size[0])
    pixdata = image.load()
    im = Image(size[0]*size[1])
    if grayScale and isinstance(pixdata[(0,0)], (list, tuple)):
        if len(pixdata[(0,0)])>3:
            print("Warning: Image.Read, image has more than 3 channels (alpha component?), one or more channels will be ignored during the grayscale conversion.")
        for i in range(size[0]*size[1]):
            im[i]=sum(pixdata[coord(i)][0:3])//3
    else:
        for i in range(size[0]*size[1]):
            im[i]=(pixdata[coord(i)])
    im.embedding = Embedding2dGrid(size[0],size[1])
    return im

def saveImage(image,filename):
    if not PILAvailable:
        raise Exception("PIL needed to save image")
    ensureDirectoryExists(filename)
    width=image.embedding.width;
    size=[width,image.embedding.height]
    if isinstance(image[0], (list, tuple)):
        im = PILimage.new("RGB", size, 0)
    else:
        im = PILimage.new("L", size, 0)
    pix = im.load()
    
    for i in range(len(image)):
        pix[i%width,i//width]=image[i]
    im.save(filename, "PNG")
    
    


def readGraph(filename):
    '''
    Read a graph (AdjacencyEdgeWeightedGraph) from ascii file.
    format (weights are optionnale, vertices are numbered from 0 to numberOfVertices-1):
    numberOfVertices numberOfEdges
    sourceVertexOfEdge_1 destinationVertexOfEdge_1 weightOfEdge_1
    sourceVertexOfEdge_2 destinationVertexOfEdge_2 weightOfEdge_2
    ...
    sourceVertexOfEdge_NumberOfEdges destinationVertexOfEdge_NumberOfEdges weightOfEdge_NumberOfEdges
    '''
    infile = open(filename, "r")
    line = infile.readline()
    while(line[0] == '#'):
        line = infile.readline()
    nbPoints = int(line.split()[0])
    lines = infile.readlines()
    infile.close()
    graph=AdjacencyEdgeWeightedGraph(nbPoints)
    
    
    for line in lines:
        tokens = line.split()
        if len(tokens)==3:
            graph.createEdge(int(tokens[0]), int(tokens[1]), int(tokens[2]))
        else:
            graph.createEdge(int(tokens[0]), int(tokens[1]))
    
    return graph

def saveGraph(graph, filename):
    '''
    Save a graph (AdjacencyEdgeWeightedGraph) in a file as plain text
    '''
    out = open(filename, "w")
    out.write(str(graph["nbPoints"]) + " " + str(len(graph["edges"])) + "\n")
    for e in graph["edges"]:
        line=" ".join([str(v) for v in e])
        out.write(line)
        out.write("\n")
    out.close()