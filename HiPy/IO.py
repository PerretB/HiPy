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
This module is dedicated to input/output function.

It allows to read/write images and graphs.
Created on 3 juin 2015

@author: Benjamin Perret
"""

import HiPy.Structures
import os
import pickle
from HiPy.Structures import Embedding2dGrid

PILAvailable = False
# for image I/O


try:
    from PIL import Image as PILimage
    PILAvailable = True
except ImportError:
    PILimage = None
    HiPy.Structures.HiPyLogger.warning("Warning: PIL library not available !")

# for drawing graphs in a pdf file
try:
    import pydot

    pydotAvailable = True
except ImportError:
    pydotAvailable = False
    HiPy.Structures.HiPyLogger.warning("Warning: pydot library not available !")


def ensureDirectoryExists(directory):
    """
    Check if a directory exists and creates it if not.
    :param directory: the directory path to check
    :return: void
    """
    d = os.path.dirname(directory)
    if d != "" and not os.path.exists(d):
        os.makedirs(d)


def readImage(filename, grayScale=False):
    """
    Reads an image from a file an returns the corresponding HiPy.Structures.Image object.

    FIXME: handle alpha channel.

    :param filename: path to the image file to read
    :param grayScale: if True the result will be forced to be gray scale
    :return: an instance of HiPy.Structures.Image
    """
    if not PILAvailable:
        raise Exception("PIL needed to load image")
    image = PILimage.open(filename)
    size = image.size
    coord = lambda x: (x % size[0], x // size[0])
    pixelData = image.load()
    im = HiPy.Structures.Image(size[0] * size[1])
    bands = 1
    if isinstance(pixelData[(0, 0)], (list, tuple)):
        bands = len(pixelData[(0, 0)])
    if grayScale and bands > 1:
        if bands > 3:
            print(
                "Warning: Image.Read, image has more than 3 channels (alpha component?), one or more channels will "
                "be ignored during the grayscale conversion.")
        for i in range(size[0] * size[1]):
            im[i] = sum(pixelData[coord(i)][0:3]) // 3
    else:
        if bands == 1:
            for i in range(size[0] * size[1]):
                im[i] = pixelData[coord(i)]
        else:
            for i in range(size[0] * size[1]):
                im[i] = list(pixelData[coord(i)])
    im.embedding = HiPy.Structures.Embedding2dGrid(size[0], size[1])
    return im


def saveImage(image, filename):
    """
    Save the given image at the given location.

    Always save in PNG format for now...

    TODO: add support to other format
    :param image: the image to save
    :param filename: path where to save the image
    :return: void
    """
    if not PILAvailable:
        raise Exception("PIL needed to save image")
    ensureDirectoryExists(filename)
    width = image.embedding.width
    size = [width, image.embedding.height]
    if isinstance(image[0], tuple):
        im = PILimage.new("RGB", size, 0)
        pix = im.load()

        for i in range(len(image)):
            pix[i % width, i // width] = image[i]
    elif isinstance(image[0], list):
        im = PILimage.new("RGB", size, 0)
        pix = im.load()

        for i in range(len(image)):
            pix[i % width, i // width] = tuple(image[i])
    else:
        im = PILimage.new("L", size, 0)
        pix = im.load()

        for i in range(len(image)):
            pix[i % width, i // width] = image[i]

    im.save(filename, "PNG")


def readEdgeGraph(filename, directed=True):
    """
    Read a graph (DirectedWeightedAdjacency) from ascii file.

    format (weights are optional, vertices are numbered from 0 to numberOfVertices-1):
    numberOfVertices numberOfEdges
    sourceVertexOfEdge_1 destinationVertexOfEdge_1 weightOfEdge_1
    sourceVertexOfEdge_2 destinationVertexOfEdge_2 weightOfEdge_2
    ...
    sourceVertexOfEdge_NumberOfEdges destinationVertexOfEdge_NumberOfEdges weightOfEdge_NumberOfEdges

    :param filename: path to the file to read

    """
    infile = open(filename, "r")
    line = infile.readline()
    while line[0] == '#':
        line = infile.readline()
    nbPoints = int(line.split()[0])
    lines = infile.readlines()
    infile.close()
    if directed:
        graph = HiPy.Structures.DirectedWeightedAdjacency(nbPoints)
    else:
        graph = HiPy.Structures.WeightedAdjacency(nbPoints)

    for line in lines:
        tokens = line.split()
        origin = int(tokens[0])
        target = int(tokens[1])
        if len(tokens) == 3:
            graph.createEdge(origin, target, int(tokens[2]))
        else:
            graph.createEdge(origin, target)

    return graph


def saveGraph(abstractAdjacency, image, filename, directed=True):
    """
    Save a graph (DirectedWeightedAdjacency) in a file as plain text (hugly) pink format

    :param abstractAdjacency: the adjacency to save
    :param filename: path to the saved file
    """
    out = open(filename, "w")
    embedding = image.embedding
    out.write("#rs " + str(embedding.width) + " cs " + str(embedding.height) + "\n")
    #rs 481 cs 321
    out.write(str(abstractAdjacency.nbPoints) + " " + str(abstractAdjacency.countEdges()) + "\n")
    out.write("val sommets\n")
    for i in image.iterateOnPixels():
        out.write(str(i) + " 1\n")
        #out.write(str(i) + " " + str(image[i]) + "\n")
    out.write("arcs values\n")
    for v in range(abstractAdjacency.nbPoints):
        for e in abstractAdjacency.getOutEdges(v):
            if directed or e[1] > e[0]:
                line = " ".join([str(v) for v in e])
                out.write(line)
                out.write("\n")
    out.close()

def readGraph(filename, directed=False):
    """
    Read a graph (WeightedAdjacency) from a file as plain text (hugly) pink format



    :param filename: path to the file to read

    """
    infile = open(filename, "r")
    line = infile.readline()
    embedding = None
    while line[0] == '#':
        if line.startswith("#rs"):
            tokens = line.split()
            embedding = Embedding2dGrid(int(tokens[1]), int(tokens[3]))
        line = infile.readline()
    nbPoints = int(line.split()[0])
    nbEdges = int(line.split()[1])
    lines = infile.readlines()
    infile.close()
    if directed:
        graph = HiPy.Structures.DirectedWeightedAdjacency(nbPoints)
    else:
        graph = HiPy.Structures.WeightedAdjacency(nbPoints)

    for i in range(2+nbPoints,len(lines)):
        line = lines[i]
        tokens = line.split()
        origin = int(tokens[0])
        target = int(tokens[1])
        if len(tokens) == 3:
            graph.createEdge(origin, target, int(tokens[2]))
        else:
            graph.createEdge(origin, target)

    return graph, embedding


def drawGraphVisualisation(name, tree, attributes=None, pixels=False):
    if not pydotAvailable:
        print("Error: this function requires the library pydot (which requires the software graphviz.")
        return

    # n = len(parent)
    G = pydot.Dot(graph_type='graph')

    nbPixels = tree.nbPixels

    for i in tree.iteratorFromPixelsToRoot(includePixels=pixels):
        node = pydot.Node(i)
        if i < nbPixels:
            node.set_shape("square")
        label = str(i)
        if attributes:
            labels = ",".join([str(tree.getAttribute(a)[i]) for a in attributes])
            label += "(" + labels + ")"
        node.set_label(label)
        G.add_node(node)

    for i in tree.iteratorFromPixelsToRoot(includePixels=pixels):
        par = tree[i]
        if par != -1:
            edge = pydot.Edge(par, i)
            edge.set_dir("forward")
            G.add_edge(edge)

    G.write_pdf(name)


def saveHiPy(data, filename):
    with open(filename, 'wb') as f:
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)


def loadHiPy(filename):
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    return data