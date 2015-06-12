'''
Created on 3 juin 2015

@author: perretb
'''


from HiPy.Structures import Image, Embedding2dGrid, AdjacencyEdgeWeightedGraph

# for image I/O
try:
    from PIL import Image as ImagePIL
except ImportError:
    print("Error: PIL library not available !")
    
    

def readImage(filename, grayScale=True): 
    image= ImagePIL.open(filename)
    size = image.size
    coord = lambda x: (x % size[0], x // size[0])
    pixdata = image.load()
    im = Image(size[0]*size[1])
    if grayScale and isinstance(pixdata[(0,0)], (list, tuple)):
        if len(pixdata[(0,0)])>3:
            print("Warning: Image.Read, image has more than 3 channels (alpha component?), one or more channels will be ignored during the grayscale conversion.")
        for i in range(size[0]*size[1]):
            im.data[i]=sum(pixdata[coord(i)][0:3])//3
    else:
        for i in range(size[0]*size[1]):
            im.data[i]=(pixdata[coord(i)])
    im.embedding = Embedding2dGrid(size[0],size[1])
    return im

def saveImage(image,filename):
    width=image.embedding.width;
    size=[width,image.embedding.height]
    if isinstance(image.data[0], (list, tuple)):
        im = ImagePIL.new("RGB", size, 0)
    else:
        im = ImagePIL.new("L", size, 0)
    pix = im.load()
    
    for i in range(len(image)):
        pix[i%width,i//width]=image.data[i]
    im.save(filename, "PNG")
    
    

# Read a graph (AdjacencyEdgeWeightedGraph) from ascii file.
# format (weights are optionnale, vertices are numbered from 0 to numberOfVertices-1):
# numberOfVertices numberOfEdges
# sourceVertexOfEdge_1 destinationVertexOfEdge_1 weightOfEdge_1
# sourceVertexOfEdge_2 destinationVertexOfEdge_2 weightOfEdge_2
# ...
# sourceVertexOfEdge_NumberOfEdges destinationVertexOfEdge_NumberOfEdges weightOfEdge_NumberOfEdges
def readGraph(filename):
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

# Save a graph (AdjacencyEdgeWeightedGraph) in a file as plain text
def saveGraph(graph, filename):
    out = open(filename, "w")
    out.write(str(graph["nbPoints"]) + " " + str(len(graph["edges"])) + "\n")
    for e in graph["edges"]:
        line=" ".join([str(v) for v in e])
        out.write(line)
        out.write("\n")
    out.close()