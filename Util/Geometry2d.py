'''
Created on 3 juin 2015

@author: perretb
'''
from HiPy.Structures import Image, Embedding2dGrid, Adjacency2d4

def crop2d(image, xmin, xmax, ymin, ymax):
    '''
    limits included
    '''
    width=xmax-xmin+1
    height=ymax-ymin+1
    im = Image(width*height)
    im.embedding = Embedding2dGrid(width,height)
    for y in range(height):
        for x in range(width):
            im.setPixelWCS(image.getPixelWCS(x+xmin,y+ymin),x,y)
    return im
    

def interpolateX2(image,interpolation = max):
    sizeOri=[image.embedding.width,image.embedding.height]
    sizeDest=[((sizeOri[0])*2-1),((sizeOri[1])*2-1)]
    imDest= Image(sizeDest[0]*sizeDest[1],sizeDest)
    imDest.embedding = Embedding2dGrid(*sizeDest)
    for y in range(sizeDest[1]):
        for x in range(sizeDest[0]):
            modx=x%2
            mody=y%2
            xori=x//2
            yori=y//2
            if modx==0 and mody==0:
                    imDest.setPixelWCS(image.getPixelWCS(xori,yori), x, y )
            elif (modx==1 and mody==0):
                imDest.setPixelWCS(interpolation(image.getPixelWCS(xori,yori),image.getPixelWCS(xori+1,yori)), x, y )
            elif (modx==0 and mody==1):
                imDest.setPixelWCS(interpolation(image.getPixelWCS(xori,yori),image.getPixelWCS(xori,yori+1)), x, y )
            else:
                imDest.setPixelWCS(interpolation(image.getPixelWCS(xori,yori),image.getPixelWCS(xori+1,yori),image.getPixelWCS(xori,yori+1),image.getPixelWCS(xori+1,yori+1)), x, y )
    return imDest 

def getKhalimsky2FacesNeighbourValues(image, x, y):
    vals=[]
    width=image.embedding.width
    height=image.embedding.height
    if x%2==0 and y%2==0: # 0 face
        if x-1>=0:
            if y-1>=0:
                vals.append(image.getPixelWCS(x-1,y-1))
            if y+1<height:
                vals.append(image.getPixelWCS(x-1,y+1))
        if x+1<width:
            if y-1>=0:
                vals.append(image.getPixelWCS(x+1,y-1))
            if y+1<height:
                vals.append(image.getPixelWCS(x+1,y+1))

    elif y%2==1 and x%2==0: # vertical 1 face
        if x-1>=0:
            vals.append(image.getPixelWCS(x-1,y))
        if x+1<width:
            vals.append(image.getPixelWCS(x+1,y))

    elif y%2==0 and x%2==1: # horizontal 1 face
        if y-1>=0:
            vals.append(image.getPixelWCS(x,y-1))
        if y+1<height:
            vals.append(image.getPixelWCS(x,y+1))

    else : # 2 face !
        vals.append(image.getPixelWCS(x,y))
    return vals
        

def interpolationMedian(*args):
    l=len(args)
    op=sorted(args)
    
    if(l%2==1):
        return op[l//2]
    else:
        return (op[l//2-1]+op[l//2])//2



def interpolatePlainMapKhalimsky(image):
    sizeOri=[image.embedding.width,image.embedding.height]
    sizeDest=[((sizeOri[0])*2+1),((sizeOri[1])*2+1)]
    imDest= Image(sizeDest[0]*sizeDest[1])
    imDest.embedding = Embedding2dGrid(*sizeDest)
    imDest.adjacency=Adjacency2d4([sizeDest[0],sizeDest[1]])
    # 2 faces copy
    for y in range(1,sizeDest[1],2):
        for x in range(1,sizeDest[0],2):
            modx=x%2
            mody=y%2
            xori=x//2
            yori=y//2
            if modx==1 and mody==1:
                    imDest.setPixelWCS(image.getPixelWCS(xori,yori),x, y)

    # <2 faces
    for y in range(sizeDest[1]):
        for x in range(sizeDest[0]):
            if x%2 != 1 or y%2 != 1:
                vals=getKhalimsky2FacesNeighbourValues(imDest, x, y)
                minV = min(vals)
                maxV = max(vals)
                imDest.setPixelWCS([minV,maxV],x,y)

    return imDest

def reduceKhalimsky(image):
    sizeBig=[image.embedding.width,image.embedding.height]
    sizeDest=[(sizeBig[0])//2,(sizeBig[1])//2]
    res = Image(sizeDest[0]*sizeDest[1])
    res.embedding = Embedding2dGrid(*sizeDest)
    for y in range(sizeDest[1]):
        for x in range(sizeDest[0]):
            res.setPixelWCS(image.getPixelWCS(x*2+1,y*2+1),x,y)
    return res

def imagePadding(image, borderValue=0):
    sizeOri=[image.embedding.width,image.embedding.height]
    sizeDest=[((sizeOri[0])+2),((sizeOri[1])+2)]
    imDest= Image(sizeDest[0]*sizeDest[1],borderValue)
    imDest.embedding = Embedding2dGrid(*sizeDest)
    
    for y in range(sizeOri[1]):
        for x in range(sizeOri[0]):
            imDest.setPixelWCS(image.getPixelWCS(x,y),x+1,y+1)
    return imDest

def removeBorder(image):
    sizeOri=[image.embedding.width,image.embedding.height]
    sizeDest=[((sizeOri[0])-2),((sizeOri[1])-2)]
    imDest= Image(sizeDest[0]*sizeDest[1])
    imDest.embedding = Embedding2dGrid(*sizeDest)
    
    for y in range(sizeDest[1]):
        for x in range(sizeDest[0]):
            imDest.setPixelWCS(image.getPixelWCS(x+1,y+1),x,y)
    return imDest