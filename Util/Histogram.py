'''
Created on 9 juin 2015

@author: perretb
'''

def imageMap(image, function):
    imDest = image.copy(False)
    for i in range(len(image)):
        imDest[i] = function(image[i])
    return imDest

def imageInverseGrayByte(image):
    return imageMap(image, lambda x:255-x)

def getMinMax(image):
    vmin=image[0]
    vmax=image[0]
    
    for i in range(len(image)):
        v=image[i]
        if v<vmin:
            vmin=v
        elif v>vmax:
            vmax=v
            
    return vmin, vmax

def rescaleGray(image, minValue=0, maxValue=1):
    vmin, vmax=getMinMax(image)  
    return imageMap(image, lambda x:(maxValue-minValue)*(x-vmin)/(vmax-vmin)+minValue)



def toInt(image):
    return imageMap(image, lambda x:int(round(x)))

def normalizeToByte(image):
    im=rescaleGray(image, 0, 255) 
    im=toInt(im)
    return im
    