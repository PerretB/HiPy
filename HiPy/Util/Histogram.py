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

from HiPy.Util.VMath import medianV
'''
Created on 9 juin 2015

@author: perretb
'''

def imageMap(image, function, marginal=False, bandNumber=None, inPlace=False ):
    '''
    Apply apointwise  function to an image.
    This method assumes that the content of the image is homogeneous (either scalars or lists but not a mix of the twos).
    
    If the image elements are not list or marginal==False then the function is simply applied to each pixel of the image:
        For all pixels i,  result[i]=function(image[i]) 
    If the image elements are lists and marginal = True the function is applied to each band independently, 
        For all pixels i, for all bands b, result[i][b]=function(image[i][b])
    Moreover if a bandNumber is specified it will only be a applid to the specified band:
        For all pixels i, for all bands b, result[i][b]= function(image[i][b]) if b==bandNumber else image[i][b]
        
    If inPlace==False the result is put in a newly allocated image, in the other case the result is put in the input image.
    '''
    imDest = image.copy(False)
    
    vec= isinstance(image[0], list)
    if not vec or not marginal:
        for i in range(len(image)):
            imDest[i] = function(image[i])
    else:
        bands=len(image[0])
        if bandNumber!=None:
            for i in range(len(image)):
                imDest[i] = imDest[i] = [function(image[i][band]) if band==bandNumber else image[i][band] for band in range(bands) ]
        else: 
            for i in range(len(image)):
                imDest[i] = [function(image[i][band]) for band in range(bands) ]
    return imDest

def imageInverseGrayByte(image):
    return imageMap(image, lambda x:255-x, marginal=True)

def rescaleGray(image, minValue=0, maxValue=1, marginal=True):
    vmin, vmax=getMinMax(image)  
    return imageMap(image, lambda x:(maxValue-minValue)*(x-vmin)/(vmax-vmin)+minValue,marginal=marginal)

def toInt(image):
    return imageMap(image, lambda x:int(round(x)), marginal=True)

def normalizeToByte(image):
    im=rescaleGray(image, 0, 255) 
    im=toInt(im)
    return im

def add(image, number, inPlace=False):
    return imageMap(image, lambda x: x+number,marginal=True,inPlace=inPlace)

def sub(image, number, inPlace=False):
    return add(image, -number, inPlace)

def mult(image, number, inPlace=False):
    return imageMap(image, lambda x: x*number,marginal=True,inPlace=inPlace)

def divide(image, number, inPlace=False):
    return imageMap(image, lambda x: x/number,marginal=True,inPlace=inPlace)

def getMinMax(image, band=None):
    '''
    Get min and max value over all band (if band==None) or a particular band
    '''
    vec= isinstance(image[0], list)
    if vec:
        if band==None:
            vmin=image[0][0]
            vmax=image[0][0]
            for i in range(len(image)):
                v=image[i]
                for vv  in v:
                    if vv<vmin:
                        vmin=vv
                    elif vv>vmax:
                        vmax=vv
            return vmin, vmax
        else:
            vmin=image[0][band]
            vmax=image[0][band]
            
            for i in range(len(image)):
                v=image[i][band]
                if v<vmin:
                    vmin=v
                elif v>vmax:
                    vmax=v
                    
            return vmin, vmax
    else:
        vmin=image[0]
        vmax=image[0]
        
        for i in range(len(image)):
            v=image[i]
            if v<vmin:
                vmin=v
            elif v>vmax:
                vmax=v
                
        return vmin, vmax

def combineBands(*args):
    '''
    Takes several images and combine them in a single multiband images.
    All the input images must have the same dimension.
    
    If args = [image1, image2, ..., imageN)
    The result is, for all pixel i, res[i]=[image1[i], image2[i], ..., imageN[i]]
    '''  
    res = args[0].copy(copyData=False)
    
    for i in res.iterateOnPixels():
        res[i]=[im[i] for im in args]
        
    return res

def median(*args):
    '''
    Compute the median of several scalar images.
    All the input images must have the same dimension.
    
    If args = [image1, image2, ..., imageN)
    The result is, for all pixel i, res[i]=median(image1[i], image2[i], ..., imageN[i])
    '''  
    res = args[0].copy(copyData=False)
    for i in res.iterateOnPixels():
        res[i]=medianV([im[i] for im in args])
    return res
    






    