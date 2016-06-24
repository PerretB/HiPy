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

import HiPy.Util.VMath
import HiPy.Structures

'''
Created on 9 juin 2015

@author: perretb
'''


def imageMap(image: "Image", function, marginal=False, bandNumber=None, inPlace=False):
    """
    Apply a pointwise  function to an image.
    This method assumes that the content of the image is homogeneous (either scalars or lists but not a mix
    of the twos).
    
    If the image elements are not list or marginal==False then the function is simply applied to each pixel
    of the image:
        For all pixels i,  result[i]=function(image[i]) 
    If the image elements are lists and marginal = True the function is applied to each band independently, 
        For all pixels i, for all bands b, result[i][b]=function(image[i][b])
    Moreover if a bandNumber is specified it will only be a applied to the specified band:
        For all pixels i, for all bands b, result[i][b]= function(image[i][b]) if b==bandNumber else image[i][b]
        
    If inPlace==False the result is put in a newly allocated image, in the other case the result is put in
    the input image.
    """
    if inPlace:
        imDest = image
    else:
        imDest = image.getCopy(False)

    vec = isinstance(image[0], list)
    if not vec or not marginal:
        for i in range(len(image)):
            imDest[i] = function(image[i])
    else:
        bands = len(image[0])
        if bandNumber is not None:
            for i in range(len(image)):
                imDest[i] = [function(image[i][band]) if band == bandNumber else image[i][band] for band in
                             range(bands)]
        else:
            for i in range(len(image)):
                imDest[i] = [function(image[i][band]) for band in range(bands)]
    return imDest


def imageInverseGrayByte(image: "Image"):
    return imageMap(image, lambda x: 255 - x, marginal=True)


def rescaleGray(image: "Image", minValue=0, maxValue=1, marginal=True):
    vmin, vmax = getMinMax(image)
    if vmin == vmax:
        HiPy.Structures.HiPyLogger.warning("rescaleGray: image is constant, cannot rescale its level.")
        return image
    return imageMap(image, lambda x: (maxValue - minValue) * (x - vmin) / (vmax - vmin) + minValue, marginal=marginal)


def toInt(image: "Image"):
    return imageMap(image, lambda x: int(round(x)), marginal=True)


def normalizeToByte(image: "Image"):
    im = rescaleGray(image, 0, 255)
    im = toInt(im)
    return im


def add(image: "Image", number, inPlace=False):
    return imageMap(image, lambda x: x + number, marginal=True, inPlace=inPlace)


def sub(image: "Image", number, inPlace=False):
    return add(image, -number, inPlace)


def mult(image: "Image", number, inPlace=False):
    return imageMap(image, lambda x: x * number, marginal=True, inPlace=inPlace)


def divide(image: "Image", number, inPlace=False):
    return imageMap(image, lambda x: x / number, marginal=True, inPlace=inPlace)


def getMinMax(image: "Image", band=None):
    '''
    Get min and max value over all band (if band==None) or a particular band
    '''
    vec = isinstance(image[0], list)
    if vec:
        if band is None:
            vmin = image[0][0]
            vmax = image[0][0]
            for i in range(len(image)):
                v = image[i]
                for vv in v:
                    if vv < vmin:
                        vmin = vv
                    elif vv > vmax:
                        vmax = vv
            return vmin, vmax
        else:
            vmin = image[0][band]
            vmax = image[0][band]

            for i in range(len(image)):
                v = image[i][band]
                if v < vmin:
                    vmin = v
                elif v > vmax:
                    vmax = v

            return vmin, vmax
    else:
        vmin = image[0]
        vmax = image[0]

        for i in range(len(image)):
            v = image[i]
            if v < vmin:
                vmin = v
            elif v > vmax:
                vmax = v

        return vmin, vmax


def combineBands(*args: "Image"):
    '''
    Takes several images and combine them in a single multiband images.
    All the input images must have the same dimension.
    
    If args = [image1, image2, ..., imageN)
    The result is, for all pixel i, res[i]=[image1[i], image2[i], ..., imageN[i]]
    '''
    res = args[0].getCopy(copyData=False)

    for i in res.iterateOnPixels():
        res[i] = [im[i] for im in args]

    return res


def median(*args: "Image"):
    '''
    Compute the median of several scalar images.
    All the input images must have the same dimension.
    
    If args = [image1, image2, ..., imageN)
    The result is, for all pixel i, res[i]=median(image1[i], image2[i], ..., imageN[i])
    '''
    res = args[0].getCopy(copyData=False)
    for i in res.iterateOnPixels():
        res[i] = HiPy.Util.VMath.medianV([im[i] for im in args])
    return res


def remapToInt(image):
    '''
    Takes a gray level image and apply the smallest strictly increasing gray level transform to the set of positive integer.

    If the input image contains n unique gray values: the transform t will have the following property:
    - t is a surjection to [0..n-1]
    - given v1, v2 in input image, v1 < v2 => t(v1) < t(v2)
    :param image: gray level image
    :return: gray level image with int values
    '''
    values = set()
    for v in image:
        values.add(v)
    values = sorted(values)
    myMap = dict()
    for i in range(len(values)):
        myMap[values[i]] = i
    return imageMap(image, lambda g: myMap[g])
