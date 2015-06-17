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
Created on 15 juin 2015

@author: perretb
'''

def add(image, number, inPlace=False):
    res = image.copy(False) if (not inPlace) else image
    for i in image.iterateOnPixels():
        res[i]=[x+number for x in image[i]]
    return res

def sub(image, number, inPlace=False):
    return add(image, -number, inPlace)

def mult(image, number, inPlace=False):
    res = image.copy(False) if (not inPlace) else image
    for i in image.iterateOnPixels():
        res[i]=[x*number for x in image[i]]
    return res

def divide(image, number, inPlace=False):
    res = image.copy(False) if (not inPlace) else image
    for i in image.iterateOnPixels():
        res[i]=[x/number for x in image[i]]
    return res

def getMinMax(image, band=None):
    '''
    Get min and max value over all band (if band==None) or a particular band
    '''
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

def rescaleColor(image, minValue=0, maxValue=1, inPlace=False):
    vmin, vmax=getMinMax(image)
    res=sub(image, vmin,inPlace)
    add(mult(divide(res,vmax-vmin,True),maxValue-minValue,True),minValue,True)
    return res
    
def toIntColor(image, inPlace=False):
    res = image.copy(False) if (not inPlace) else image
    for i in image.iterateOnPixels():
        res[i]=[int(round(x)) for x in image[i]]
    return res
    

def normalizeToByteColor(image):
    im=rescaleColor(image, 0, 255) 
    im=toIntColor(im)
    return im    

def convertRGBtoLAB(image):
    '''
    RGB in [0,255] -> LAB in [0,1]
    '''
    res = image.copy(False)
    
    for i in image.iterateOnPixels():
        v=image[i]
        r1=v[0]
        g1=v[1]
        b1=v[2]
        
        var_R1 = ( r1 / 255. )     
        var_G1 = ( g1 / 255. )    
        var_B1 = ( b1 / 255. )  
        
        if ( var_R1 > 0.04045 ): 
            var_R1 = (( var_R1 + 0.055 ) / 1.055)**2.4
        else:
            var_R1 = var_R1 / 12.92
        if ( var_G1 > 0.04045 ):
            var_G1 = ( ( var_G1 + 0.055 ) / 1.055)**2.4
        else:
            var_G1 = var_G1 / 12.92
        if ( var_B1 > 0.04045 ) :
            var_B1 = ( ( var_B1 + 0.055 ) / 1.055)**2.4
        else :
            var_B1 = var_B1 / 12.92

        var_R1 = var_R1 * 100.
        var_G1 = var_G1 * 100.
        var_B1 = var_B1 * 100.

        #Observer. = 2°, Illuminant = D65
        var_X1 = (var_R1 * 0.4124 + var_G1 * 0.3576 + var_B1 * 0.1805) / 95.047          #ref_X =  95.047   Observer= 2°, Illuminant= D65
        var_Y1 = (var_R1 * 0.2126 + var_G1 * 0.7152 + var_B1 * 0.0722) / 100.000         #ref_Y = 100.000
        var_Z1 = (var_R1 * 0.0193 + var_G1 * 0.1192 + var_B1 * 0.9505) / 108.883         #ref_Z = 108.883

        if ( var_X1 > 0.008856 ) :
            var_X1 = var_X1**(1./3.)
        else:
            var_X1 = ( 7.787 * var_X1 ) + ( 16. / 116. )
        if ( var_Y1 > 0.008856 ): 
            var_Y1 = var_Y1**(1./3.)
        else:
            var_Y1 = ( 7.787 * var_Y1 ) + ( 16. / 116. )
        if ( var_Z1 > 0.008856 ) :
            var_Z1 = var_Z1**(1./3.)
        else:
            var_Z1 = ( 7.787 * var_Z1 ) + ( 16. / 116. )

        L1 = ( 116. * var_Y1 ) - 16.
        A1 = 500. * ( var_X1 - var_Y1 )
        B1 = 200. * ( var_Y1 - var_Z1 )
    
        res[i] = (L1,A1,B1)
    
    return res


