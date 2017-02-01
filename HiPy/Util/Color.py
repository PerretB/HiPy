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
Created on 15 juin 2015

@author: Benjamin Perret
"""


def convertRGBtoLAB(image):
    """
    RGB in [0,255] -> LAB in [0,1]
    Observer. = 2°, Illuminant = D65
    """
    res = image.getCopy(False)
    if (not isinstance(image[0], list)) or len(image[0]) != 3:
        raise Exception("RGBtoLAB: input must be an image with 3 bands")

    for i in image.iterateOnPixels():
        res[i] = _RGBToLab(*image[i])
    return res

def _RGBToLab(r,g,b):
    # Observer. = 2°, Illuminant = D65
    var_R1 = (r / 255.)
    var_G1 = (g / 255.)
    var_B1 = (b / 255.)

    # RGB -> XYZ
    if var_R1 > 0.0404:
        var_R1 = ((var_R1 + 0.055) / 1.055) ** 2.4
    else:
        var_R1 /= 12.92
    if var_G1 > 0.04045:
        var_G1 = ((var_G1 + 0.055) / 1.055) ** 2.4
    else:
        var_G1 /= 12.92
    if var_B1 > 0.04045:
        var_B1 = ((var_B1 + 0.055) / 1.055) ** 2.4
    else:
        var_B1 /= 12.92

    var_R1 *= 100.
    var_G1 *= 100.
    var_B1 *= 100.

    var_X1 = (var_R1 * 0.4124 + var_G1 * 0.3576 + var_B1 * 0.1805) / 95.047   # ref_X =  95.047   Observer= 2°, Illuminant= D65
    var_Y1 = (var_R1 * 0.2126 + var_G1 * 0.7152 + var_B1 * 0.0722) / 100.000  # ref_Y = 100.000
    var_Z1 = (var_R1 * 0.0193 + var_G1 * 0.1192 + var_B1 * 0.9505) / 108.883  # ref_Z = 108.883

    # XYZ -> L*ab
    if var_X1 > 0.008856:
        var_X1 **= (1. / 3.)
    else:
        var_X1 = (7.787 * var_X1) + (16. / 116.)
    if var_Y1 > 0.008856:
        var_Y1 **= (1. / 3.)
    else:
        var_Y1 = (7.787 * var_Y1) + (16. / 116.)
    if var_Z1 > 0.008856:
        var_Z1 **= (1. / 3.)
    else:
        var_Z1 = (7.787 * var_Z1) + (16. / 116.)

    L1 = (116. * var_Y1) - 16.
    A1 = 500. * (var_X1 - var_Y1)
    B1 = 200. * (var_Y1 - var_Z1)

    # Normalization
    return [L1/100.0, (A1+86.185)/184.44, (B1+107.864)/202.35] # see _testRange() for normalization constants



def _testRange():
    minL = 99999
    mina = 99999
    minb = 99999
    maxL = -99999
    maxa = -99999
    maxb = -99999

    for r in range(256):
        print(r)
        for g in range(256):
            for b in range(256):
                [l,a,b] = _RGBToLab(r,g,b)
                minL = min(l, minL)
                maxL = max(l, maxL)
                mina = min(a, mina)
                maxa = max(a, maxa)
                minb = min(b, minb)
                maxb = max(b, maxb)

    print("L: " +str(minL) + " - " +str(maxL)) # L: 0.0 - 100.0
    print("a: " + str(mina) + " - " + str(maxa)) #a: -86.18463649762525 - 98.25421868616114
    print("b: " + str(minb) + " - " + str(maxb)) #b: -107.86368104495168 - 94.48248544644461


