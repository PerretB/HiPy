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
Created on 26 juin 2015

@author: perretb
'''

from math import * #@UnusedWildImport

def mapV(vector, function):
    return [function(v) for v in vector]

def mapVV(v1, v2, function):
    return [function(v1[i],v2[i]) for i in range(len(v1))]

def diffV(v1,v2):
    return mapVV(v1,v2,lambda x,y:x-y)

def multV(v1,v2):
    return mapVV(v1,v2,lambda x,y:x*y)

def addV(v1,v2):
    return mapVV(v1,v2,lambda x,y:x+y)

def divV(v1,v2):
    return mapVV(v1,v2,lambda x,y:x/y)

def norm(v):
    return sqrt(sum(multV(v,v)))

def mult(v):
    res=1
    for val in v:
        res*=val
    return res

def euclideanDistance(v1,v2):
    return norm(diffV(v1,v2))

def medianV(v):
    s=sorted(v)
    nbE = len(s)
    nbEd2 = nbE//2
    if nbE%2==1:
        return s[nbEd2]
    else:
        return (s[nbEd2]+s[nbEd2+1])/2