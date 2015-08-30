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
Created on 3 juin 2015

@author: perretb
"""


def initialiseStructure(nbElem):
    ufParent = [i for i in range(nbElem)]
    ufRank = [0]*nbElem
    return ufParent, ufRank


def makeSet(ufParent, ufRank):
    ufParent.append(len(ufParent))
    ufRank.append(0)

# Union find
############################################################################## 
# find canonical node of elem with path compression
def findTarjan(elem, Par):
    i = elem
    while Par[i] != i:
        i = Par[i]
    while Par[elem] != i:
        temp = elem
        elem = Par[elem]
        Par[temp] = i
    return i


# union of nodes i and j with ranking, the chosen canonical node is returned
# and the successors of the other nodes are transfered to the canonical one
def unionTarjan(i, j, Par, Rank):
    if Rank[i] > Rank[j]:
        i, j = j, i
    elif Rank[i] == Rank[j]:
        Rank[j] += 1
    Par[i] = j
    return j, i
