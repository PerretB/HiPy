'''
Created on 3 juin 2015

@author: perretb
'''

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
        Rank[j] = Rank[j] + 1
    Par[i] = j
    return j,i

