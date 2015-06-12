# Copyright ESIEE (2015) 
# 
# benjamin.perret@esiee.fr
# 
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
Created on 11 juin 2015

@author: perretb
'''

from HiPy.IO import * #@UnusedWildImport
from HiPy.Hierarchies.ComponentTree import * #@UnusedWildImport
from HiPy.Hierarchies.TreeOfShape import * #@UnusedWildImport
from HiPy.Util.Histogram import * #@UnusedWildImport
from HiPy.Structures import Adjacency2d4
from HiPy.Util.Geometry2d import crop2d

def dummyDemoShapping():
    print("Reading image...")
    image = readImage("../samples/blood1.png")
    image.adjacency = Adjacency2d4(image.embedding.size)
    
    print("Building tree...")
    tree=constructTreeOfShapes(image,verbose=False)
    addAttributeArea(tree)
    addAttributeChildren(tree)

    print("Building tree of tree on area")
    area=prepareForShapping(tree,tree.area)
    tree2 = constructComponentTree(area,verbose=False)

    print("Filtering tree of tree on area based on level")
    tree2.filterDirect(lambda _,x:tree2.level[x]<5000)
    print("Shaping first tree based on previous filter result")
    res = tree.reconstructImage("level",shapingCriterion(tree2))
    
    print("Filtering first tree based on area")
    res2 = tree.reconstructImage("level", lambda x: tree.area[x]>=5000)
    
    if res.equals(res2):
        print("Hurray: its the same")
    else:
        print("Yeark: does not work")
    

class prettyfloat(float):
    def __repr__(self):
        return "%10.2f" % self

def demoNonIncreasingFilter(): 
    print("Reading image...")
    image = readImage("../samples/blood1.png")
    image = crop2d(image,100,200,100,200)
    image.adjacency = Adjacency2d4(image.embedding.size)
    
    print("Building tree...")
    tree=constructTreeOfShapes(image,verbose=False)
    
#     print(tree.nbLeaves)
#     print(len(tree.data))
#     print(list(map(prettyfloat, tree.data)))
    

    
    addAttributeCompactness(tree)
#     print(list(map(prettyfloat, tree.compactness.data)))
#     print(list(map(prettyfloat, tree.area.data)))
    
    
    tree.filterDirect(lambda _,x:tree.area[x]<20)
    updateAttributeAfterFiltering(tree, "compactness")
    print(list(map(prettyfloat, tree.compactness.data)))
#     nbLeaves = tree.nbLeaves
#     print(list(map(prettyfloat, tree.area.data[nbLeaves:-1])))
#     print(list(map(prettyfloat, tree.perimeter.data[nbLeaves:-1])))
#     print(list(map(prettyfloat, tree.compactness.data[nbLeaves:-1])))
    
    print("Building tree of tree on Compactness...")
    compactness=prepareForShapping(tree,tree.compactness)
    tree2 = constructComponentTree(compactness)
    addAttributeHighest(tree2)
    addAttributeDynamics(tree2,"highest")
    addAttributeExtrema(tree2)
    nbLeaves=tree2.nbLeaves
    nbNodes=len(tree2)
    print(list(map(prettyfloat, tree2.level.data[nbLeaves:nbNodes])))
    print(list(map(prettyfloat, tree2.highest.data[nbLeaves:nbNodes])))
    print(list(map(prettyfloat, tree2.dynamics.data[nbLeaves:nbNodes])))
    print(list(map(prettyfloat, tree2.extrema.data[nbLeaves:nbNodes])))
    print("Filtering tree of tree on Compactness based on dynamics..")
    tree2.filterDirect(lambda _,x:not(tree2.level[x]>0.5))# and tree2.extrema[x]))# 
    print("Shaping first tree based on previous filter result...")
    res = tree.reconstructImage("level",lambda x: shapingCriterion(tree2)(x) and tree.area[x]>10)
    
    print("Saving result...")
    saveImage(res, "res.png")

    
def main():
    demoNonIncreasingFilter()
    #dummyDemoShapping()
    
if __name__ == '__main__':
    main()