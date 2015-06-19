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
Created on 3 juin 2015

@author: perretb
'''

from HiPy.Structures import  Adjacency2d4
from HiPy.Util.Histogram import imageInverseGrayByte
from HiPy.Util.Geometry2d import imagePadding, reduceKhalimsky, removeBorder
from HiPy.Processing.Attributes import addAttributeArea
from HiPy.Hierarchies.TreeOfShape import constructTreeOfShapes
from HiPy.IO import readImage, saveImage


def testTreeIsomorphism(tree1,tree2):
    '''
    Test if tree1 and tree2 are isomorph assuming that leaves are ordered
    '''
    if(len(tree1)!=len(tree2) or tree1.nbPixels!=tree2.nbPixels):
        return False
    
    #both tree have same size so we need to find an injection m from the nodes of t1 to the nodes of t2 such that 
    # for any node n of t1 m(parent(n))=parent(m(n))
    mapT1T2 = [-1]*len(tree1)
    
    
    for i in range(len(mapT1T2)-1): #root is root !
        #pixel mapping is constant
        if i<tree1.nbPixels:
            mapT1T2[i]=i
        #parent(n)
        pT1=tree1[i]

        #parent(m(n))
        pT2=tree2[mapT1T2[i]]
        if mapT1T2[pT1]==-1:
            mapT1T2[pT1]=pT2
        elif mapT1T2[pT1]!=pT2:
            return False
        
        
    return True




def testSelfDuality():
    images =["monsters.png","macaws.png","mandrill.png","remotesensing1.png","stample.png","stample2.png","stample4.png","blobs-ndg.png","lennaGray256.png","blood1.png","detection_test.png","spot5.png"]
    for imName in images:
     
        im = readImage('../samples/' +imName)
        print("image " +imName)
        
        imInv = imageInverseGrayByte(im)
         
        tree1 =constructTreeOfShapes(im,None,False)
        tree2 =constructTreeOfShapes(imInv,None,False)
           
        if testTreeIsomorphism(tree1, tree2):
            print("Self duality verified")
        else:
            print("Arg: self duality broken")


def testAreaFilter():
    im = readImage('../samples/lenna.png')
    im = imagePadding(im, 0)
    
    im.adjacency = Adjacency2d4([im.embedding.width,im.embedding.height])
    tree= constructTreeOfShapes(im,None)
    addAttributeArea(tree)

    print("Reconstruction")
    reconstr = tree.reconstructImage("level",lambda x : tree.area[x]<=1000)
    
    reconstr=reduceKhalimsky(reconstr)
    reconstr=removeBorder(reconstr)
    
    resultName = 'Results/reconstructionAreaFilter_TreeOfShapes.png'
    print("Image save: " +  resultName)
    saveImage(reconstr, resultName)


def main():
    print("--- Grain filter on the tree of shapes")
    testAreaFilter()
    print("--- Experimental assesment of the self duality")
    testSelfDuality()
    
    
if __name__ == '__main__':
    main()

