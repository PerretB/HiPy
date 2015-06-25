# Copyright ESIEE (2015) 
# 
# benjamin.perret@esiee.fr
# 
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
Created on 17 juin 2015

@author: perretb
'''

from HiPy.Hierarchies.ComponentTree import constructComponentTree,\
    ComponentTreeType
from HiPy.Structures import Adjacency2d4

from HiPy.IO import * #@UnusedWildImport
from HiPy.Hierarchies import * #@UnusedWildImport
from HiPy.Processing.Attributes import * #@UnusedWildImport
from HiPy.Processing.Markov import * #@UnusedWildImport
from HiPy.Util.Histogram import * #@UnusedWildImport
from HiPy.Util.Color import * #@UnusedWildImport

def testShapeClassif():
    print('Reading image...')
    image = readImage("../samples/squares and circles.png", True)
    image.adjacency = Adjacency2d4(image.embedding.size)
    
    print('Pre-filtering image...')
    tree = constructComponentTree(image, ComponentTreeType.MaxTree)
    addAttributeArea(tree)
    addAttributeHeight(tree)
    imageFiltered = tree.reconstructImage(criterion=lambda n:tree.area[n]<5 or tree.height[n]<25)
    
    print('Tree construction...')
    tree = constructComponentTree(imageFiltered, ComponentTreeType.MaxTree)
    addAttributeArea(tree)
    addAttributePerimeter(tree)
    addAttributeSimpleMoments2d(tree)
    
    tree.addAttribute("observation")
    observation = tree.observation
    area = tree.area
    perimeter = tree.perimeter
    moments = tree.moments
    for i in tree.iteratorFromPixelsToRoot():
        comp = area[i]/(perimeter[i]**2)
        hu1 = computeHuInvariant(computeScaleInvarariantMoments2d(computeCentralMoments2d(moments[i])))[0]
        observation[i]=[comp,hu1]
        
        
    print('Markov model initialization...')
    markovModel=getInitialGuess(tree, 3)
    print('Markov model Estimation...')
    markovModel=MPMSegmentTree(markovModel,tree,maxIteration=10,verbose=False)
    
    reOrderLabels(markovModel, tree)
    
    res = tree.reconstructImage("estimLabel", lambda x:tree.estimLabel[x]==0)
    
    print('Saving label reconstructions to "Results/squares and circles classification.png"')
    saveImage(normalizeToByte(res), "Results/squares and circles classification.png")
 
    

def main():
    
    testShapeClassif()

if __name__ == '__main__':
    main()
