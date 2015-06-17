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
from HiPy.Hierarchies.ComponentTree import constructComponentTree,\
    ComponentTreeType
from HiPy.Structures import Adjacency2d4

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
from HiPy.IO import * #@UnusedWildImport
from HiPy.Hierarchies import * #@UnusedWildImport
from HiPy.Processing.Markov import * #@UnusedWildImport
from HiPy.Util.Histogram import * #@UnusedWildImport
from HiPy.Util.Color import * #@UnusedWildImport

def simuDirect():
    print('Reading image...')
    image = readImage("../samples/blobs-ndg.png", True)
    image.adjacency = Adjacency2d4(image.embedding.size)
    
    print('Constructing tree...')
    maxTree = constructComponentTree(image, ComponentTreeType.MaxTree)
    
    print('Defining direct model...')
    initialProb= [0.2,0.2,0.6]
    transitionProb=[[0.9,0.08,0.02],[0.01,0.95,0.04],[0.05,0.05,0.9]]
    class1 = MultivariateGaussian([0.9,0.2,0.1],[[0.2,0.001,0.0],[0.001,0.1,0.0],[0.0,0.0,0.1]])
    class2 = MultivariateGaussian([0.3,0.7,0.3],[[0.1,0.0,0.02],[0.0,0.1,0.0],[0.02,0.0,0.2]])
    class3 = MultivariateGaussian([0.0,0.2,0.7],[[0.2,0.0,0.0],[0.0,0.1,0.0],[0.0,0.0,0.1]])
    classes=[class1,class2,class3]
    model = MarkovModel(initialProb,transitionProb,classes)
   
    

    print('Simulating direct model...')
    simulateDirectModel(model, maxTree)
    
    print('Generating images...')
    imLabel=maxTree.reconstructImage("label")
    imObservation=maxTree.reconstructImage("observation")
    

    saveImage(normalizeToByte(imLabel), "Results/Markov simulation labels.png")
    saveImage(normalizeToByteColor(imObservation), "Results/Markov simulation observations.png")
    
    

def main():
    
    simuDirect()

if __name__ == '__main__':
    main()
