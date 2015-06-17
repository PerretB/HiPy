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

import sys
import random
import numpy as np
from scipy.stats import multivariate_normal
# system dependant epsilon: why not !
epsilon = sys.float_info.epsilon*10.0

# determinist random generator
rnd = random.Random() 
rnd.seed(1)

class MultivariateGaussian(object):
    def __init__(self, mean, varCovarMatrix):
        self.sciPyPDF= multivariate_normal(mean, varCovarMatrix)
        self.mean=mean
        self.varCovarMatrix=varCovarMatrix
        
    def rand(self):
        # it's always a pleasure to interract with numpy from pop
        #return np.random.multivariate_normal(self.mean,self.varCovarMatrix,1)[0].tolist()
        return self.sciPyPDF.rvs().tolist()
        

class MarkovModel(object):
    def __init__(self,initialProb, transitionProb, classes):
        self.initialProb = initialProb
        self.transitionProb = transitionProb
        self.classes = classes
        self.numLabels=len(initialProb)
        
    

def getRandomValue(discretePDF):
    '''
    Get a random value from a discrete pdf stored in an array
    '''
    s=discretePDF[0]
    v=rnd.random()
    i=0
    while v>s:
        i = i+1
        if i>=len(discretePDF):
            return i-1
        s = s+discretePDF[i]
    return i


def simulateDirectModel(model, tree, attributeLabelName="label",attributeObservationName="observation"):
    '''
    Simulate the given model on the given tree, add two attributes to the tree:
    - simulated labels in the attribute "attributeLabelName" 
    - simulated observations in the attribute "attributeObservationName" 
    '''
    transitionProb=model.transitionProb;
    classes=model.classes;
    label=tree.addAttribute(attributeLabelName,0,True)
    observation=tree.addAttribute(attributeObservationName,0,True)
    
    root=len(tree)-1
    lbl=getRandomValue(model.initialProb)
    label[root]=lbl
    observation[root]=classes[lbl].rand()
    
    for i in range(root-1,tree.nbLeaves-1,-1):
        pLbl = label[tree[i]]
        lbl=getRandomValue(transitionProb[pLbl])
        label[i]=lbl
        observation[i]=classes[lbl].rand()
        

    


