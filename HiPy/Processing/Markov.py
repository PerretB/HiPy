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
import copy
import numpy as np
from sklearn import mixture
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
    
    def density(self,value):
        return self.sciPyPDF.pdf(value)
        

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
    
    for i in tree.iteratorFromRootToLeaves(True,False): #range(root-1,tree.nbPixels-1,-1):
        pLbl = label[tree[i]]
        lbl=getRandomValue(transitionProb[pLbl])
        label[i]=lbl
        observation[i]=classes[lbl].rand()
        
def gmmEstimate(data, nbClasses, estimator= mixture.GMM):
    '''
    Performs an GMM estimation on the data with nbClasses components.
       
    Raises an exception if the algorithm fails to converge
    '''
  
    clf = estimator(n_components=nbClasses, covariance_type='full')
    clf.fit(data)
    
    if not clf.converged_:
        raise Exception("GMM estimation failed, check your data and number of components!")

    classes=[]
    for i in range(nbClasses):
        classes.append(MultivariateGaussian(clf.means_[i], clf.covars_[i]))

    return classes


def computeLabelsPrior(markovModel ,tree):
    '''
    Compute p(x_s = w) the apriori probability that label x of site s equals w
    '''
    initialProb =markovModel.initialProb
    transitionProb=markovModel.transitionProb
    nbClasses=markovModel.numLabels
    prior = tree.prior
    
    prior[-1]=copy.copy(initialProb)
    for n in tree.iterateFromRootToLeaves(includeRoot=False):
        par = tree[n]
        for i in range(nbClasses):
            p = 0
            for j in range(nbClasses):
                p = p + prior[par][j]*transitionProb[j][i]
            prior[n][i]=p
    
def getInverseTransitionProb(markovModel, priorParent, i, priorChild, j):
    '''
    Compute the inverse transition probability p(x_s- = i | x_s = j) = p(x_s = j | x_s- = i)*p(x_s- = i)/p(x_s = j)
    if p(x_s = j)~=0 the result if 1.0/nClasses
     markovModelthe model
     priorParent : apriori probability for the parent x_s-
     i label of the parent
     priorChild : apriori probability for the child x_s
     j label of the child
    '''
    if priorChild[j]>epsilon:
        return markovModel.transitionProb[i][j]*priorParent[i]/priorChild[j]
    else:
        return 1.0/markovModel.numLabels


def computeNonNormalizedForwardProbability(markovModel, tree, node, masked):
    '''
    Compute p(x_n = i | y>n) the forward probability for all labels i on a non leaf node n. 
    The probabilities are not normalized ! The norm is returned instead.
    markovModel: the model
    tree: the tree
    node: the node where to compute probabilities (not a leaf !)
    masked: masked is false then the likelyhood given by the model is used otherwise we assume a 
            non informative likelyhood (uniform law on R^d) and thus p(y | x_n=i) = 1 for all y and all i
    '''
    norm=0
    observationAtNode = tree.observation[node]
    nbClasses = markovModel.numLabels
    prior = tree.prior
    pforward = tree.pforward
    children = tree.children
    
    # consistancy and infinity check!
    numericExceptions=[]
    
    for i in range(nbClasses):
        try:
            if prior[node][i]<epsilon:
                pforward[node][i]=0
            else:
                p=1.0
                for c in children[node]:
                    s=0.0
                    for j in range(nbClasses):
                        pTransitionInverse=getInverseTransitionProb(markovModel,prior[node],i,prior[c],j)
                        s = s + pTransitionInverse*pforward[c][j]
                    p = p * s / prior[node][i] #we do not factorize /mh.prior[i] in /mh.prior[i]^(nChildren-1) to avoid numerical issues when nChildren is large
                if not masked:
                    p = p * markovModel.classes[i].density(observationAtNode)
                p = p * prior[node][i]
                pforward[node][i] = p
                norm = norm + p
        except ZeroDivisionError:
            numericExceptions.append(i)
    
    if len(numericExceptions)>0:
        # we will distribute evenly the probabilities on the numeric exceptions
        print("Warning infinite value detected... trying to clean the mess.")
        nv = 1.0/nbClasses;
        for i in range(nbClasses):
            if i in numericExceptions:
                pforward[node][i] = 0
            else:
                pforward[node][i] = nv
        norm=1.0 # yahou at least we have a nice norm !
    
    return norm            
 
def forwardPass(markovModel, tree):
    '''
    Computes P(x_s | y_s> )
    the posterior prob of the labels conditionnally to all the observations of the descendants of the site s (s included)
    '''
    classes = markovModel.classes
    nbClasses = markovModel.numLabel
    
    computeLabelsPrior(markovModel , tree)
    prior = tree.prior
    observation = tree.observation
    pForward = tree.pForward
    
    for n in tree.iterateOnLeaves():
        norm=0
        y = observation[n]
        for i in range(nbClasses):
            v = prior[n][i]*classes[i].density(y)
            pForward[n][i] = v
            norm = norm + v
        #normalization
        if(norm<epsilon):
            print("Warning very low norm detected... maybe an outlier ?")
            #forget about likelyhood... @FIXME is this really what should be done ?
            for i in range(nbClasses):
                pForward[n][i]=prior[n][i]
        else:
            for i in range(nbClasses):
                pForward[n][i]=pForward[n][i]/norm
                
    for n in tree.iterateFromLeavesToRoot(False):
        norm=computeNonNormalizedForwardProbability(markovModel,tree, n, False);
        if norm > epsilon: #fine
            for i in range(nbClasses):
                pForward[n][i]=pForward[n][i]/norm
        else:
            print("Warning very low norm detected... maybe an outlier ?")
            norm=computeNonNormalizedForwardProbability(markovModel,tree, n, True);
            if norm > epsilon: #fine
                for i in range(nbClasses):
                    pForward[n][i]=pForward[n][i]/norm
            else:
                raise Exception("Numerical problem during forward pass, check your data and number of classes")


def backwardPass(markovModel, tree):
    '''
    Computes p(x_s=i|y) for all label i
    the marginal posterior distribution of the label x of s conditonnaly to the whole observation y
    '''
    nbClasses = markovModel.numLabel
    pBackward = tree.pBackward
    pForward = tree.pForward
    
    pBackward[-1]=copy.copy(pForward[-1])
    
    for n in tree.iterateFromRootToLeaves(False,False):
        
/**
     * Computes p(x_s=i|y) for all label i
     * the marginal posterior distribution of the label x of s conditonnaly to the whole observation y
     * @param m
     * @param tree
     */
    private static void backwardPass(Model m, ComponentTree<double []> tree)
    {
        double [] initialProb =m.initialProb;
        int nClasses=initialProb.length;
        for( ComponentNode<double []> n: tree.iterateFromRootToLeaf())
        {
            MarkovHelper mh=n.getAttributeValue(AttributeMarkovEstimation.class);
            if(n==tree.root)
            {
                mh.pbackward=mh.pforward.clone();
            }else{
                for(int i=0;i<nClasses;i++)
                {
                    double p=0.0;
                    
                    MarkovHelper mhp=n.parent.getAttributeValue(AttributeMarkovEstimation.class);
                    for(int j=0;j<nClasses;j++)
                    {
                        p+=getJointPost(m,tree,i,j,mh,mhp);    
                    }
                    mh.pbackward[i]=p;
                }
            }
        }
    }