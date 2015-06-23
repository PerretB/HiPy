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
import time 
from HiPy.Processing.Attributes import addAttributeChildrenLogical
from math import * #@UnusedWildImport


def timeMili():
    return time.time()*1000.0
    
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
    
    def __str__(self):
        stri ="Mean:\t\t" + str(self.mean) + "\n"
        stri+="Covar:"
        for i in self.varCovarMatrix:
            stri+= "\t\t"+ str(i) + "\n"
        return stri

class MarkovModel(object):
    def __init__(self,initialProb, transitionProb, classes):
        self.initialProb = initialProb
        self.transitionProb = transitionProb
        self.classes = classes
        self.numLabels=len(initialProb)
    
    def __str__(self):
        stri=[]
        stri.append("Initial Probability:\n\t\t" + str(self.initialProb) + "\n")
        stri.append("Transition Probability:\n")
        for i in self.transitionProb:
            stri+= "\t\t"+ str(i) + "\n"
        for i in range(self.numLabels):
            stri.append("Classe " +str(i) +":\n")
            stri.append(str(self.classes[i]))
        return "".join(stri)
    

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

def MPMSegmentTree(markovModel, tree, convergenceRatio=0.01, maxIteration=10, cleanTemporaryAttributes=True, verbose=False):
    '''
    Performs an unsupervised iterative MPM segmentation of the tree
    The initial model parameters
    '''
    nbSites = tree.nbNodes()
    nbClasses = markovModel.numLabels
    if verbose:
        print("Number of sites (nodes) :" + str(nbSites))
    changeRate=1.0
    count=1
    prepareTreeForMPMSegment(tree, nbClasses)
    if verbose:
        print("Initial Markov model:")
        print(markovModel)
    
    while changeRate>convergenceRatio and count<=maxIteration:
        if verbose:
            print("Iteration: " +str(count))
            t1=timeMili()
            
        if count>1:
            markovModel=updateModel(markovModel, tree)
            if verbose:
                print("New Markov model:")
                print(markovModel)
        
        browseTree(markovModel, tree)
        changeRate = finalMPMSegment(tree)/nbSites

        if verbose:
            t2=timeMili()
            print("End of iteration: " +str(count) + " time: " + str(t2-t1) +" ms")
            print("Modification rate: " + str(changeRate*100) +"%")
            
        count+=1
    if verbose:
        print("Final Markov model:")
        print(markovModel)
    if cleanTemporaryAttributes:
        cleanTreeAfterMPMSegment(tree)
    return markovModel



def prepareTreeForMPMSegment(tree, nbClasses):
    addAttributeChildrenLogical(tree)
    tree.addAttribute("pForward",[0]*nbClasses)
    tree.addAttribute("pBackward",[0]*nbClasses)
    tree.addAttribute("prior",[0]*nbClasses)
    tree.addAttribute("estimLabel",-1)

def cleanTreeAfterMPMSegment(tree):
    tree.deleteAttribute("pForward")
    tree.deleteAttribute("pBackward")
    tree.deleteAttribute("prior")
    tree.deleteAttribute("childrenLogical")


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
    for n in tree.iteratorFromRootToLeaves(includeRoot=False):
        par = tree[n]
        for i in range(nbClasses):
            p = 0
            for j in range(nbClasses):
                p += prior[par][j]*transitionProb[j][i]
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
    pForward = tree.pForward
    children = tree.childrenLogical
    
    # consistancy and infinity check!
    numericExceptions=[]
    
    for i in range(nbClasses):
        try:
            if prior[node][i]<epsilon:
                pForward[node][i]=0
            else:
                p=1.0
                for c in children[node]:
                    s=0.0
                    pForwardC=pForward[c]
                    for j in range(nbClasses):
                        pTransitionInverse=getInverseTransitionProb(markovModel,prior[node],i,prior[c],j)
                        s = s + pTransitionInverse*pForwardC[j]
                    p = p * s / prior[node][i] #we do not factorize /mh.prior[i] in /mh.prior[i]^(nChildren-1) to avoid numerical issues when nChildren is large
                if not masked:
                    p = p * markovModel.classes[i].density(observationAtNode)
                p = p * prior[node][i]
                pForward[node][i] = p
                norm = norm + p
        except ZeroDivisionError:
            numericExceptions.append(i)
    
    if len(numericExceptions)>0:
        # we will distribute evenly the probabilities on the numeric exceptions
        print("Warning infinite value detected... trying to clean the mess.")
        nv = 1.0/nbClasses;
        for i in range(nbClasses):
            if i in numericExceptions:
                pForward[node][i] = 0
            else:
                pForward[node][i] = nv
        norm=1.0 # yahou at least we have a nice norm !
    
    return norm            
 
def forwardPass(markovModel, tree):
    '''
    Computes P(x_s | y_s> )
    the posterior prob of the labels conditionnally to all the observations of the descendants of the site s (s included)
    '''
    classes = markovModel.classes
    nbClasses = markovModel.numLabels
    
    computeLabelsPrior(markovModel , tree)
    prior = tree.prior
    observation = tree.observation
    pForward = tree.pForward
    
    for n in tree.iteratorOnLeaves():
        norm=0
        
        y = observation[n]
        for i in range(nbClasses):
            v = prior[n][i]*classes[i].density(y)
            pForward[n][i] = v
            norm +=  v
        #normalization
        if(norm<epsilon):
            print("Warning very low norm detected on a leave... maybe an outlier ?")
            #forget about likelyhood... @FIXME is this really what should be done ?
            for i in range(nbClasses):
                pForward[n][i]=prior[n][i]
        else:
            for i in range(nbClasses):
                pForward[n][i]=pForward[n][i]/norm
       
    for n in tree.iteratorFromLeavesToRoot(False):
        norm=computeNonNormalizedForwardProbability(markovModel,tree, n, False);
        if norm > epsilon: #fine
            for i in range(nbClasses):
                pForward[n][i]=pForward[n][i]/norm
        else:
            print("Warning very low norm detected on node " +str(n) + "... maybe an outlier ?")
            norm=computeNonNormalizedForwardProbability(markovModel,tree, n, True);
            if norm > epsilon: #fine
                for i in range(nbClasses):
                    pForward[n][i]=pForward[n][i]/norm
            else:
                raise Exception("Numerical problem during forward pass, check your data and number of classes")


def getJointPost(markovModel,tree,n,p,i,j):
    '''
    P(x_s,x_s- |y)
    the joint posterior distribution of the label of s and its parent s- knowing all the observation y
    @param markovModel Current model
    @param tree 
    @param n node xs
    @param p node xs-
    @param i class of xs 
    @param j class of xs-
    @return
    '''
    s=0
    transitionProb = markovModel.transitionProb
    nbClasses = markovModel.numLabels
    
    prior=tree.prior
    pForward = tree.pForward
    pBackward = tree.pBackward
    
    if prior[n][i]<epsilon: #@FIXME is this really what should be done !
        return 0.0
    
    for k in range(nbClasses):
        if prior[n][k]>epsilon:
            pTransitionInverse=getInverseTransitionProb(markovModel, prior[p],j, prior[n],k)
            s+=pTransitionInverse*pForward[n][k]
    
    if s<epsilon: #@FIXME is this really what should be done !
        return 0.0
    
    pj = transitionProb[j][i]*prior[p][j]*pForward[n][i]*pBackward[p][j]/prior[n][i]
    return pj/s
    
  
   

def backwardPass(markovModel, tree):
    '''
    Computes p(x_s=i|y) for all label i
    the marginal posterior distribution of the label x of s conditonnaly to the whole observation y
    '''
    nbClasses = markovModel.numLabels
    pBackward = tree.pBackward
    pForward = tree.pForward
    
    pBackward[-1]=copy.copy(pForward[-1])

    for n in tree.iteratorFromRootToLeaves(True,False):
        par=tree[n]
        for i in range(nbClasses):
            p=0.0
            for j in range(nbClasses):
                p+=getJointPost(markovModel,tree,n,par,i,j)

            pBackward[n][i]=p


def browseTree(markovModel,tree):
    '''
    browse forward and then backward
    '''
    forwardPass(markovModel, tree)
    backwardPass(markovModel, tree)


def finalMPMSegment(tree):
    '''
    Gives the label of max marginal posterior prob at each node : arg (w) max p(x_s=w|y)
    Return the number of node whose labels have changed
    '''
    diff=0
    pBackward = tree.pBackward
    estimLabel = tree.estimLabel
    for n in tree.iteratorFromLeavesToRoot():
        lbl = pBackward[n].index(max(pBackward[n]))
        if lbl != estimLabel[n]:
            diff += 1
        estimLabel[n] = lbl
    return diff

def updateInitialProb(markovModel, newMarkovModel, tree):
    '''
    Update initial probabilities  pi_w = p(xr=w)
    '''
    s=0
    pBackwardRoot = tree.pBackward[-1]
    for i in range(markovModel.numLabels):
        v = pBackwardRoot[i]
        if s+v > 1.0:
            v = max(1.0-s,0)
        s += v
        newMarkovModel.initialProb[i]=v
    
    if abs(s-1)>0.001:
            raise Exception("Update initial probability failed: integral different of 1 ! " + str(s));


def updateTransitionProb(markovModel, newMarkovModel, tree):
    '''
    Update transition probabilities a_ij = p( x=j | x-=i) (x- := parent of x)
    '''
    nbClasses = markovModel.numLabels
    pBackward = tree.pBackward
    tmp = [0] * nbClasses
    transitionProb = newMarkovModel.transitionProb
    
    for n in tree.iteratorFromRootToLeaves(includeRoot=False):
        par = tree[n]
        for i in range(nbClasses):
            tmp[i]+=pBackward[par][i]
            for j in range(nbClasses):
                transitionProb[i][j]+=getJointPost(markovModel, tree, n, par, j, i)
                
    for i in range(nbClasses):
        if tmp[i]<epsilon:
            raise Exception("Update transition probability failed: sum to low for class "+ str(i))
        s=0
        for j in range(nbClasses):
            v = transitionProb[i][j]/tmp[i]
            if s+v>1.0:
                v=max(1.0-s,0)
            s+=v
            transitionProb[i][j]=v
        if abs(s-1.0)>0.001:
            raise Exception("Update transition probability failed: integral different of 1 ! " + str(s));

def updateClassGaussian(markovModel, newMarkovModel, tree, classe):
    '''
    Update parametres of gaussian likelihood for classe "classe"
    '''
    dim = len(markovModel.classes[0].mean)
    observation = tree.observation
    pBackward = tree.pBackward
    
    mean = [0]*dim
    covar = [[0]*dim for _ in range(dim)] 

    #mean update
    s=0
    for n in tree.iteratorFromLeavesToRoot():
        y = observation[n]
        marginalPost = pBackward[n][classe]
        s+=marginalPost
        for i in range(dim):
            mean[i]+=y[i]*marginalPost
            
    if s<epsilon:
        raise Exception("Update gaussian parameters failed: sum to low for class " + str(classe))
    for i in range(dim):
            mean[i]/=s
            
    #covar update
    tmp =[0]*dim
    for n in tree.iteratorFromLeavesToRoot():
        y = observation[n]
        marginalPost = pBackward[n][classe]
        for i in range(dim):
            tmp[i]= y[i]-mean[i]
        for i in range(dim):
            for j in range(dim):
                covar[i][j]+=marginalPost*tmp[i]*tmp[j]
                
    for i in range(dim):
        for j in range(dim):
            covar[i][j]/=s
            
    newMarkovModel.classes[classe] = MultivariateGaussian(mean,covar)
    
def updatePDFParamGaussian(markovModel, newMarkovModel, tree):
    '''
    Update parametres of gaussian likelihood for all classes
    '''
    for i in range(markovModel.numLabels):
        updateClassGaussian(markovModel, newMarkovModel, tree, i)

def updateModel(markovModel, tree):
    '''
    Update all model parameters with EM
          - initial pobabilities pi_w = p(xr=w)
          - transition probabilities a_ij = p( x=j | x-=i) (x- := parent of x)
          - parameters of the likelyhood (means and covariance matrices of the gaussians)
          
    Return the new updated model
    '''
    nbClasses = markovModel.numLabels
    initialProb = [0]*nbClasses
    transitionProb = [[0]*nbClasses for _ in range(nbClasses)] 
    classes = [None]*nbClasses
    newMarkovModel = MarkovModel(initialProb, transitionProb, classes)
    updateInitialProb(markovModel,newMarkovModel,tree)
    updateTransitionProb(markovModel,newMarkovModel,tree)
    updatePDFParamGaussian(markovModel, newMarkovModel, tree)
    return newMarkovModel


def getInitialGuess(tree, nbClasses, maxSamples=40000):
    '''
    Fits an initial model with nClasses classes to the observation Y: 
        - the initial probability pi_i is equal for all i: pi_i=1/nClasses
        - the transition probability a_ij is 0.75 if i=j and 0.25/(nClasses-1) otherwise
        - parameters of the gaussian likelyhood are estimated by fitting a gaussian mixture model on the (not) the whole data 
        
    @todo: support maxSamples
    '''
    initialProb = [1.0/nbClasses]*nbClasses
    transitionProb = [[0]*nbClasses for _ in range(nbClasses)] 
    oProb = 0.25/(nbClasses-1.0)
    for i in range(nbClasses):
        for j in range(nbClasses):
            transitionProb[i][j]= 0.75 if i==j else oProb
            
    nbNodes=tree.nbNodes()
    lenTree=len(tree)
    data=tree.observation[lenTree-nbNodes:lenTree+1]
    classes=gmmEstimate(data,nbClasses)
    return MarkovModel(initialProb, transitionProb, classes)

def norm(vec):
    s=0
    for i in vec:
        s+=i**2
    return sqrt(s)

def reOrderLabels(markovModel, tree, rankingFunction=None):
    '''
    Reorder labels according to a ranking function.
    
    If no ranking function is provided
    '''
    estimLabel = tree.estimLabel
    if rankingFunction==None:
        rankingFunction=lambda x:norm(x.mean)
        
    classes = markovModel.classes
    nbClasses = markovModel.numLabels
    
    labels = [i for i in range(nbClasses)]
    ranks = [rankingFunction(classes[i]) for i in range(nbClasses)]
    labels=sorted(labels, key=lambda x:ranks[x])
    inverse = [0]*nbClasses
    for i in range(nbClasses):
        inverse[labels[i]]=i
    for n in tree.iteratorFromLeavesToRoot():
        estimLabel[n]=inverse[estimLabel[n]]
        
    nClasses=[]
    for i in range(nbClasses):
        nClasses.append(classes[labels[i]])
    markovModel.classes=nClasses
