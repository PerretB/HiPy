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

from HiPy.Hierarchies.WatershedHierarchy import constructAltitudeBPT, \
    reweightMSTByAttribute, transformAltitudeBPTtoWatershedHierarchy
from HiPy.Processing.Attributes import addAttributeChildren, addAttributeArea


def constructExactRandomSeedsWatershed(adjacency, computeAttributeFunction=addAttributeArea,
                                       seedMeasureAttribute="area", verbose=False):
    if verbose:
        print("Constructing Watershed hierarchy...")
    bpt = transformAltitudeBPTtoWatershedHierarchy(constructAltitudeBPT(adjacency))
    if verbose:
        print("Computing attribute...")
    computeAttributeFunction(bpt)
    if verbose:
        print("Computing contour pdf...")
    addAttributeRandomSeedBoundaryProbability(bpt, seedMeasureAttribute=seedMeasureAttribute, numberOfSeeds=2,
                                              attributeName="randomSeedPdf")
    newadj = reweightMSTByAttribute(bpt, "randomSeedPdf", extinctionValue=False)
    if verbose:
        print("Re-Constructing BPT...")
    return constructAltitudeBPT(newadj)


def addAttributeRandomSeedBoundaryProbability(bpt, seedMeasureAttribute="area", numberOfSeeds=2,
                                              attributeName="randomSeedPdf"):
    attr, created = bpt.addAttribute(attributeName, True)
    if not created:
        return attr
    addAttributeChildren(bpt)
    children = bpt.children
    measure = bpt.getAttribute(seedMeasureAttribute)
    measure0 = measure[-1]
    for i in bpt.iteratorFromLeavesToRoot(includeLeaves=False):
        attr[i] = 1 - (1 - measure[children[i][0]] / measure0) ** numberOfSeeds - (1 - measure[
            children[i][1]] / measure0) ** numberOfSeeds + (1 - measure[i] / measure0) ** numberOfSeeds
    return attr