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

"""
Created on 6 juil. 2015

@author: perretb
"""

import copy

from HiPy.Util import VMath
from HiPy.Util.VMath import medianV


class AbstractAccumulator(object):
    def reset(self):
        raise NotImplementedError("Unsupported method" + " reset")

    def accumulate(self, *args):
        raise NotImplementedError("Unsupported method" + " accumulate")

    def result(self):
        raise NotImplementedError("Unsupported method" + " result")


class BasicAccumulator(AbstractAccumulator):
    def __init__(self, accumulateFunction, resultFunction, resetValues):
        self.resetValues = resetValues
        self.accumulateFunction = accumulateFunction
        self.resultFunction = resultFunction
        self.values = None

    def reset(self):
        self.values = copy.copy(self.resetValues)

    def accumulate(self, *args):
        self.accumulateFunction(self.values, *args)

    def result(self):
        return self.resultFunction(self.values)

    def copy(self):
        return BasicAccumulator(self.accumulateFunction, self.resultFunction, self.resetValues)


    @staticmethod
    def getCounterAccumulator():
        def accFun(values, newValue, *_):
            values[0] += 1
        return BasicAccumulator(accFun, lambda values: values[0], [0])

    @staticmethod
    def getSumAccumulator():
        def accFun(values, newValue, *_):
            values[0] = values[0] + newValue

        return BasicAccumulator(accFun, lambda values: values[0], [0])

    @staticmethod
    def getWeightedSumAccumulator():
        def accFun(values, newValue, weight, *_):
            values[0] = values[0] + weight * newValue

        return BasicAccumulator(accFun, lambda values: values[0], [0])

    @staticmethod
    def getMinAccumulator(initValue=99999999):
        def accFun(values, newValue, *_):
            values[0] = min(values[0], newValue)

        return BasicAccumulator(accFun, lambda values: values[0], [initValue])

    @staticmethod
    def getMaxAccumulator(initValue=-99999999):
        def accFun(values, newValue, *_):
            values[0] = max(values[0], newValue)

        return BasicAccumulator(accFun, lambda values: values[0], [initValue])

    @staticmethod
    def getMeanAccumulator():
        def accFun(values, newValue, *_):
            values[0] += newValue
            values[1] += 1

        return BasicAccumulator(accFun, lambda values: values[0] / values[1], [0, 0])

    @staticmethod
    def getMeanVAccumulator(dim):
        def accFun(values, newValue, *_):
            values[0] = VMath.addV(values[0], newValue)
            values[1] += 1

        return BasicAccumulator(accFun, lambda values: VMath.divS(values[0], values[1]), [[0]*dim, 0])


    @staticmethod
    def getMedianAccumulator():
        def accFun(values, newValue, *_):
            values.append(newValue)

        return BasicAccumulator(accFun, lambda values: medianV(values), [])

    @staticmethod
    def getWeightedMeanAccumulator(normalize=False):
        def accFun(values, newValue, weight, *_):
            values[0] += newValue * weight
            values[1] += weight

        return BasicAccumulator(accFun, lambda values: values[0] / values[1] if normalize else values[0], [0, 0])
