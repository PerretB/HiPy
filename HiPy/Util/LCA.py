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
Optimized lowest common ancestor

Although not the best one, this one performs a O(n.log(n)) pre-processing and has a O(log(n)) query time.
'''



from math import log, ceil


class LCA(object):

    def __init__(self, tree):
        self.tree = tree
        nbNodes = len(tree)
        P = [None]*nbNodes
        logNbNodes = ceil(log(nbNodes))
        for i in range(len(P)):
            P[i] = [-1]*logNbNodes
        self.logNbNodes = logNbNodes
        self.P = P


    def __preprocess(self):
        P = self.P
        for i in range(len(P)):
            P[i][0] = self.tree[i]

        for j in range(1, self.logNbNodes):
            for i in range(len(P)):
                if P[i][j-1] != -1:
                    P[i][j] = P[P[i][j-1]][j-1]



__author__ = 'perretb'




