from .PMM_base import PMM_model
from math import *
from laml_libs import DEFAULT_max_mu, DEFAULT_max_nu
from .utils import *
from .AlleleTable import AlleleTable
from .CharMtrx import CharMtrx

class PMM_charMtrx(PMM_model):
    def __init__(self,treeList,data,prior,**kw_params):
    # data is a dictionary mapping a data module's name to values
    # data must have 'charMtrx'; data['charMtrx'] must be an instance of class CharMtrx
    # kw_params: a dictionary mapping param's name to value. Here it should inclue 'mu', 'nu', 'phi'
        Q = prior['Q']
        charMtrx = data['charMtrx'] # this is an instance of class CharMtrx
        K = charMtrx.K
        J = charMtrx.J
        #alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])
        allele_table = charMtrx_2_alleleTable(charMtrx)#,alphabet)
        data['alleleTable'] = allele_table

        super(PMM_charMtrx,self).__init__(treeList,data,prior,**kw_params)
