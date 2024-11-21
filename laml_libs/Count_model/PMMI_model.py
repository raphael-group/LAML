from .PMM_base_model import PMM_base_model
from .IntensityData import IntensityData 
from .Param import Param
from math import *
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import numpy as np

DEFAULT_max_mu = 10
DEFAULT_max_nu = 1
DEFAULT_min_rho = 0.5

class PMMI_model(PMM_base_model):
    """
    This class represents the PMMI model as a generative process for the dynamic lineage tracing (DLT) data.
    The DLT data of this model must be an instance of the IntensityData class
    This class inherits all attributes and methods from PMM_base_model. 
    This class overrides logGamma of PMM_base_model
    """
    # PMMI = probabilistic mixed-type missing with noise 
    def __init__(self,treeList,data,prior,kw_params={}):
    # data is a dictionary; must have 'DLT_data' and data['DLT_data'] must be an instance of the IntensityData class
        mu = kw_params['mu'] if 'mu' in kw_params else 1 # arbitrarily chosen to set as a default value
        nu = kw_params['nu'] if 'nu' in kw_params else 0 # arbitrarily chosen to set as a default value
        phi = kw_params['phi'] if 'phi' in kw_params else 0.1 # arbitrarily chosen to set as a default value

        DEFAULT_min_eta = dict() 
        DEFAULT_max_eta = dict() 

        self.generative_emission_model = None
        # create the KDE model from training_data, assume it's a dictionary with label mapping to sampled feature vectors
        
        if 'training_data' in data:
            training_data = data['training_data']
            kde_dict = dict()
            for class_id in training_data:
                params = {'bandwidth': np.logspace(-1, 1, 20)}
                grid = GridSearchCV(KernelDensity(), params)
                grid.fit(training_data[class_id])
                # controls the smoothness of the density estimate
                kde_dict[class_id] = grid.best_estimator_
            self.generative_emission_model = kde_dict
        else:
            print("Asked to use PMMI model without providing training data...")

        params = Param(['mu','nu','phi'],[mu,nu,phi],[0,0,0],[DEFAULT_max_mu,DEFAULT_max_nu,1])
        super(PMMI_model,self).__init__(treeList,data,prior,params)
    
    def logGamma(self,k,x,c): 
        """
        Compute the emission probability in Layer 2
        This method overrides that of the Base_model class
            x is a cassette state of cassette k (data type: tuple of length J)
            c is the feature_vector corresponding to cassette k (data type: tuple of vectors, tuple of length J)
        Log-transforms the probabilities before returning.
        """    
        def score_sample(kde, samples): 
            ### From the kde documentation, this yields the log-likelihood of each sample under the model
            # this will be low for high-dimensional data
            log_density = kde.score_samples(samples)
            density = np.exp(log_density)
            return density[0]

        J = self.data['DLT_data'].J
        M = self.data['DLT_data'].alphabet.get_M(k)
        x_is_silenced = (x == tuple([-1]*J))
        phi = self.params.get_value('phi')
        rho = self.params.get_value('rho') # will be None if 'rho' not in self.params
        missing_state = tuple(['?']*J)   
        self.data['DLT_data']

        if x_is_silenced:
            p = 1 if c == missing_state else 0
        else:
            if c == missing_state:
                p = phi
            else:

                if self.generative_emission_model is not None:
                    p = 1

                    generating_state = x[0]
                    observation_vec = c 
                    score = score_sample(self.generative_emission_model[generating_state], [observation_vec])
                    p *= score
                else:
                    p = (1-phi)*rho if c==x else (1-phi)*(1-rho)/(M-2) 

        return log(p) if p>0 else None

