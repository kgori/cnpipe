from hmmlearn import hmm
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

def model_builder(prob_change_state_01=1e-3, prob_change_state_10=1e-6,
                  mean_state0=0.0, mean_state1=1.0, stdev=0.3):
    model = hmm.GaussianHMM(n_components=2, covariance_type='diag',
                            init_params='s', params='smc',
                            transmat_prior=1e-6)
    model.transmat_ = np.array([
        [1 - prob_change_state_01, prob_change_state_01],
        [prob_change_state_10,     1 - prob_change_state_10]
    ])
    model.means_ = np.array([[mean_state0], [mean_state1]])
    model.covars_ = np.array([[stdev], [stdev]])
    return model
