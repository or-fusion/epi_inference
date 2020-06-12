import pytest
import numpy as np
import pandas as pd
import csv
from epi_inference.formulations.simulation import simulate_discrete_seiiir_deterministic
from epi_inference.util import roundall
from epi_inference.formulations.decay_lsq import run_decay_lsq
from epi_inference.formulations.decay_blike import run_decay_blike
from epi_inference.formulations.multinode_decay_lsq import run_multinode_decay_lsq
import matplotlib.pyplot as plt

def generate_datafile():
    """
    Test the decay inference using data from a simulation with the seiiir deterministic model
    """
    N = 1000000
    y0={'S': N, 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5
    gamma = 1/4
    rho = 1
    report_delay = 7
    tf = 25
    tx = [0]*tf
    tx[10] = 1

    dfdict = dict()
    dfpop = dict({'FIPS':[],'pop':[]})
    for beta in [0.25, 0.5, 0.75, 1.0, 1.25, 1.5]:
        Cdates,C,dates,T,S,E,I1,I2,I3,R = simulate_discrete_seiiir_deterministic(y0, tf, beta=beta,
                                                                        sigma=sigma, gamma=gamma,
                                                                        rho=rho, N=N,
                                                                        report_delay=report_delay,
                                                                        tx=tx)
        if len(dfdict) == 0:
            dfdict['Date'] = Cdates
        FIPSstr = '{:03d}'.format(int(100*beta))
        dfdict[FIPSstr] = C
        dfpop['FIPS'].append(FIPSstr)
        dfpop['pop'].append(N)

    dfC = pd.DataFrame(dfdict)
    dfC.to_csv('simulated_independent_county_different_beta.csv', index=False, quoting=csv.QUOTE_NONNUMERIC, date_format="%Y-%m-%d")

    dfC = dfC.round()
    dfC.to_csv('simulated_independent_county_different_beta_int.csv', index=False, quoting=csv.QUOTE_NONNUMERIC, date_format="%Y-%m-%d")

    dfP = pd.DataFrame(dfpop)
    dfP.to_csv('simulated_independent_county_different_beta_pop.csv', index=False, quoting=csv.QUOTE_NONNUMERIC)

if __name__ == '__main__':
    generate_datafile()
    

