import pandas as pd
import os
from datetime import datetime
from epi_inference.formulations.simulation import simulate_discrete_seiiir_stochastic, simulate_continuous_seiiir_deterministic, simulate_discrete_seiiir_deterministic
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
"""
This module runs a deterministic model and compares the resutls with a number
of realizations of the stochastic simulation. There is an obvious discrepency
that is a result of the day discretization (instead of using a smaller discretization).
ToDo: rewrite the stochastic simulation to operate over smaller discretizations.
Note that these simulations are used only in testing, and not part of the results
produced by the package.
"""

def generate_simulation_figures(tf, fname):
    output_path = os.path.dirname(fname)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    N=100000
    y0={'S': N, 'E': 0, 'I1': 5, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5.2
    gamma = 1/4.3
    beta = 2.2*gamma
    rho = 10
    report_delay = 8
    tx = None
    #[0]*tf
    #tx[30] = 1

    Cdates,Ccases,dates,T,S,E,I1,I2,I3,R = simulate_discrete_seiiir_deterministic(y0, tf, beta=beta,
                                                                         sigma=sigma, gamma=gamma,
                                                                         rho=rho, N=N,
                                                                         report_delay=report_delay,
                                                                         tx=tx)

    Ccases = np.round(Ccases).astype(int)
    dfsim_S = pd.DataFrame({'dates': dates, 'values':S}).set_index('dates')
    dfsim_T = pd.DataFrame({'dates': dates, 'values':T}).set_index('dates')
    dfsim_E = pd.DataFrame({'dates': dates, 'values':E}).set_index('dates')
    dfsim_I1 = pd.DataFrame({'dates': dates, 'values':I1}).set_index('dates')
    dfsim_I2 = pd.DataFrame({'dates': dates, 'values':I2}).set_index('dates')
    dfsim_I3 = pd.DataFrame({'dates': dates, 'values':I3}).set_index('dates')
    dfsim_R = pd.DataFrame({'dates': dates, 'values':R}).set_index('dates')

    dfstoch_S = None
    dfstoch_T = None
    dfstoch_E = None
    dfstoch_I1 = None
    dfstoch_I2 = None
    dfstoch_I3 = None
    dfstoch_R = None
    for real in range(50):
#        Cdates, Ccases, dates, T, S, E, I1, I2, I3, R = simulate_continuous_seiiir_deterministic(y0, tf, beta=beta,
#                                                                                                 sigma=sigma, gamma=gamma,
#                                                                                                 rho=rho, N=N,
#                                                                                                 report_delay=report_delay,
#                                                                                                 tx=tx)
        Cdates,Ccases,dates,T,S,E,I1,I2,I3,R = simulate_discrete_seiiir_stochastic(y0, tf, beta=beta,
                                                                          sigma=sigma, gamma=gamma,
                                                                          rho=rho, N=N,
                                                                          report_delay=report_delay,
                                                                          tx=tx)
        #rdates, rcases, dates, T, S, E, I1, I2, I3, R = stochastic_reconstruction(Cdates, Ccases, N)

        if dfstoch_S is None:
            dfstoch_S = pd.DataFrame({'dates': pd.to_datetime(dates)}).set_index('dates')
            dfstoch_T = pd.DataFrame({'dates': pd.to_datetime(dates)}).set_index('dates')
            dfstoch_E = pd.DataFrame({'dates': pd.to_datetime(dates)}).set_index('dates')
            dfstoch_I1 = pd.DataFrame({'dates': pd.to_datetime(dates)}).set_index('dates')
            dfstoch_I2 = pd.DataFrame({'dates': pd.to_datetime(dates)}).set_index('dates')
            dfstoch_I3 = pd.DataFrame({'dates': pd.to_datetime(dates)}).set_index('dates')
            dfstoch_R = pd.DataFrame({'dates': pd.to_datetime(dates)}).set_index('dates')
            
        dfstoch_S['{}'.format(real)] = S
        dfstoch_T['{}'.format(real)] = T
        dfstoch_E['{}'.format(real)] = E
        dfstoch_I1['{}'.format(real)] = I1
        dfstoch_I2['{}'.format(real)] = I2
        dfstoch_I3['{}'.format(real)] = I3
        dfstoch_R['{}'.format(real)] = R

    with PdfPages(fname) as pdf:
        ax = dfstoch_S.plot(color='silver', legend=False)
        dfsim_S[dfsim_S.index.isin(dfstoch_S.index)].plot(ax=ax, color='black', legend='Simulated S')
        plt.title('S comparison')
        pdf.savefig()
        plt.close()

        ax = dfstoch_T.plot(color='silver', legend=False)
        dfsim_T[dfsim_T.index.isin(dfstoch_T.index)].plot(ax=ax, color='black', legend='Simulated T')
        plt.title('Comparison of daily transmissions')
        pdf.savefig()
        plt.close()
        
        ax = dfstoch_E.plot(color='silver', legend=False)
        dfsim_E[dfsim_E.index.isin(dfstoch_E.index)].plot(ax=ax, color='black', legend='Simulated E')
        plt.title('E comparison')
        pdf.savefig()
        plt.close()

        ax = dfstoch_I1.plot(color='silver', legend=False)
        dfsim_I1[dfsim_I1.index.isin(dfstoch_I1.index)].plot(ax=ax, color='black', legend='Simulated I1')
        plt.title('I1 comparison')
        pdf.savefig()
        plt.close()
        
        ax = dfstoch_I2.plot(color='silver', legend=False)
        dfsim_I2[dfsim_I2.index.isin(dfstoch_I2.index)].plot(ax=ax, color='black', legend='Simulated I2')
        plt.title('I2 comparison')
        pdf.savefig()
        plt.close()
        
        ax = dfstoch_I3.plot(color='silver', legend=False)
        dfsim_I3[dfsim_I3.index.isin(dfstoch_I3.index)].plot(ax=ax, color='black', legend='Simulated I3')
        plt.title('I3 comparison')
        pdf.savefig()
        plt.close()

        ax = dfstoch_R.plot(color='silver', legend=False)
        dfsim_R[dfsim_R.index.isin(dfstoch_R.index)].plot(ax=ax, color='black', legend='Simulated R')
        plt.title('R comparison')
        pdf.savefig()
        plt.close()

        #ax = dfstoch_T.cumsum().plot(color='silver', legend=False)
        #dfsim_R.plot(ax=ax, color='black', legend='Simulated R')
        #plt.title('R comparison')
        #pdf.savefig()
        #plt.close()

    return

if __name__ == '__main__':
    #np.random.seed(1975)
    generate_simulation_figures(365, './figures/simulations.pdf')
    

        
