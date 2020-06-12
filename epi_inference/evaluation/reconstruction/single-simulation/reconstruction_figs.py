import pandas as pd
import os
from datetime import datetime
from epi_inference.reconstruction.common import reported_cases_from_cumulative
from epi_inference.reconstruction.stochastic import stochastic_reconstruction
from epi_inference.simulation.simulation import simulate_discrete_seiiir_stochastic
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
"""
This module runs reconstructions on data from a stochastic simulation
and produces some figures showing the results.
"""

def compare_simulation_and_reconstruction(tf, fname):
    output_path = os.path.dirname(fname)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    N=100000
    y0={'S': N, 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5.2
    gamma = 1/4.3
    beta = 2.2*gamma
    rho = 10
    report_delay = 8
    tx = [0]*tf
    tx[30] = 1

    Cdates,Ccases,dates,T,S,E,I1,I2,I3,R = \
        simulate_discrete_seiiir_stochastic(y0, tf, beta=beta,
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
    for real in range(100):
        rdates, reported_cases_per_day = reported_cases_from_cumulative(Cdates, Ccases)
        dates, T, S, E, I1, I2, I3, R = stochastic_reconstruction(rdates, reported_cases_per_day, N, 1)

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

def generate_reconstruction_figures(Cdates, Ccases, N, fname, comment):
    output_path = os.path.dirname(fname)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    dfstoch_S = None
    dfstoch_T = None
    dfstoch_E = None
    dfstoch_I1 = None
    dfstoch_I2 = None
    dfstoch_I3 = None
    dfstoch_R = None
    for real in range(100):
        rdates, reported_cases_per_day = reported_cases_from_cumulative(Cdates, Ccases)
        dates, T, S, E, I1, I2, I3, R = stochastic_reconstruction(rdates, reported_cases_per_day, N, 1)

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
        plt.title('S {}'.format(comment))
        pdf.savefig()
        plt.close()

        ax = dfstoch_T.plot(color='silver', legend=False)
        plt.title('daily transmissions {}'.format(comment))
        pdf.savefig()
        plt.close()
        
        ax = dfstoch_E.plot(color='silver', legend=False)
        plt.title('E {}'.format(comment))
        pdf.savefig()
        plt.close()

        ax = dfstoch_I1.plot(color='silver', legend=False)
        plt.title('I1 {}'.format(comment))
        pdf.savefig()
        plt.close()
        
        ax = dfstoch_I2.plot(color='silver', legend=False)
        plt.title('I2 {}'.format(comment))
        pdf.savefig()
        plt.close()
        
        ax = dfstoch_I3.plot(color='silver', legend=False)
        plt.title('I3 {}'.format(comment))
        pdf.savefig()
        plt.close()

        ax = dfstoch_R.plot(color='silver', legend=False)
        plt.title('R {}'.format(comment))
        pdf.savefig()
        plt.close()


if __name__ == '__main__':
    np.random.seed(1975)
    compare_simulation_and_reconstruction(120*3, './figures/reconstruction-comparison-long-time-horizon.pdf')
    np.random.seed(1975)
    compare_simulation_and_reconstruction(120, './figures/reconstruction-comparison-med-time-horizon.pdf')
    np.random.seed(1975)
    compare_simulation_and_reconstruction(60, './figures/reconstruction-comparison-short-time-horizon.pdf')
    
    ##
    # what do results look like when the case counts are very low
    ##
    Cdates = pd.date_range(end=datetime(year=2020, month=5, day=15), periods=30).to_pydatetime().tolist()
    Ccases = [0]*len(Cdates)
    for i in range(10,len(Cdates)):
        Ccases[i]=1
    generate_reconstruction_figures(Cdates, Ccases, 100000, './figures/reconstruction-low-cases.pdf', '(one reported case on {})'.format(Cdates[10]))
    

        
