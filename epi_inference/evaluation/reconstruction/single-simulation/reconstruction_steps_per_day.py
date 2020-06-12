import pandas as pd
import os
from datetime import datetime
from epi_inference.reconstruction.common import reported_cases_from_cumulative
from epi_inference.reconstruction.stochastic import stochastic_reconstruction
from epi_inference.reconstruction.deterministic import reconstruct_states_deterministic_decay
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

    sim = simulate_discrete_seiiir_stochastic(y0, tf, beta=beta,
                                            sigma=sigma, gamma=gamma,
                                            rho=rho, N=N,
                                            report_delay=report_delay,
                                            tx=tx)

    dfsim_S = pd.DataFrame({'dates': sim.SEIIIR.dates, 'values':sim.SEIIIR.S}).set_index('dates')
    dfsim_T = pd.DataFrame({'dates': sim.SEIIIR.dates, 'values':sim.SEIIIR.transmissions}).set_index('dates')
    dfsim_E = pd.DataFrame({'dates': sim.SEIIIR.dates, 'values':sim.SEIIIR.E}).set_index('dates')
    dfsim_I1 = pd.DataFrame({'dates': sim.SEIIIR.dates, 'values':sim.SEIIIR.I1}).set_index('dates')
    dfsim_I2 = pd.DataFrame({'dates': sim.SEIIIR.dates, 'values':sim.SEIIIR.I2}).set_index('dates')
    dfsim_I3 = pd.DataFrame({'dates': sim.SEIIIR.dates, 'values':sim.SEIIIR.I3}).set_index('dates')
    dfsim_R = pd.DataFrame({'dates': sim.SEIIIR.dates, 'values':sim.SEIIIR.R}).set_index('dates')

    Ccases = np.round(sim.cumulative_reported_cases.values).astype(int)
    Cdates = sim.cumulative_reported_cases.dates

    dfstoch_S = None
    dfstoch_T = None
    dfstoch_E = None
    dfstoch_I1 = None
    dfstoch_I2 = None
    dfstoch_I3 = None
    dfstoch_R = None
    bunch_reported_cases_per_day = reported_cases_from_cumulative(dates=Cdates,
                                                                      cumulative_reported_cases=Ccases)
    for real in range(100):
        bunch_recon = stochastic_reconstruction(
            dates=bunch_reported_cases_per_day.dates,
            reported_cases_per_day=bunch_reported_cases_per_day.values,
            population=N,
            n_steps_per_day=1,
            reporting_delay_mean=8,
            reporting_delay_dev=1.35,
            reporting_multiplier=10,
            fixed_incubation=5.2,
            infectious_lower=2.6,
            infectious_upper=6.0,
        )

        if dfstoch_S is None:
            dates = bunch_recon.dates
            dfstoch_S = pd.DataFrame({'dates': pd.to_datetime(dates)}).set_index('dates')
            dfstoch_T = pd.DataFrame({'dates': pd.to_datetime(dates)}).set_index('dates')
            dfstoch_E = pd.DataFrame({'dates': pd.to_datetime(dates)}).set_index('dates')
            dfstoch_I1 = pd.DataFrame({'dates': pd.to_datetime(dates)}).set_index('dates')
            dfstoch_I2 = pd.DataFrame({'dates': pd.to_datetime(dates)}).set_index('dates')
            dfstoch_I3 = pd.DataFrame({'dates': pd.to_datetime(dates)}).set_index('dates')
            dfstoch_R = pd.DataFrame({'dates': pd.to_datetime(dates)}).set_index('dates')
            
        dfstoch_S['{}'.format(real)] = bunch_recon.S
        dfstoch_T['{}'.format(real)] = bunch_recon.transmissions
        dfstoch_E['{}'.format(real)] = bunch_recon.E
        dfstoch_I1['{}'.format(real)] = bunch_recon.I1
        dfstoch_I2['{}'.format(real)] = bunch_recon.I2
        dfstoch_I3['{}'.format(real)] = bunch_recon.I3
        dfstoch_R['{}'.format(real)] = bunch_recon.R

    # do a deterministic reconstruction
    det_recon = reconstruct_states_deterministic_decay(
        dates=bunch_reported_cases_per_day.dates,
        reported_cases_per_day=bunch_reported_cases_per_day.values,
        population=N,
        sigma=1/5.2,
        gamma=1/4.3,
        reporting_factor=10,
        report_delay=8
        )

    with PdfPages(fname) as pdf:
        ax = dfstoch_S.plot(color='silver', legend=False)
        dfsim_S[dfsim_S.index.isin(dfstoch_S.index)].plot(ax=ax, color='black', legend='Simulated S')
        detdf = pd.DataFrame({'dates': pd.to_datetime(det_recon.dates), 'det_recon_S':det_recon.S}).set_index('dates')
        detdf[detdf.index.isin(dfstoch_S.index)].plot(ax=ax, color='red', legend='det_recon')
        plt.title('S comparison')
        pdf.savefig()
        plt.close()

        ax = dfstoch_T.plot(color='silver', legend=False)
        dfsim_T[dfsim_T.index.isin(dfstoch_T.index)].plot(ax=ax, color='black', legend='Simulated T')
        detdf = pd.DataFrame({'dates': pd.to_datetime(det_recon.dates), 'det_recon_T':det_recon.transmissions}).set_index('dates')
        detdf[detdf.index.isin(dfstoch_S.index)].plot(ax=ax, color='red', legend='det_recon')
        plt.title('Comparison of daily transmissions')
        pdf.savefig()
        plt.close()
        
        ax = dfstoch_E.plot(color='silver', legend=False)
        dfsim_E[dfsim_E.index.isin(dfstoch_E.index)].plot(ax=ax, color='black', legend='Simulated E')
        detdf = pd.DataFrame({'dates': pd.to_datetime(det_recon.dates), 'det_recon_E':det_recon.E}).set_index('dates')
        detdf[detdf.index.isin(dfstoch_S.index)].plot(ax=ax, color='red', legend='det_recon')
        plt.title('E comparison')
        pdf.savefig()
        plt.close()

        ax = dfstoch_I1.plot(color='silver', legend=False)
        dfsim_I1[dfsim_I1.index.isin(dfstoch_I1.index)].plot(ax=ax, color='black', legend='Simulated I1')
        detdf = pd.DataFrame({'dates': pd.to_datetime(det_recon.dates), 'det_recon_I1':det_recon.I1}).set_index('dates')
        detdf[detdf.index.isin(dfstoch_S.index)].plot(ax=ax, color='red', legend='det_recon')
        plt.title('I1 comparison')
        pdf.savefig()
        plt.close()
        
        ax = dfstoch_I2.plot(color='silver', legend=False)
        dfsim_I2[dfsim_I2.index.isin(dfstoch_I2.index)].plot(ax=ax, color='black', legend='Simulated I2')
        detdf = pd.DataFrame({'dates': pd.to_datetime(det_recon.dates), 'det_recon_I2':det_recon.I2}).set_index('dates')
        detdf[detdf.index.isin(dfstoch_S.index)].plot(ax=ax, color='red', legend='det_recon')
        plt.title('I2 comparison')
        pdf.savefig()
        plt.close()
        
        ax = dfstoch_I3.plot(color='silver', legend=False)
        dfsim_I3[dfsim_I3.index.isin(dfstoch_I3.index)].plot(ax=ax, color='black', legend='Simulated I3')
        detdf = pd.DataFrame({'dates': pd.to_datetime(det_recon.dates), 'det_recon_I3':det_recon.I3}).set_index('dates')
        detdf[detdf.index.isin(dfstoch_S.index)].plot(ax=ax, color='red', legend='det_recon')
        plt.title('I3 comparison')
        pdf.savefig()
        plt.close()

        ax = dfstoch_R.plot(color='silver', legend=False)
        dfsim_R[dfsim_R.index.isin(dfstoch_R.index)].plot(ax=ax, color='black', legend='Simulated I3')
        detdf = pd.DataFrame({'dates': pd.to_datetime(det_recon.dates), 'det_recon_R':det_recon.R}).set_index('dates')
        detdf[detdf.index.isin(dfstoch_S.index)].plot(ax=ax, color='red', legend='det_recon')
        plt.title('R comparison')
        pdf.savefig()
        plt.close()


        dfsim_R = dfsim_R[:dfstoch_R.index[-1]].astype(float)
        dfstoch_R = dfstoch_R[dfsim_R.index[0]:]
        dferr_R = dfstoch_R.subtract(dfsim_R['values'], axis=0)
        ax = dferr_R.mean(axis=1).plot(color='silver', legend=False)
        print(dferr_R.mean(axis=1))
        plt.title('R Errors')
        pdf.savefig()
        plt.close()

        #ax = dfstoch_T.cumsum().plot(color='silver', legend=False)
        #dfsim_R.plot(ax=ax, color='black', legend='Simulated R')
        #plt.title('R comparison')
        #pdf.savefig()
        #plt.close()

    return

def generate_reconstruction_figures(Cdates, Ccases, fname, comment):
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

"""
def compare_stepsize_and_not_stepsize_reconstruction(tf):
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

    np.random.seed(42)
    rdates, rcases, dates, T, S, E, I1, I2, I3, R = stochastic_reconstruction(Cdates, Ccases, N, 1)
    np.random.seed(42)
    rdatesb, rcasesb, datesb, Tb, Sb, Eb, I1b, I2b, I3b, Rb = stochastic_reconstruction(Cdates, Ccases, N)

    for i,v in enumerate(rcases):
        assert rdates[i] == rdatesb[i]
        assert abs(v - rcasesb[i]) <= 1e-10

    for i in range(len(dates)):
        assert dates[i] == datesb[i]
        assert abs(T[i] - Tb[i]) <= 1e-10
        assert abs(S[i] - Sb[i]) <= 1e-10
        assert abs(E[i] - Eb[i]) <= 1e-10
        assert abs(I1[i] - I1b[i]) <= 1e-10
        assert abs(I2[i] - I2b[i]) <= 1e-10
        assert abs(I3[i] - I3b[i]) <= 1e-10
        assert abs(R[i] - Rb[i]) <= 1e-10
"""

if __name__ == '__main__':
    np.random.seed(1975)
    # compare_stepsize_and_not_stepsize_reconstruction(tf=60)
    compare_simulation_and_reconstruction(120, './figures/reconstruction-steps_per_day.pdf')
    

        
