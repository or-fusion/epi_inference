import pandas as pd
import os
from epi_inference.reconstruction.common import reported_cases_from_cumulative
from epi_inference.reconstruction.stochastic import stochastic_reconstruction
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
"""
This module runs reconstructions on data from a stochastic simulation
and produces some figures showing the results.
"""

def compare_florida_reconstruction(seiiir_fname, reported_cases_fname, geodata_fname, output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # read the true data from the SEIIIR simulation
    seirdf = pd.read_csv(seiiir_fname, parse_dates=['Date'])
    seirdf['dates'] = pd.to_datetime(seirdf['Date'])
    seirdf = seirdf.set_index('dates')

    # read the cumulative reported cases
    crcdf = pd.read_csv(reported_cases_fname, parse_dates=['Date'])
    crcdf['dates'] = pd.to_datetime(crcdf['Date'])
    cdcdf = crcdf.set_index('dates')

    # read the populations
    popdf = pd.read_csv(geodata_fname)
    popdf = popdf.set_index('geoid')
    populations = popdf['pop2010'].to_dict()
    populations = {str(int(k)):v for k,v in populations.items()}
    
    # get the list of counties
    counties = set(seirdf.columns.to_list())
    counties.remove('Date')
    counties.remove('comp')
    counties = sorted(counties)

    # loop through all the counties and perform the reconstruction
    # based on the reported cases
    pdf = PdfPages(os.path.join(output_path, 'reconstruction-comparison-florida.pdf'))
    for c in counties:
        print('...', c)
        Cdates = crcdf['dates'].tolist()
        Ccases = crcdf[c].astype(int).tolist()
        dfsim_S = seirdf[seirdf['comp'] == 'S'][c]
        dfsim_E = seirdf[seirdf['comp'] == 'E'][c]
        dfsim_I1 = seirdf[seirdf['comp'] == 'I1'][c]
        dfsim_I2 = seirdf[seirdf['comp'] == 'I2'][c]
        dfsim_I3 = seirdf[seirdf['comp'] == 'I3'][c]
        dfsim_R = seirdf[seirdf['comp'] == 'R'][c]

        dfrecon_S = None
        dfrecon_T = None
        dfrecon_E = None
        dfrecon_I1 = None
        dfrecon_I2 = None
        dfrecon_I3 = None
        dfrecon_R = None

        for real in range(100):
            bunch_reported_cases_per_day = reported_cases_from_cumulative(dates=Cdates,
                                                                            cumulative_reported_cases=Ccases)
            
            bunch_recon = stochastic_reconstruction(
                dates=bunch_reported_cases_per_day.dates,
                reported_cases_per_day=bunch_reported_cases_per_day.values,
                population=populations[c],
                n_steps_per_day=4,
                reporting_delay_mean=8,
                reporting_delay_dev=1.35,
                fixed_incubation=5.2,
                infectious_lower=2.6,
                infectious_upper=6.0,
                )

            if dfrecon_S is None:
                dfrecon_S = pd.DataFrame({'dates': pd.to_datetime(bunch_recon.dates)}).set_index('dates')
                dfrecon_T = pd.DataFrame({'dates': pd.to_datetime(bunch_recon.dates)}).set_index('dates')
                dfrecon_E = pd.DataFrame({'dates': pd.to_datetime(bunch_recon.dates)}).set_index('dates')
                dfrecon_I1 = pd.DataFrame({'dates': pd.to_datetime(bunch_recon.dates)}).set_index('dates')
                dfrecon_I2 = pd.DataFrame({'dates': pd.to_datetime(bunch_recon.dates)}).set_index('dates')
                dfrecon_I3 = pd.DataFrame({'dates': pd.to_datetime(bunch_recon.dates)}).set_index('dates')
                dfrecon_R = pd.DataFrame({'dates': pd.to_datetime(bunch_recon.dates)}).set_index('dates')

            dfrecon_S['{}'.format(real)] = bunch_recon.S
            dfrecon_T['{}'.format(real)] = bunch_recon.transmissions
            dfrecon_E['{}'.format(real)] = bunch_recon.E
            dfrecon_I1['{}'.format(real)] = bunch_recon.I1
            dfrecon_I2['{}'.format(real)] = bunch_recon.I2
            dfrecon_I3['{}'.format(real)] = bunch_recon.I3
            dfrecon_R['{}'.format(real)] = bunch_recon.R

        ax = dfrecon_S.plot(color='silver', legend=False)
        dfsim_S[dfsim_S.index.isin(dfrecon_S.index)].plot(ax=ax, color='black', legend='Simulated S')
        plt.title('S comparison')
        pdf.savefig()
        plt.close()

        ax = dfrecon_E.plot(color='silver', legend=False)
        dfsim_E[dfsim_E.index.isin(dfrecon_E.index)].plot(ax=ax, color='black', legend='Simulated E')
        plt.title('E comparison')
        pdf.savefig()
        plt.close()

        ax = dfrecon_I1.plot(color='silver', legend=False)
        dfsim_I1[dfsim_I1.index.isin(dfrecon_I1.index)].plot(ax=ax, color='black', legend='Simulated I1')
        plt.title('I1 comparison')
        pdf.savefig()
        plt.close()

        ax = dfrecon_I2.plot(color='silver', legend=False)
        dfsim_I2[dfsim_I2.index.isin(dfrecon_I2.index)].plot(ax=ax, color='black', legend='Simulated I2')
        plt.title('I2 comparison')
        pdf.savefig()
        plt.close()

        ax = dfrecon_I3.plot(color='silver', legend=False)
        dfsim_I3[dfsim_I3.index.isin(dfrecon_I3.index)].plot(ax=ax, color='black', legend='Simulated I3')
        plt.title('I3 comparison')
        pdf.savefig()
        plt.close()

        ax = dfrecon_R.plot(color='silver', legend=False)
        dfsim_R[dfsim_R.index.isin(dfrecon_R.index)].plot(ax=ax, color='black', legend='Simulated R')
        lower_percentile = dfrecon_R.quantile(0.025, axis=1)[dfrecon_R.index[-1]]
        upper_percentile = dfrecon_R.quantile(0.975, axis=1)[dfrecon_R.index[-1]]
        sim_value = dfsim_R[dfrecon_R.index[-1]]
        msg = ''
        if sim_value < lower_percentile or sim_value > upper_percentile:
            msg = '*'
            print('Simulated R outside of 95th percentiles for count {}: ({} ({}) {}) {}'.format(c, lower_percentile, sim_value, upper_percentile, msg))

        plt.title('R comparison ({} ({}) {}) {}'.format(lower_percentile, sim_value, upper_percentile, msg))
        pdf.savefig()
        plt.close()

    pdf.close()



if __name__ == '__main__':
    np.random.seed(1975)
    seiiir_fname = './data/FL_SEIIIR_R0_2.25_short_realization_4.csv'
    reported_cases_fname = './data/FL_cumulative_reported_cases_R0_2.25_short_realization_4.csv'
    geodata_fname = './data/geodata.csv'
    output_path = './figures/'
    
    compare_florida_reconstruction(seiiir_fname, reported_cases_fname, geodata_fname, output_path)
    
