import pandas as pd
pd.set_option("display.max_rows", None)
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

def compare_florida_reconstruction(seiiir_fname, det_recon_fname, recon_folder_name, geodata_fname, output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # read the true data from the SEIIIR simulation
    seirdf = pd.read_csv(seiiir_fname, parse_dates=['Date'])
    seirdf['Date'] = pd.to_datetime(seirdf['Date'])
    seirdf = seirdf.set_index('Date')

    # read the deterministic reconstruction
    detrecondf = pd.read_csv(det_recon_fname, parse_dates=['time'])
    detrecondf['Date'] = pd.to_datetime(detrecondf['time'])
    detrecondf = detrecondf.set_index('Date')

    # read the populations
    popdf = pd.read_csv(geodata_fname)
    popdf = popdf.set_index('geoid')
    populations = popdf['pop2010'].to_dict()
    populations = {str(int(k)):v for k,v in populations.items()}
    
    # get the list of counties
    counties = set(seirdf.columns.to_list())
    if 'Date' in counties:
        counties.remove('Date')
    counties.remove('comp')
    counties = sorted(counties)

    # get a list of all the files in the recon folder
    recon_files = list()
    for f in os.listdir(recon_folder_name):
        fname = os.path.join(recon_folder_name, f)
        if os.path.isfile(fname):
            filenameonly, extension = os.path.splitext(f)
            if extension == '.csv':
                recon_files.append(fname)

    # loop through all the counties and perform the reconstruction
    # based on the reported cases
    pdf = PdfPages(os.path.join(output_path, 'reconstruction-comparison-florida.pdf'))
    for c in counties:
        # this is expensive - we are opening the files for each county
        print('...', c)
        dfsim_S = pd.DataFrame(seirdf[seirdf['comp'] == 'S'][c])
        dfsim_E = seirdf[seirdf['comp'] == 'E'][c]
        dfsim_I1 = seirdf[seirdf['comp'] == 'I1'][c]
        dfsim_I2 = seirdf[seirdf['comp'] == 'I2'][c]
        dfsim_I3 = seirdf[seirdf['comp'] == 'I3'][c]
        dfsim_R = seirdf[seirdf['comp'] == 'R'][c]

        dfdetrecon_S = detrecondf[detrecondf['comp']=='S'][c]
        dfdetrecon_E = detrecondf[detrecondf['comp']=='E'][c]
        dfdetrecon_I1 = detrecondf[detrecondf['comp']=='I1'][c]
        dfdetrecon_I2 = detrecondf[detrecondf['comp']=='I2'][c]
        dfdetrecon_I3 = detrecondf[detrecondf['comp']=='I3'][c]
        dfdetrecon_R = detrecondf[detrecondf['comp']=='R'][c]

        dfrecon_S = None
        dfrecon_E = None
        dfrecon_I1 = None
        dfrecon_I2 = None
        dfrecon_I3 = None
        dfrecon_R = None

        idx = 0
        for rfname in recon_files:
            recondf = pd.read_csv(rfname)
            recondf = recondf.rename({'time':'Date'}, axis='columns')
            recondf.set_index('Date')
            
            tempdf = recondf[recondf['comp']=='S'].set_index('Date')
            if dfrecon_S is None:
                dfrecon_S = pd.DataFrame({'Date': tempdf.index}).set_index('Date')
                dfrecon_E = pd.DataFrame({'Date': tempdf.index}).set_index('Date')
                dfrecon_I1 = pd.DataFrame({'Date': tempdf.index}).set_index('Date')
                dfrecon_I2 = pd.DataFrame({'Date': tempdf.index}).set_index('Date')
                dfrecon_I3 = pd.DataFrame({'Date': tempdf.index}).set_index('Date')
                dfrecon_R = pd.DataFrame({'Date': tempdf.index}).set_index('Date')

            tempdf = recondf[recondf['comp']=='S'].set_index('Date')
            dfrecon_S['realization_{}'.format(idx)] = tempdf[c]
            tempdf = recondf[recondf['comp']=='E'].set_index('Date')
            dfrecon_E['realization_{}'.format(idx)] = tempdf[c]
            tempdf = recondf[recondf['comp']=='I1'].set_index('Date')
            dfrecon_I1['realization_{}'.format(idx)] = tempdf[c]
            tempdf = recondf[recondf['comp']=='I2'].set_index('Date')
            dfrecon_I2['realization_{}'.format(idx)] = tempdf[c]
            tempdf = recondf[recondf['comp']=='I3'].set_index('Date')
            dfrecon_I3['realization_{}'.format(idx)] = tempdf[c]
            tempdf = recondf[recondf['comp']=='R'].set_index('Date')
            dfrecon_R['realization_{}'.format(idx)] = tempdf[c]
            idx += 1
            
        dfrecon_S.index = pd.to_datetime(dfrecon_S.index)
        ax = dfrecon_S.plot(color='silver', legend=False)
        dfsim_S[dfsim_S.index.isin(dfrecon_S.index)].plot(ax=ax, color='black', legend='Simulated S')
        dfdetrecon_S[dfdetrecon_S.index.isin(dfrecon_S.index)].plot(ax=ax, color='red', legend='Deterministic_Recon S')
        plt.title('S comparison')
        pdf.savefig()
        plt.close()

        dfrecon_E.index = pd.to_datetime(dfrecon_S.index)
        ax = dfrecon_E.plot(color='silver', legend=False)
        dfsim_E[dfsim_E.index.isin(dfrecon_E.index)].plot(ax=ax, color='black', legend='Simulated E')
        dfdetrecon_E[dfdetrecon_E.index.isin(dfrecon_E.index)].plot(ax=ax, color='red', legend='Deterministic_Recon E')
        plt.title('E comparison')
        pdf.savefig()
        plt.close()

        dfrecon_I1.index = pd.to_datetime(dfrecon_S.index)
        ax = dfrecon_I1.plot(color='silver', legend=False)
        dfsim_I1[dfsim_I1.index.isin(dfrecon_I1.index)].plot(ax=ax, color='black', legend='Simulated I1')
        dfdetrecon_I1[dfdetrecon_I1.index.isin(dfrecon_I1.index)].plot(ax=ax, color='red', legend='Deterministic_Recon I1')
        plt.title('I1 comparison')
        pdf.savefig()
        plt.close()

        dfrecon_I2.index = pd.to_datetime(dfrecon_S.index)
        ax = dfrecon_I2.plot(color='silver', legend=False)
        dfsim_I2[dfsim_I2.index.isin(dfrecon_I2.index)].plot(ax=ax, color='black', legend='Simulated I2')
        dfdetrecon_I2[dfdetrecon_I2.index.isin(dfrecon_I2.index)].plot(ax=ax, color='red', legend='Deterministic_Recon I2')
        plt.title('I2 comparison')
        pdf.savefig()
        plt.close()

        dfrecon_I3.index = pd.to_datetime(dfrecon_S.index)
        ax = dfrecon_I3.plot(color='silver', legend=False)
        dfsim_I3[dfsim_I3.index.isin(dfrecon_I3.index)].plot(ax=ax, color='black', legend='Simulated I3')
        dfdetrecon_I3[dfdetrecon_I3.index.isin(dfrecon_I3.index)].plot(ax=ax, color='red', legend='Deterministic_Recon I3')
        plt.title('I3 comparison')
        pdf.savefig()
        plt.close()

        dfrecon_R.index = pd.to_datetime(dfrecon_R.index)
        ax = dfrecon_R.plot(color='silver', legend=False)
        dfsim_R[dfsim_R.index.isin(dfrecon_R.index)].plot(ax=ax, color='black', legend='Simulated R')
        dfdetrecon_R[dfdetrecon_R.index.isin(dfrecon_R.index)].plot(ax=ax, color='red', legend='Deterministic_Recon R')
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
    seiiir_fname = './input_data/FL_SEIIIR_R0_2.25_short_realization_4.csv'
    det_recon_fname = './data/FL/recon-det/recon_deterministic_delay_10.csv'
    recon_folder_name = './data/FL/recon-stoch'
    geodata_fname = './input_data/geodata.csv'
    output_path = './figures/'
    
    compare_florida_reconstruction(seiiir_fname, det_recon_fname, recon_folder_name, geodata_fname, output_path)
    
