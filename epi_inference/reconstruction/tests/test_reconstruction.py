import pytest
import pandas as pd
pd.set_option("display.max_rows", None)
import numpy as np
from datetime import datetime, timedelta
import epi_inference.reconstruction.common as rcommon
import epi_inference.reconstruction.deterministic as recond
import epi_inference.reconstruction.stochastic as recons
from epi_inference.simulation.simulation import simulate_discrete_seiiir_deterministic
from epi_inference.util import roundall
import matplotlib.pyplot as plt

"""
def test_force_keyword_args():
    @force_keyword_args
    def func(a,b):
        return '{}{}'.format(a,b)

    ans = func(a='A', b='B')
    assert ans == 'AB'

    with pytest.Raises(SyntaxError):
        func('A','B')
"""

#ToDo: switch the np stochastic reconstruction test to use simulated data from pipeline with our reporting model

def test_reported_cases_from_cumulative():
    # test handling of all zeros and that dates are correct
    dates = pd.date_range(end=datetime(year=2020, month=4, day=12), periods=4).to_pydatetime().tolist()
    assert dates[0] == datetime(year=2020, month=4, day=9)
    
    cc = [0, 0, 0, 0]    
    reported_cases_per_day \
        = rcommon.reported_cases_from_cumulative(dates=dates,
                                                 cumulative_reported_cases=cc)
    
    assert len(reported_cases_per_day.values) == len(cc)-1
    assert len(reported_cases_per_day.values) == 3
    assert len(reported_cases_per_day.dates) == len(dates)-1
    assert len(reported_cases_per_day.dates) == 3
    assert all(r == 0 for r in reported_cases_per_day.values)
    assert reported_cases_per_day.dates[0] == datetime(year=2020, month=4, day=10)
    expected_rc = [0, 0, 0]
    assert all([r == er for r, er in zip(reported_cases_per_day.values, expected_rc)])

    # test that we obtain expected numbers
    dates = pd.date_range(end=datetime(year=2020, month=4, day=3), periods=7).to_pydatetime().tolist()
    assert dates[0] == datetime(year=2020, month=3, day=28)
    cc = [0, 1, 2, 3, 3, 3, 5]
    reported_cases_per_day = \
        rcommon.reported_cases_from_cumulative(dates=dates,
                                               cumulative_reported_cases=cc)
    assert len(reported_cases_per_day.dates) == 6
    assert len(reported_cases_per_day.values) == 6
    expected_rc = [1, 1, 1, 0, 0, 2]
    assert all([r == er for r,er in zip(reported_cases_per_day.values, expected_rc)])
    assert reported_cases_per_day.dates[0] == datetime(year=2020, month=3, day=29)
    assert reported_cases_per_day.dates[-1] == datetime(year=2020, month=4, day=3)

    # test that we throw an error if dates are not the same length as cumulative reported cases
    dates = pd.date_range(end=datetime(year=2020, month=4, day=3), periods=5).to_pydatetime().tolist()
    cc = [0, 1, 2, 3, 3, 3, 5]
    with pytest.raises(ValueError):
        reported_cases_per_day = \
            rcommon.reported_cases_from_cumulative(dates=dates,
                                                   cumulative_reported_cases=cc)

    # test that we throw an error if cumulative reported cases does
    # not start at zero
    dates = pd.date_range(end=datetime(year=2020, month=4, day=3), periods=4).to_pydatetime().tolist()
    cc = [4, 5, 6, 7]
    with pytest.raises(ValueError):
        reported_cases_per_day = \
            rcommon.reported_cases_from_cumulative(dates=dates,
                                                   cumulative_reported_cases=cc)

    # that that we throw an error if cumulative cases are non-increasing
    dates = pd.date_range(end=datetime(year=2020, month=4, day=3), periods=5).to_pydatetime().tolist()
    cc = [0, 1, 2, 1, 5]
    with pytest.raises(ValueError):
        reported_cases_per_day = \
            rcommon.reported_cases_from_cumulative(dates=dates,
                                                   cumulative_reported_cases=cc)

def test_df_reported_cases_from_cumulative():
    # test handling of all zeros and that dates are correct
    dates = pd.date_range(end=datetime(year=2020, month=4, day=12), periods=4).to_pydatetime().tolist()
    assert dates[0] == datetime(year=2020, month=4, day=9)

    counties = ['12001', '12002', '12003']

    cumulative_reported_cases = dict()
    for i,c in enumerate(counties):
        cumulative_reported_cases[c] = [i*j for j in range(4)]

    cumulative_reported_cases['date'] = dates
    df_cumulative_reported_cases = pd.DataFrame(cumulative_reported_cases)

    df_reported_cases = rcommon.df_reported_cases_from_cumulative(df_cumulative_reported_cases)

    # test that we obtain the expected numbers
    reported_cases = dict()
    reported_cases['date'] = dates[1:]
    for i,c in enumerate(counties):
        reported_cases[c] = [i for j in range(3)]
    df_expected_reported_cases = pd.DataFrame(reported_cases).set_index('date')
    pd.testing.assert_frame_equal(df_reported_cases, df_expected_reported_cases.astype(float))

    # test that we throw an error if cumulative reported cases does
    # not start at zero
    cumulative_reported_cases = dict()
    for i,c in enumerate(['12001', '12002', '12003']):
        cumulative_reported_cases[c] = [i*j+1 for j in range(4)]

    cumulative_reported_cases['date'] = dates
    df_cumulative_reported_cases = pd.DataFrame(cumulative_reported_cases)

    with pytest.raises(ValueError):
        df_reported_cases = rcommon.df_reported_cases_from_cumulative(df_cumulative_reported_cases)

    # test that we throw an error if cumulative reported have decreasing
    # numbers
    cumulative_reported_cases = dict()
    for i,c in enumerate(['12001', '12002', '12003']):
        cumulative_reported_cases[c] = [i*j for j in range(4)]
    cumulative_reported_cases['12003'][2]=1

    cumulative_reported_cases['date'] = dates
    df_cumulative_reported_cases = pd.DataFrame(cumulative_reported_cases)

    with pytest.raises(ValueError):
        df_reported_cases = rcommon.df_reported_cases_from_cumulative(df_cumulative_reported_cases)

    """
    # test that we throw an error if dates are not the same length as cumulative reported cases
    dates = pd.date_range(end=datetime(year=2020, month=4, day=3), periods=5).to_pydatetime().tolist()
    cc = [0, 1, 2, 3, 3, 3, 5]
    with pytest.raises(ValueError):
        rdates, rc = rcommon.reported_cases_from_cumulative(dates, cc)
    """

def test_np_reported_cases_from_cumulative():
    # test handling of all zeros and that dates are correct
    dates = np.arange(datetime(2020,1,1), datetime(2020,1,5), timedelta(days=1)).astype(datetime)
    counties = np.asarray(['12001', '12002', '12003']).astype(object)

    cumulative_reported_cases = np.zeros((len(dates), len(counties)))
    for i,c in enumerate(counties):
        cumulative_reported_cases[:,i] = np.asarray([i*j for j in range(4)])

    reported_cases_per_day = \
        rcommon.np_reported_cases_from_cumulative(dates=dates,
                                                  counties=counties,
                                                  cumulative_reported_cases=cumulative_reported_cases)
    
    # test that we obtain the expected numbers
    np.testing.assert_array_equal(reported_cases_per_day.dates,dates[1:])
    np.testing.assert_array_equal(reported_cases_per_day.counties,counties)
    np_expected_reported_cases = np.zeros((len(reported_cases_per_day.dates), len(counties)))
    for i,c in enumerate(counties):
        np_expected_reported_cases[:,i] = np.asarray([i for j in range(3)])
    np.testing.assert_array_equal(reported_cases_per_day.values, np_expected_reported_cases)

    # test that we throw an error if dates are not the same length as cumulative reported cases
    with pytest.raises(ValueError):
        reported_cases_per_day = \
            rcommon.np_reported_cases_from_cumulative(dates=dates[2:],
                                                      counties=counties,
                                                      cumulative_reported_cases=cumulative_reported_cases)

    # test that we throw an error if cumulative reported cases do
    # not start at zero
    cumulative_reported_cases = np.zeros((len(dates), len(counties)))
    for i,c in enumerate(counties):
        cumulative_reported_cases[:,i] = np.asarray([i*j+1 for j in range(4)])

    with pytest.raises(ValueError):
        reported_cases_per_day = \
            rcommon.np_reported_cases_from_cumulative(dates=dates[2:],
                                                      counties=counties,
                                                      cumulative_reported_cases=cumulative_reported_cases)

    # test that we throw an error if cumulative reported have decreasing
    # numbers
    cumulative_reported_cases = np.zeros((len(dates), len(counties)))
    for i,c in enumerate(counties):
        cumulative_reported_cases[:,i] = np.asarray([i*j for j in range(4)])
    cumulative_reported_cases[2,2]=1

    with pytest.raises(ValueError):
        reported_cases_per_day = \
            rcommon.np_reported_cases_from_cumulative(dates=dates,
                                                      counties=counties,
                                                      cumulative_reported_cases=cumulative_reported_cases)

def test_transmissions_from_reported_cases():
    # run a simulation and test that we recover the transmission
    # vector from the cumulative reported cases
    N = 1000
    y0={'S': N, 'E': 5, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5.2
    gamma = 1/4
    beta = 2.2*gamma
    reporting_factor = 8
    report_delay = 7
    tf=100
    
    sim = simulate_discrete_seiiir_deterministic(y0, tf, beta, sigma, gamma, reporting_factor, N, report_delay, tx=None)
    SEIIIR = sim.SEIIIR

    reported_cases_per_day = rcommon.reported_cases_from_cumulative(dates=sim.cumulative_reported_cases.dates,
                                                                    cumulative_reported_cases=sim.cumulative_reported_cases.values
    )

    assert len(reported_cases_per_day.dates) == len(sim.cumulative_reported_cases.dates)-1
    assert len(reported_cases_per_day.dates) == len(sim.SEIIIR.dates)
    assert len(reported_cases_per_day.values) == len(sim.cumulative_reported_cases.values)-1
    
    transmissions = recond._transmissions_from_reported_cases(dates=reported_cases_per_day.dates,
                                                              reported_cases_per_day=reported_cases_per_day.values,
                                                              reporting_factor=reporting_factor,
                                                              report_delay=report_delay
    )
    assert len(transmissions.dates) == len(reported_cases_per_day.dates)
    assert len(transmissions.values) == len(reported_cases_per_day.values)
    
    np.testing.assert_allclose(np.asarray(SEIIIR.transmissions), np.asarray(transmissions.values))
    for i in range(len(transmissions.dates)):
        assert transmissions.dates[i] == SEIIIR.dates[i]

def test_np_deterministic():
    # run some simulations and test that we recover the transmissions
    # from the cumulative reported cases
    N = 1000
    sigma = 1/5.2
    gamma = 1/4
    beta = 2.2*gamma
    reporting_factor = 8
    report_delay = 7
    tf=100
    tx=[0]*tf
    tx[30]=1

    counties = np.asarray(['12001', '12002', '12003'], dtype=np.object)
    populations = np.zeros(len(counties))
    cumulative_reported_cases = np.zeros((tf+1,len(counties)))
    simT = np.zeros((tf,len(counties)))
    simS = np.zeros((tf,len(counties)))
    simE = np.zeros((tf,len(counties)))
    simI1 = np.zeros((tf,len(counties)))
    simI2 = np.zeros((tf,len(counties)))
    simI3 = np.zeros((tf,len(counties)))
    simR = np.zeros((tf,len(counties)))
    for i,c in enumerate(counties):
        y0={'S': N*(i+1), 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
        populations[i] = N*(i+1)
        sim = simulate_discrete_seiiir_deterministic(y0, tf, beta, sigma, gamma, reporting_factor, populations[i], report_delay, tx=tx)
        cumulative_reported_cases[:,i] = sim.cumulative_reported_cases.values
        simT[:,i] = sim.SEIIIR.transmissions
        simS[:,i] = sim.SEIIIR.S
        simE[:,i] = sim.SEIIIR.E
        simI1[:,i] = sim.SEIIIR.I1
        simI2[:,i] = sim.SEIIIR.I2
        simI3[:,i] = sim.SEIIIR.I3
        simR[:,i] = sim.SEIIIR.R
        cumulative_reported_cases_dates = np.asarray(sim.cumulative_reported_cases.dates)
        dates = np.asarray(sim.SEIIIR.dates)

    reported_cases_per_day = \
        rcommon.np_reported_cases_from_cumulative(dates=cumulative_reported_cases_dates,
                                                  counties=counties,
                                                  cumulative_reported_cases=cumulative_reported_cases)

    assert len(reported_cases_per_day.dates) == len(cumulative_reported_cases_dates)-1
    assert len(reported_cases_per_day.dates) == len(dates)
    assert len(reported_cases_per_day.values) == len(dates)
    np.testing.assert_array_equal(reported_cases_per_day.dates, cumulative_reported_cases_dates[1:])
    np.testing.assert_array_equal(reported_cases_per_day.counties, counties)

    ### test the intermediate method
    transmissions = recond._np_transmissions_from_reported_cases(dates=reported_cases_per_day.dates,
                                                                 counties=reported_cases_per_day.counties,
                                                                 reported_cases_per_day=reported_cases_per_day.values,
                                                                 reporting_factor=reporting_factor,
                                                                 report_delay=report_delay
    )
    assert len(transmissions.dates) == len(reported_cases_per_day.dates)
    assert reported_cases_per_day.values.shape == transmissions.values.shape
    
    np.testing.assert_allclose(simT, transmissions.values)
    np.testing.assert_array_equal(transmissions.dates, dates)

    ### test the reconstruction
    recon = recond.np_reconstruct_states_deterministic_decay(dates=reported_cases_per_day.dates,
                                                             counties=reported_cases_per_day.counties,
                                                             reported_cases_per_day=reported_cases_per_day.values,
                                                             populations=populations,
                                                             sigma=sigma,
                                                             gamma=gamma,
                                                             reporting_factor=reporting_factor,
                                                             report_delay=report_delay
                                                             )

    np.testing.assert_array_equal(dates, recon.dates)
    np.testing.assert_array_almost_equal(simS,recon.S)
    np.testing.assert_array_almost_equal(simE,recon.E)
    np.testing.assert_array_almost_equal(simI1,recon.I1)
    np.testing.assert_array_almost_equal(simI2,recon.I2)
    np.testing.assert_array_almost_equal(simI3,recon.I3)
    np.testing.assert_array_almost_equal(simR,recon.R)

def test_reconstruct_states_deterministic_decay():
    # run a simulation and test that we reconstruct the states
    # from the cumulative reported cases
    N = 1000
    y0={'S': N, 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5.2
    gamma = 1/4
    beta = 2.2*gamma
    reporting_factor = 8
    report_delay = 7
    tf=100
    tx=[0]*tf
    tx[10] = 1
    
    # run the simulation
    sim = simulate_discrete_seiiir_deterministic(y0, tf, beta, sigma, gamma, reporting_factor, N, report_delay, tx=tx)
    SEIIIR = sim.SEIIIR

    # reconstruct the states
    reported_cases_per_day = rcommon.reported_cases_from_cumulative(dates=sim.cumulative_reported_cases.dates,
                                                                    cumulative_reported_cases=sim.cumulative_reported_cases.values)
    
    recon = recond.reconstruct_states_deterministic_decay(dates=reported_cases_per_day.dates,
                                                          reported_cases_per_day=reported_cases_per_day.values,
                                                          population=N, sigma=sigma,
                                                          gamma=gamma, reporting_factor=reporting_factor,
                                                          report_delay=report_delay)

    assert len(reported_cases_per_day.dates) == len(sim.SEIIIR.dates)
    assert len(recon.dates) == len(sim.SEIIIR.dates)
    assert len(reported_cases_per_day.values) == len(sim.SEIIIR.dates)
    assert len(recon.transmissions) == len(sim.SEIIIR.transmissions)
    assert len(recon.S) == len(SEIIIR.S)
    assert len(recon.E) == len(SEIIIR.E)
    assert len(recon.I1) == len(SEIIIR.I1)
    assert len(recon.I2) == len(SEIIIR.I2)
    assert len(recon.I3) == len(SEIIIR.I3)
    assert len(recon.R) == len(SEIIIR.R)

    np.testing.assert_allclose(np.asarray(SEIIIR.transmissions), np.asarray(recon.transmissions))
    np.testing.assert_allclose(np.asarray(SEIIIR.S), np.asarray(recon.S))
    np.testing.assert_allclose(np.asarray(SEIIIR.E), np.asarray(recon.E))
    np.testing.assert_allclose(np.asarray(SEIIIR.I1), np.asarray(recon.I1))
    np.testing.assert_allclose(np.asarray(SEIIIR.I2), np.asarray(recon.I2))
    np.testing.assert_allclose(np.asarray(SEIIIR.I3), np.asarray(recon.I3))
    np.testing.assert_allclose(np.asarray(SEIIIR.R), np.asarray(recon.R))

    for i in range(len(recon.dates)):
        assert recon.dates[i] == SEIIIR.dates[i]

    for i in range(len(reported_cases_per_day.dates)):
        assert reported_cases_per_day.dates[i] == sim.cumulative_reported_cases.dates[i+1]

def test_stochastic_reconstruction():
    # run a simulation and test that we reconstruct the states
    # from the cumulative reported cases
    N = 1000
    y0={'S': N, 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5.2
    gamma = 1/4.3
    beta = 2.2*gamma
    reporting_factor = 10
    report_delay = 8
    tf=100
    tx=[0]*tf
    tx[10] = 1
    
    # run the simulation
    sim = simulate_discrete_seiiir_deterministic(y0, tf, beta, sigma, gamma, reporting_factor, N, report_delay, tx=tx)

    # reconstruct the states
    dfT = None
    dfS = None
    dfE = None
    dfI1 = None
    dfI2 = None
    dfI3 = None
    dfR = None
    np.random.seed(1975)
    for real in range(50):
        reported_cases_per_day = rcommon.reported_cases_from_cumulative(dates=sim.cumulative_reported_cases.dates,
                                                                        cumulative_reported_cases=sim.cumulative_reported_cases.values)

        recon = recons.stochastic_reconstruction(dates=reported_cases_per_day.dates,
                                                 reported_cases_per_day=reported_cases_per_day.values,                                       
                                                 population=N,
                                                 n_steps_per_day=4,
                                                 reporting_delay_mean=8,
                                                 reporting_delay_dev=1.35,
                                                 reporting_multiplier=reporting_factor,
                                                 fixed_incubation=5.2,
                                                 infectious_lower=2.6,
                                                 infectious_upper=6.0)

        if dfT is None:
            dfT = pd.DataFrame({'dates': recon.dates}).set_index('dates')
            dfS = pd.DataFrame({'dates': recon.dates}).set_index('dates')
            dfE = pd.DataFrame({'dates': recon.dates}).set_index('dates')
            dfI1 = pd.DataFrame({'dates': recon.dates}).set_index('dates')
            dfI2 = pd.DataFrame({'dates': recon.dates}).set_index('dates')
            dfI3 = pd.DataFrame({'dates': recon.dates}).set_index('dates')
            dfR = pd.DataFrame({'dates': recon.dates}).set_index('dates')

        dfT['{}'.format(real)] = recon.transmissions
        dfS['{}'.format(real)] = recon.S
        dfE['{}'.format(real)] = recon.E
        dfI1['{}'.format(real)] = recon.I1
        dfI2['{}'.format(real)] = recon.I2
        dfI3['{}'.format(real)] = recon.I3
        dfR['{}'.format(real)] = recon.R

    dfsimT = pd.DataFrame({'dates':sim.SEIIIR.dates, 'sim':sim.SEIIIR.transmissions}).set_index('dates')
    assert_if_mean_significantly_different(dfT, dfsimT, 4, 0.2)
    dfsimS = pd.DataFrame({'dates':sim.SEIIIR.dates, 'sim':sim.SEIIIR.S}).set_index('dates')
    assert_if_mean_significantly_different(dfS, dfsimS, 4, 0.25)
    dfsimE = pd.DataFrame({'dates':sim.SEIIIR.dates, 'sim':sim.SEIIIR.E}).set_index('dates')
    assert_if_mean_significantly_different(dfE, dfsimE, 4, 0.2)
    dfsimI1 = pd.DataFrame({'dates':sim.SEIIIR.dates, 'sim':sim.SEIIIR.I1}).set_index('dates')
    assert_if_mean_significantly_different(dfI1, dfsimI1, 4, 0.2)
    dfsimI2 = pd.DataFrame({'dates':sim.SEIIIR.dates, 'sim':sim.SEIIIR.I2}).set_index('dates')
    assert_if_mean_significantly_different(dfI2, dfsimI2, 4, 0.2)
    dfsimI3 = pd.DataFrame({'dates':sim.SEIIIR.dates, 'sim':sim.SEIIIR.I3}).set_index('dates')
    assert_if_mean_significantly_different(dfI3, dfsimI3, 4, 0.2)
    dfsimR = pd.DataFrame({'dates':sim.SEIIIR.dates, 'sim':sim.SEIIIR.R}).set_index('dates')
    assert_if_mean_significantly_different(dfR, dfsimR, 4, 0.2)

@pytest.mark.skip('Not updating numpy version for now')
def test_np_stochastic_reconstruction():
    # run a simulation and test that we reconstruct the states
    # from the cumulative reported cases
    N = 10000
    sigma = 1/5.2
    gamma = 1/4.3
    beta = 2.2*gamma
    reporting_factor = 10
    report_delay = 8
    tf=100
    tx=[0]*tf
    tx[10] = 1

    counties = np.asarray(['12001', '12002', '12003'], dtype=np.object)
    #counties = np.asarray(['12001'], dtype=np.object)
    populations = np.zeros(len(counties))
    cumulative_reported_cases = np.zeros((tf+1,len(counties)))
    simT = np.zeros((tf,len(counties)))
    simS = np.zeros((tf,len(counties)))
    simE = np.zeros((tf,len(counties)))
    simI1 = np.zeros((tf,len(counties)))
    simI2 = np.zeros((tf,len(counties)))
    simI3 = np.zeros((tf,len(counties)))
    simR = np.zeros((tf,len(counties)))
    for i,c in enumerate(counties):
        y0={'S': N*(i+1), 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
        populations[i] = N*(i+1)
        sim = simulate_discrete_seiiir_deterministic(y0, tf, beta, sigma, gamma, reporting_factor, populations[i], report_delay, tx=tx)
        cumulative_reported_cases[:,i] = sim.cumulative_reported_cases.values
        simT[:,i] = sim.SEIIIR.transmissions
        simS[:,i] = sim.SEIIIR.S
        simE[:,i] = sim.SEIIIR.E
        simI1[:,i] = sim.SEIIIR.I1
        simI2[:,i] = sim.SEIIIR.I2
        simI3[:,i] = sim.SEIIIR.I3
        simR[:,i] = sim.SEIIIR.R
        cumulative_reported_cases_dates = np.asarray(sim.cumulative_reported_cases.dates)
        dates = np.asarray(sim.SEIIIR.dates)

    # reconstruct the states
    np.random.seed(42)
    for real in range(100):
        reported_cases_per_day = rcommon.np_reported_cases_from_cumulative(dates=cumulative_reported_cases_dates,
                                                                           counties=counties,
                                                                           cumulative_reported_cases=cumulative_reported_cases)

        recon = recons.np_stochastic_reconstruction(dates=reported_cases_per_day.dates,
                                                    counties=reported_cases_per_day.counties,
                                                    reported_cases_per_day=reported_cases_per_day.values,
                                                    populations=populations,
                                                    n_steps_per_day=4)

        realization_df = _long_dataframe_from_recon(recon)
        realization_df['realization'] = real
        if real == 0:
            recon_df = realization_df
        else:
            recon_df = pd.concat([recon_df, realization_df])

    abstol=4
    dfsimS = pd.DataFrame(index=dates, columns=counties, data=simS)
    _special_tolerance_check(dfsimS, recon_df, 'S', abstol, 0.05)
    dfsimE = pd.DataFrame(index=dates, columns=counties, data=simE)
    _special_tolerance_check(dfsimE, recon_df, 'E', abstol, 0.15)
    dfsimI1 = pd.DataFrame(index=dates, columns=counties, data=simI1)
    _special_tolerance_check(dfsimI1, recon_df, 'I1', abstol, 0.2)
    dfsimI2 = pd.DataFrame(index=dates, columns=counties, data=simI2)
    _special_tolerance_check(dfsimI2, recon_df, 'I2', abstol, 0.2)
    dfsimI3 = pd.DataFrame(index=dates, columns=counties, data=simI3)
    _special_tolerance_check(dfsimI3, recon_df, 'I3', abstol, 0.2)
    dfsimR = pd.DataFrame(index=dates, columns=counties, data=simR)
    _special_tolerance_check(dfsimR, recon_df, 'R', abstol, 0.2)
    
def assert_if_mean_significantly_different(df1, df2, abstol, reltol):
    mean = df1.mean(axis=1)
    mean.name = 'mean'
    dferrors = df2.join(mean, how='inner')
    dferrors['errors'] = (dferrors['mean']-dferrors['sim']).abs()
    dferrors['relerrors'] = dferrors['errors'].divide((dferrors['sim']+1))
    # print(dferrors['relerrors'])
    dferrors['flag'] = dferrors['errors'].gt(abstol) & dferrors['relerrors'].gt(reltol)
    nerrs = dferrors['flag'].sum()
    if nerrs > 0:
        print(df1)
        print(df2)
        print(dferrors)
    assert nerrs == 0

def _special_tolerance_check(dfsim, recon_df, comp, abstol, reltol):
    dfsim.index.name = 'dates'
    mean = recon_df[recon_df['comp']==comp].groupby(['dates','county']).mean().reset_index().pivot(index='dates', columns='county', values='count')
    diff = (mean-dfsim).dropna(how='all')
    reldiff = diff.divide(dfsim[dfsim.index.isin(diff.index)])
    # print('...', comp)
    # print(diff)
    # print(reldiff)
    nerrs = (diff.gt(abstol) & reldiff.gt(reltol)).sum().sum()
    assert nerrs == 0

def _quantile_check(dfsim, recon_df, comp, counties, lowerq, upperq):
    dfsim.index.name = 'dates'
    # loop over all the counties
    for i,c in enumerate(counties):
        df = recon_df[recon_df['comp']==comp]
        df = df[df['county']==c]
        df = df.pivot(index='dates', columns='realization', values='count')
        lower = df.quantile(lowerq, axis=1)[df.index.isin(dfsim.index)]
        upper = df.quantile(upperq, axis=1)[df.index.isin(dfsim.index)]
        sim = dfsim[c]
        sim = sim[sim.index.isin(lower.index)]
        # print(lower)
        # print(sim)
        # print(upper)
        # print((lower.gt(sim) | sim.gt(upper)))
        nerrs = (lower.gt(sim) | sim.gt(upper)).sum()
        print('Checking compartment:', comp, ' for county: ', c)
        assert nerrs == 0

def _wide_dataframe_from_recon(recon):
    S = pd.DataFrame(index=recon.dates, columns=recon.counties, data=recon.S)
    S['comp'] = 'S'
    E = pd.DataFrame(index=recon.dates, columns=recon.counties, data=recon.E)
    E['comp'] = 'E'
    I1 = pd.DataFrame(index=recon.dates, columns=recon.counties, data=recon.I1)
    I1['comp'] = 'I1'
    I2 = pd.DataFrame(index=recon.dates, columns=recon.counties, data=recon.I2)
    I2['comp'] = 'I2'
    I3 = pd.DataFrame(index=recon.dates, columns=recon.counties, data=recon.I3)
    I3['comp'] = 'I3'
    R = pd.DataFrame(index=recon.dates, columns=recon.counties, data=recon.R)
    R['comp'] = 'R'
    return pd.concat([S, E, I1, I2, I3, R])

def _long_dataframe_from_recon(recon):
    df = None
    for i,c in enumerate(recon.counties):
        S = pd.DataFrame(index=recon.dates, columns=['count'], data=recon.S[:,i])
        S.index.name = 'dates'
        S.reset_index(inplace=True)
        S['comp'] = 'S'
        S['county'] = c
        E = pd.DataFrame(index=recon.dates, columns=['count'], data=recon.E[:,i])
        E.index.name = 'dates'
        E.reset_index(inplace=True)
        E['comp'] = 'E'
        E['county'] = c
        I1 = pd.DataFrame(index=recon.dates, columns=['count'], data=recon.I1[:,i])
        I1.index.name = 'dates'
        I1.reset_index(inplace=True)
        I1['comp'] = 'I1'
        I1['county'] = c
        I2 = pd.DataFrame(index=recon.dates, columns=['count'], data=recon.I2[:,i])
        I2.index.name = 'dates'
        I2.reset_index(inplace=True)
        I2['comp'] = 'I2'
        I2['county'] = c
        I3 = pd.DataFrame(index=recon.dates, columns=['count'], data=recon.I3[:,i])
        I3.index.name = 'dates'
        I3.reset_index(inplace=True)
        I3['comp'] = 'I3'
        I3['county'] = c
        R = pd.DataFrame(index=recon.dates, columns=['count'], data=recon.R[:,i])
        R.index.name = 'dates'
        R.reset_index(inplace=True)
        R['comp'] = 'R'
        R['county'] = c

        if df is None:
            df = pd.concat([S, E, I1, I2, I3, R])
        else:
            df = pd.concat([df, S, E, I1, I2, I3, R])
    return df
