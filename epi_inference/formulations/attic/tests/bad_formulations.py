import pytest
import numpy as np
import pandas as pd
from datetime import datetime

from epi_inference.simulation.simulation import simulate_discrete_seiiir_deterministic
from epi_inference.util import roundall
from epi_inference.formulations.decay_lsq import run_decay_lsq
from epi_inference.formulations.decay_blike import run_decay_blike
#from epi_inference.formulations.multinode_decay_lsq import run_multinode_decay_lsq
#import matplotlib.pyplot as plt

# Todo: Add tests for other formulations

def test_decay_lsq_standard():
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

    for beta in [0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5]:
        Cdates,C,dates,T,S,E,I1,I2,I3,R = simulate_discrete_seiiir_deterministic(y0, tf, beta=beta,
                                                                        sigma=sigma, gamma=gamma,
                                                                        rho=rho, N=N,
                                                                        report_delay=report_delay,
                                                                        tx=tx)

        results = run_decay_lsq(cm_rep_cases=C,
                                population=N,
                                sigma=sigma,
                                gamma=gamma,
                                report_delay=report_delay,
                                analysis_window=dict(),
                                reporting_factor=rho,
                                Cdates=Cdates
        )

        assert beta == pytest.approx(results['est_beta'])
        assert results['status'] == 'ok'
        assert results['msg'] == 'Optimal solution found'

def test_decay_lsq_analysis_window():
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

    for beta in [0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5]:
        Cdates,C,dates,T,S,E,I1,I2,I3,R = simulate_discrete_seiiir_deterministic(y0, tf, beta=beta,
                                                                        sigma=sigma, gamma=gamma,
                                                                        rho=rho, N=N,
                                                                        report_delay=report_delay,
                                                                        tx=tx)

        results = run_decay_lsq(cm_rep_cases=C,
                                population=N,
                                sigma=sigma,
                                gamma=gamma,
                                report_delay=report_delay,
                                analysis_window={'days':10},
                                reporting_factor=rho,
                                Cdates=Cdates
        )

        assert beta == pytest.approx(results['est_beta'])
        assert results['status'] == 'ok'
        assert results['msg'] == 'Optimal solution found'

def test_decay_lsq_nonunity_reporting_factor():
    # test with non-unity reporting factor
    N = 1000000
    y0={'S': N, 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5
    gamma = 1/4
    rho = 8
    report_delay = 7
    tf = 25
    tx = [0]*tf
    tx[10] = 1
    for beta in [0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5]:
        Cdates,C,dates,T,S,E,I1,I2,I3,R = simulate_discrete_seiiir_deterministic(y0, tf, beta=beta,
                                                                        sigma=sigma, gamma=gamma,
                                                                        rho=rho, N=N,
                                                                        report_delay=report_delay,
                                                                        tx=tx)

        results = run_decay_lsq(cm_rep_cases=C,
                                population=N,
                                sigma=sigma,
                                gamma=gamma,
                                analysis_window=dict(),
                                report_delay=report_delay,
                                reporting_factor=rho,
                                Cdates=Cdates
        )

        assert beta == pytest.approx(results['est_beta'])
        assert results['status'] == 'ok'
        assert results['msg'] == 'Optimal solution found'

def test_decay_lsq_incorrect_reporting_factor():
    # test with incorrect reporting factor
    # here, we expect some error, especially with larger beta values (since the S/N
    # term is inaccurate)
    N = 1000000
    y0={'S': N, 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5
    gamma = 1/4
    rho_assumed = 2
    rho_actual = 8
    report_delay = 7
    tf = 25
    tx = [0]*tf
    tx[10] = 1

    for beta in [0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5]:
        Cdates,C,dates,T,S,E,I1,I2,I3,R = simulate_discrete_seiiir_deterministic(y0, tf, beta=beta,
                                                                        sigma=sigma, gamma=gamma,
                                                                        rho=rho_actual, N=N,
                                                                        report_delay=report_delay,
                                                                        tx=tx)
    
        results = run_decay_lsq(cm_rep_cases=C,
                                population=N,
                                sigma=sigma,
                                gamma=gamma,
                                analysis_window=dict(),
                                report_delay=report_delay,
                                reporting_factor=rho_assumed,
                                Cdates=Cdates
        )

        assert beta == pytest.approx(results['est_beta'],1e-2)
        assert results['status'] == 'ok'
        assert results['msg'] == 'Optimal solution found'

def test_decay_lsq_zero_data():
    # test that the status flags are correct when the solve does not set a value for beta
    N = 1000000
    y0={'S': N, 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5
    gamma = 1/4
    rho = 1
    report_delay = 7
    tf = 25
    tx = [0]*tf
    tx[10] = 1

    # test the inference with no data
    zeroC = [0]*(tf+1)
    Cdates = pd.date_range(end=datetime(year=2020, month=4, day=12), periods=tf+1).to_pydatetime().tolist()
    results = run_decay_lsq(cm_rep_cases=zeroC,
                            population=N,
                            sigma=sigma,
                            gamma=gamma,
                            report_delay=report_delay,
                            analysis_window=dict(),
                            reporting_factor=rho,
                            Cdates=Cdates
        )

    assert results['est_beta'] is None
    assert results['status'] == 'stale'
    assert results['msg'] == 'Transmission parameter beta not solved (stale).'

def test_decay_blike_standard():
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

    for beta in [0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5]:
        Cdates,C,dates,T,S,E,I1,I2,I3,R = simulate_discrete_seiiir_deterministic(y0, tf, beta=beta,
                                                                        sigma=sigma, gamma=gamma,
                                                                        rho=rho, N=N,
                                                                        report_delay=report_delay,
                                                                        tx=tx)

        results = run_decay_blike(cm_rep_cases=C,
                                  population=N,
                                  sigma=sigma,
                                  gamma=gamma,
                                  report_delay=report_delay,
                                  analysis_window=dict(),
                                  reporting_factor=rho,
                                  Cdates=Cdates
        )

        # Note, we expect errors here since the models are not the same
        assert beta == pytest.approx(results['est_beta'], 1e-2)
        assert results['status'] == 'ok'
        assert results['msg'] == 'Optimal solution found'

def test_decay_blike_analysis_window():
    N = 1000000
    y0={'S': N, 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5
    gamma = 1/4
    rho = 1
    report_delay = 7
    tf = 25
    tx = [0]*tf
    tx[10] = 1

    for beta in [0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5]:
        Cdates,C,dates,T,S,E,I1,I2,I3,R = simulate_discrete_seiiir_deterministic(y0, tf, beta=beta,
                                                                        sigma=sigma, gamma=gamma,
                                                                        rho=rho, N=N,
                                                                        report_delay=report_delay,
                                                                        tx=tx)

        results = run_decay_blike(cm_rep_cases=C,
                                  population=N,
                                  sigma=sigma,
                                  gamma=gamma,
                                  report_delay=report_delay,
                                  analysis_window={'days':10},
                                  reporting_factor=rho,
                                  Cdates=Cdates
        )

        # Note, we expect errors here since the models are not the same
        assert beta == pytest.approx(results['est_beta'], 1e-2)
        assert results['status'] == 'ok'
        assert results['msg'] == 'Optimal solution found'

def test_decay_blike_nonunity_reporting_factor():
    # test with non-unity reporting factor
    N = 1000000
    y0={'S': N, 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5
    gamma = 1/4
    rho = 8
    report_delay = 7
    tf = 25
    tx = [0]*tf
    tx[10] = 1
    for beta in [0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5]:
        Cdates,C,dates,T,S,E,I1,I2,I3,R = simulate_discrete_seiiir_deterministic(y0, tf, beta=beta,
                                                                        sigma=sigma, gamma=gamma,
                                                                        rho=rho, N=N,
                                                                        report_delay=report_delay,
                                                                        tx=tx)
    
        results = run_decay_blike(cm_rep_cases=C,
                                  population=N,
                                  sigma=sigma,
                                  gamma=gamma,
                                  analysis_window=dict(),
                                  report_delay=report_delay,
                                  reporting_factor=rho,
                                  Cdates=Cdates
        )
        
        # Note, we expect errors here since the models are not the same
        assert beta == pytest.approx(results['est_beta'], 1e-2)
        assert results['status'] == 'ok'
        assert results['msg'] == 'Optimal solution found'

def test_decay_blike_incorrect_reporting_factor():
    # test with incorrect reporting factor
    # here, we expect some error, especially with larger beta values (since the S/N
    # term is inaccurate)
    N = 1000000
    y0={'S': N, 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5
    gamma = 1/4
    rho_assumed = 2
    rho_actual = 8
    report_delay = 7
    tf = 25
    tx = [0]*tf
    tx[10] = 1

    for beta in [0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5]:
        Cdates,C,dates,T,S,E,I1,I2,I3,R = simulate_discrete_seiiir_deterministic(y0, tf, beta=beta,
                                                                        sigma=sigma, gamma=gamma,
                                                                        rho=rho_actual, N=N,
                                                                        report_delay=report_delay,
                                                                        tx=tx)
    
        results = run_decay_lsq(cm_rep_cases=C,
                                population=N,
                                sigma=sigma,
                                gamma=gamma,
                                analysis_window=dict(),
                                report_delay=report_delay,
                                reporting_factor=rho_assumed,
                                Cdates=Cdates
        )

        # Note, we expect errors here since the models are not the same
        assert beta == pytest.approx(results['est_beta'], 1e-2)
        assert results['status'] == 'ok'
        assert results['msg'] == 'Optimal solution found'

def test_decay_blike_zero_data():
    # test that the status flags are correct when the solve does not set a value for beta
    N = 1000000
    y0={'S': N, 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5
    gamma = 1/4
    rho = 1
    report_delay = 7
    tf = 25
    tx = [0]*tf
    tx[10] = 1

    # test the inference with no data
    zeroC = [0]*(tf+1)
    Cdates = pd.date_range(end=datetime(year=2020, month=4, day=12), periods=tf+1).to_pydatetime().tolist()
    results = run_decay_blike(cm_rep_cases=zeroC,
                              population=N,
                              sigma=sigma,
                              gamma=gamma,
                              report_delay=report_delay,
                              analysis_window=dict(),
                              reporting_factor=rho,
                              Cdates=Cdates
    )

    assert results['est_beta'] is None
    assert results['status'] == 'failed'
    assert results['msg'] == 'Transmission parameter beta not solved (stale) in least-squares initialization.'

def test_multinode_decay_lsq_standard():
    """
    Test the multinode decay inference using data from a simulation with the seiiir deterministic model
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

    for beta in [0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5]:
        Cdates,C,dates,T,S,E,I1,I2,I3,R = simulate_discrete_seiiir_deterministic(y0, tf, beta=beta,
                                                                        sigma=sigma, gamma=gamma,
                                                                        rho=rho, N=N,
                                                                        report_delay=report_delay,
                                                                        tx=tx)

        # create a dataframe with the column of cases
        df_cm_rep_cases = pd.DataFrame({'0': list(C), '1': list(C), '3': list(C)})
        df_populations = pd.DataFrame({'county': ['0', '1', '2'], 'population':[N, N, N]})
        df_populations = df_populations.set_index('county')
        results = run_multinode_decay_lsq(
            cm_rep_cases=df_cm_rep_cases,
            populations=df_populations,
            sigma=sigma,
            gamma=gamma,
            report_delay=report_delay,
            analysis_window=dict(),
            reporting_factor=rho,
            Cdates=Cdates
        )

        # here, all three columns are the same, so we would expect beta to match
        assert beta == pytest.approx(results['est_beta'])
        assert results['status'] == 'ok'
        assert results['msg'] == 'Optimal solution found'

@pytest.mark.skip('analysis_window not yet implemented for multinode decay')
def test_multinode_decay_analysis_window():
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

    for beta in [0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5]:
        Cdates,C,dates,T,S,E,I1,I2,I3,R = simulate_discrete_seiiir_deterministic(y0, tf, beta=beta,
                                                                        sigma=sigma, gamma=gamma,
                                                                        rho=rho, N=N,
                                                                        report_delay=report_delay,
                                                                        tx=tx)

        # create a dataframe with the column of cases
        df_cm_rep_cases = pd.DataFrame({'0': list(C), '1': list(C), '3': list(C)})
        df_populations = pd.DataFrame({'county': ['0', '1', '2'], 'population':[N, N, N]})
        df_populations = df_populations.set_index('county')
        results = run_multinode_decay_lsq(
            cm_rep_cases=df_cm_rep_cases,
            populations=df_populations,
            sigma=sigma,
            gamma=gamma,
            report_delay=report_delay,
            analysis_window={'days':10},
            reporting_factor=rho,
            Cdates=Cdates
        )

        assert beta == pytest.approx(results['est_beta'])
        assert results['status'] == 'ok'
        assert results['msg'] == 'Optimal solution found'

def test_multinode_decay_lsq_nonunity_reporting_factor():
    # test with non-unity reporting factor
    N = 1000000
    y0={'S': N, 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5
    gamma = 1/4
    rho = 8
    report_delay = 7
    tf = 25
    tx = [0]*tf
    tx[10] = 1
    for beta in [0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5]:
        Cdates,C,dates,T,S,E,I1,I2,I3,R = simulate_discrete_seiiir_deterministic(y0, tf, beta=beta,
                                                                        sigma=sigma, gamma=gamma,
                                                                        rho=rho, N=N,
                                                                        report_delay=report_delay,
                                                                        tx=tx)

        # create a dataframe with the column of cases
        df_cm_rep_cases = pd.DataFrame({'0': list(C), '1': list(C), '3': list(C)})
        df_populations = pd.DataFrame({'county': ['0', '1', '2'], 'population':[N, N, N]})
        df_populations = df_populations.set_index('county')
        results = run_multinode_decay_lsq(
            cm_rep_cases=df_cm_rep_cases,
            populations=df_populations,
            sigma=sigma,
            gamma=gamma,
            report_delay=report_delay,
            analysis_window=dict(),
            reporting_factor=rho,
            Cdates=Cdates
        )

        assert beta == pytest.approx(results['est_beta'])
        assert results['status'] == 'ok'
        assert results['msg'] == 'Optimal solution found'

def test_multinode_decay_lsq_incorrect_reporting_factor():
    # test with incorrect reporting factor
    # here, we expect some error, especially with larger beta values (since the S/N
    # term is inaccurate)
    N = 1000000
    y0={'S': N, 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5
    gamma = 1/4
    rho_assumed = 2
    rho_actual = 8
    report_delay = 7
    tf = 25
    tx = [0]*tf
    tx[10] = 1

    for beta in [0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5]:
        Cdates,C,dates,T,S,E,I1,I2,I3,R = simulate_discrete_seiiir_deterministic(y0, tf, beta=beta,
                                                                        sigma=sigma, gamma=gamma,
                                                                        rho=rho_actual, N=N,
                                                                        report_delay=report_delay,
                                                                        tx=tx)
    
        # create a dataframe with the column of cases
        df_cm_rep_cases = pd.DataFrame({'0': list(C), '1': list(C), '3': list(C)})
        df_populations = pd.DataFrame({'county': ['0', '1', '2'], 'population':[N, N, N]})
        df_populations = df_populations.set_index('county')
        results = run_multinode_decay_lsq(
            cm_rep_cases=df_cm_rep_cases,
            populations=df_populations,
            sigma=sigma,
            gamma=gamma,
            report_delay=report_delay,
            analysis_window=dict(),
            reporting_factor=2,
            Cdates=Cdates
        )

        assert beta == pytest.approx(results['est_beta'],1e-2)
        assert results['status'] == 'ok'
        assert results['msg'] == 'Optimal solution found'

def test_multinode_decay_lsq_zero_data():
    # test that the status flags are correct when the solve does not set a value for beta
    N = 1000000
    y0={'S': N, 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5
    gamma = 1/4
    rho = 1
    report_delay = 7
    tf = 25
    tx = [0]*tf
    tx[10] = 1

    # test the inference with no data
    zeroC = [0]*(tf+1)
    Cdates = pd.date_range(end=datetime(year=2020, month=4, day=12), periods=tf+1).to_pydatetime().tolist()
    # create a dataframe with the column of cases
    df_cm_rep_cases = pd.DataFrame({'0': list(zeroC), '1': list(zeroC), '3': list(zeroC)})
    df_populations = pd.DataFrame({'county': ['0', '1', '2'], 'population':[N, N, N]})
    df_populations = df_populations.set_index('county')
    results = run_multinode_decay_lsq(
        cm_rep_cases=df_cm_rep_cases,
        populations=df_populations,
        sigma=sigma,
        gamma=gamma,
        report_delay=report_delay,
        analysis_window=dict(),
        reporting_factor=rho,
        Cdates=Cdates
    )

    assert results['est_beta'] is None
    assert results['status'] == 'failed'
    assert results['msg'] == 'Transmission parameter beta not solved (stale).'

def test_multinode_decay_lsq_diff_beta():
    # test with different beta values for each column
    N = 1000000
    y0={'S': N, 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5
    gamma = 1/4
    rho = 8
    report_delay = 7
    tf = 25
    tx = [0]*tf
    tx[10] = 1

    cm_rep_cases = list()
    for beta in [0.5, 0.75, 1.0]:
        Cdates,C,dates,T,S,E,I1,I2,I3,R = simulate_discrete_seiiir_deterministic(y0, tf, beta=beta,
                                                                        sigma=sigma, gamma=gamma,
                                                                        rho=rho, N=N,
                                                                        report_delay=report_delay,
                                                                        tx=tx)
        cm_rep_cases.append(C)

    # create a dataframe with the column of cases
    df_cm_rep_cases = pd.DataFrame({'0': list(cm_rep_cases[0]), '1': list(cm_rep_cases[1]), '3': list(cm_rep_cases[2])})
    df_populations = pd.DataFrame({'county': ['0', '1', '2'], 'population':[N, N, N]})
    df_populations = df_populations.set_index('county')
    results = run_multinode_decay_lsq(
        cm_rep_cases=df_cm_rep_cases,
        populations=df_populations,
        sigma=sigma,
        gamma=gamma,
        report_delay=report_delay,
        analysis_window=dict(),
        reporting_factor=rho,
        Cdates=Cdates
    )

    # this is the result we obtain - beta seems weighted to the column with the most cases
    assert 0.888 == pytest.approx(results['est_beta'],1e-2)
    assert results['status'] == 'ok'
    assert results['msg'] == 'Optimal solution found'
