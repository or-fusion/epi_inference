import pytest
import numpy as np
from epi_inference.simulation.simulation import simulate_discrete_seiiir_deterministic
from epi_inference.util import roundall

def test_simulate_discrete_seiiir_deterministic():
    """
    Test a short simulation with the seiiir deterministic model
    """
    N = 1000000
    y0={'S': N, 'E': 0, 'I1': 1, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5
    gamma = 1/4
    beta = 1
    rho = 1
    report_delay = 7
    tf=10

    results  = simulate_discrete_seiiir_deterministic(y0, tf, beta=beta,
                                                      sigma=sigma, gamma=gamma,
                                                      rho=rho, N=N,
                                                      report_delay=report_delay,
                                                      tx=None)

    C = results.cumulative_reported_cases.values
    Cdates = results.cumulative_reported_cases.dates
    SEIIIR = results.SEIIIR
    
    assert len(Cdates) == len(SEIIIR.dates)+1
    delta = Cdates[1] - SEIIIR.dates[0]
    assert delta.days == report_delay
    assert SEIIIR.S[0] == 1000000
    assert SEIIIR.E[0] == 0
    assert SEIIIR.I1[0] == 1
    assert SEIIIR.I2[0] == 0
    assert SEIIIR.I3[0] == 0
    assert SEIIIR.R[0] == 0
    assert SEIIIR.transmissions[0] == 1
    assert C[0] == 0

    assert SEIIIR.S[1] == 1000000-1
    assert SEIIIR.E[1] == 1
    assert SEIIIR.I1[1] == pytest.approx(0.25)
    assert SEIIIR.I2[1] == pytest.approx(0.75)
    assert SEIIIR.I3[1] == 0
    assert SEIIIR.R[1] == 0
    assert SEIIIR.transmissions[1] == pytest.approx(1.0*999999/1000000)

    expS = [1000000, 999999.0, 999998.000001, 999996.8000034, 999995.661882242, 999994.3121701872, 999992.5550166771, 999990.288982791, 999987.4325757417, 999983.8563387465]
    expE = [0, 1.0, 1.799999, 2.6399968000011995, 3.2501185980054688, 3.949806933169381, 4.916999056689872, 6.199633131502013, 7.8161135545123095, 9.829127838741329]
    expI1 = [1, 0.25, 0.2625, 0.4256247999999999, 0.63440556000024, 0.808625109601154, 0.9921176640341648, 1.2314292273465157, 1.5477839331370316, 1.9501686941867198]
    expI2 = [0, 0.75, 0.375, 0.290625, 0.3918748499999999, 0.5737728825001799, 0.7499120528259104, 0.9315662612321013, 1.1564634858179121, 1.4499538213072518]
    expI3 = [0, 0.0, 0.5625, 0.421875, 0.32343750000000004, 0.37476551249999995, 0.5240210400001349, 0.6934392996194666, 0.8720345208289425, 1.0853562445706697]
    expR = [0, 0.0, 0.0, 0.421875, 0.73828125, 0.9808593750000001, 1.261933509375, 1.654949289375101, 2.175028764089701, 2.829054654711408]
    expT = [1.0, 0.999999, 1.1999976000011998, 1.1381211580045094, 1.349712054765006, 1.7571535101543665, 2.2660338861501166, 2.856407049310699, 3.5762369951314814, 4.485406348014979]
    expC = [0, 1.0, 1.9999989999999999, 3.1999966000011995, 4.338117758005708, 5.6878298127707145, 7.444983322925081, 9.711017209075198, 12.567424258385897, 16.14366125351738, 20.62906760153236]

    np.testing.assert_allclose(np.asarray(expS), np.asarray(SEIIIR.S))
    np.testing.assert_allclose(np.asarray(expE), np.asarray(SEIIIR.E))
    np.testing.assert_allclose(np.asarray(expI1), np.asarray(SEIIIR.I1))
    np.testing.assert_allclose(np.asarray(expI2), np.asarray(SEIIIR.I2))
    np.testing.assert_allclose(np.asarray(expI3), np.asarray(SEIIIR.I3))
    np.testing.assert_allclose(np.asarray(expR), np.asarray(SEIIIR.R))
    np.testing.assert_allclose(np.asarray(expT), np.asarray(SEIIIR.transmissions))
    np.testing.assert_allclose(np.asarray(expC), np.asarray(results.cumulative_reported_cases.values))

    # test the reporting factor
    rho = 8
    results = simulate_discrete_seiiir_deterministic(y0, tf, beta, sigma, gamma, rho, N, report_delay, tx=None)
    SEIIIR = results.SEIIIR
    np.testing.assert_allclose(np.asarray(expS), np.asarray(SEIIIR.S))
    np.testing.assert_allclose(np.asarray(expE), np.asarray(SEIIIR.E))
    np.testing.assert_allclose(np.asarray(expI1), np.asarray(SEIIIR.I1))
    np.testing.assert_allclose(np.asarray(expI2), np.asarray(SEIIIR.I2))
    np.testing.assert_allclose(np.asarray(expI3), np.asarray(SEIIIR.I3))
    np.testing.assert_allclose(np.asarray(expR), np.asarray(SEIIIR.R))
    np.testing.assert_allclose(np.asarray(expT), np.asarray(SEIIIR.transmissions))
    np.testing.assert_allclose(1/8*np.asarray(expC), np.asarray(results.cumulative_reported_cases.values))

def test_simulate_discrete_seiiir_deterministic_tx():
    """
    Test a short simulation with the seiiir deterministic model including
    transmission from outside
    """
    N = 1000000
    y0={'S': N, 'E': 0, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5
    gamma = 1/4
    beta = 1
    rho = 1
    report_delay = 7
    tf=30
    tx = [0]*tf
    tx[10] = 1
    
    results = simulate_discrete_seiiir_deterministic(y0, tf, beta=beta,
                                                     sigma=sigma, gamma=gamma,
                                                     rho=rho, N=N,
                                                     report_delay=report_delay,
                                                     tx=tx)
    Cdates = results.cumulative_reported_cases.dates
    C = results.cumulative_reported_cases.values
    SEIIIR = results.SEIIIR

    assert len(Cdates) == len(SEIIIR.dates)+1
    delta = Cdates[1] - SEIIIR.dates[0]
    assert delta.days == report_delay
    assert SEIIIR.S[0] == 1000000
    assert SEIIIR.E[0] == 0
    assert SEIIIR.I1[0] == 0
    assert SEIIIR.I2[0] == 0
    assert SEIIIR.I3[0] == 0
    assert SEIIIR.R[0] == 0
    assert SEIIIR.transmissions[0] == 0
    assert C[0] == 0

    assert SEIIIR.transmissions[10] == 1
    assert SEIIIR.S[11] == 1000000-1
    assert SEIIIR.E[11] == 1
    assert SEIIIR.I1[11] == 0
    assert SEIIIR.I2[11] == 0
    assert SEIIIR.I3[11] == 0
    assert SEIIIR.R[11] == 0


    assert SEIIIR.transmissions[11] == 0
    assert SEIIIR.S[12] == 1000000-1
    assert SEIIIR.E[12] == pytest.approx(0.8)
    assert SEIIIR.I1[12] == pytest.approx(0.2)
    assert SEIIIR.I2[12] == 0
    assert SEIIIR.I3[12] == 0
    assert SEIIIR.R[12] == 0

    t = 0.2*999999/1000000
    assert SEIIIR.transmissions[12] == pytest.approx(t)
    assert SEIIIR.S[13] == pytest.approx(1000000 - 1 - t)
    assert SEIIIR.E[13] == pytest.approx(0.8 - 0.8*1/5 + 0.2)
    assert SEIIIR.I1[13] == pytest.approx(0.2 + 0.8*1/5 - 0.2*3*1/4)
    assert SEIIIR.I2[13] == pytest.approx(0.2*3*1/4)
    assert SEIIIR.I3[13] == 0
    assert SEIIIR.R[13] == 0

    expT = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.1999998, 0.35999956800007205, 0.5279991363203962, 0.6500234843491176, 0.7899611612071538, 0.9833997294942862, 1.2399269425355373, 1.5632238454281848, 1.9658282341630997, 2.4697471508748774, 3.103237816906961, 3.9001484736462553, 4.902027026250278, 6.161128231301663, 7.743451993146161, 9.732090520466189, 12.231436853410257, 15.372620496106643]
    expS = [1000000, 1000000.0, 1000000.0, 1000000.0, 1000000.0, 1000000.0, 1000000.0, 1000000.0, 1000000.0, 1000000.0, 1000000.0, 999999.0, 999999.0, 999998.8000002, 999998.4400006321, 999997.9120014957, 999997.2619780113, 999996.4720168501, 999995.4886171207, 999994.2486901782, 999992.6854663327, 999990.7196380985, 999988.2498909476, 999985.1466531307, 999981.246504657, 999976.3444776307, 999970.1833493994, 999962.4398974063, 999952.7078068858, 999940.4763700324]
    expE = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0.8, 0.8399998000000001, 1.0319994080000723, 1.353598662720454, 1.7329024145254808, 2.1762830928275383, 2.724426203756317, 3.419467905540591, 4.298798169860658, 5.404866770051626, 6.793640566916178, 8.538150270439903, 10.730668689998177, 13.48656197824882, 16.95037781390072, 21.303754244266738, 26.77509391587958, 33.651511986113924]
    expI1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2, 0.21000000000000002, 0.22049996, 0.2615248716000145, 0.3361009504440945, 0.4306057205161198, 0.5429080486945377, 0.6806122529248979, 0.8540466443393429, 1.0732712950569674, 1.349291177774567, 1.6960509078268777, 2.1316427810447003, 2.6790444332608105, 3.367073503964967, 4.2318439387713855, 5.3187118335461925, 6.6846967415624645]
    expI2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.15000000000000002, 0.19500000000000003, 0.21412497000000005, 0.2496748962000109, 0.3144944368830736, 0.4015778996078582, 0.5075755114228678, 0.6373530675493905, 0.7998732501418546, 1.0049217838281892, 1.2631988292879723, 1.5878378881921515, 1.9956915578315628, 2.508206214403499, 3.1523566815746, 3.9619721244721893, 4.979526906277692]
    expI3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.11250000000000002, 0.17437500000000003, 0.20418747750000005, 0.2383030415250082, 0.2954465880435573, 0.375045071716783, 0.4744429014963466, 0.5966255260361295, 0.7490613191154234, 0.9409566676499976, 1.1826382888784785, 1.4865379883637333, 1.8684031654646056, 2.3482554521687757, 2.9513313742231433, 3.709311936909928]
    expR = [0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.084375, 0.21515625000000002, 0.368296858125, 0.5470241392687561, 0.768609080301424, 1.0498928840890112, 1.405725060211271, 1.853194204738368, 2.4149901940749356, 3.1207076948124337, 4.007686411471292, 5.122589902744092, 6.523892276842547, 8.28508386596913, 10.498582396636486]
    expC = [0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.1999998, 1.559999368000072, 2.0879985043204683, 2.738021988669586, 3.5279831498767398, 4.511382879371026, 5.751309821906563, 7.314533667334748, 9.280361901497848, 11.750109052372725, 14.853346869279687, 18.75349534292594, 23.655522369176218, 29.81665060047788, 37.560102593624045, 47.29219311409023, 59.52362996750049, 74.89625046360713]

    np.testing.assert_allclose(np.asarray(expT), np.asarray(SEIIIR.transmissions))
    np.testing.assert_allclose(np.asarray(expS), np.asarray(SEIIIR.S))
    np.testing.assert_allclose(np.asarray(expE), np.asarray(SEIIIR.E))
    np.testing.assert_allclose(np.asarray(expI1), np.asarray(SEIIIR.I1))
    np.testing.assert_allclose(np.asarray(expI2), np.asarray(SEIIIR.I2))
    np.testing.assert_allclose(np.asarray(expI3), np.asarray(SEIIIR.I3))
    np.testing.assert_allclose(np.asarray(expR), np.asarray(SEIIIR.R))
    np.testing.assert_allclose(np.asarray(expC), np.asarray(C))

    # test the reporting factor
    rho = 8
    results = simulate_discrete_seiiir_deterministic(y0, tf, beta, sigma, gamma, rho, N, report_delay, tx=tx)
    SEIIIR = results.SEIIIR
    np.testing.assert_allclose(np.asarray(expT), np.asarray(SEIIIR.transmissions))
    np.testing.assert_allclose(np.asarray(expS), np.asarray(SEIIIR.S))
    np.testing.assert_allclose(np.asarray(expE), np.asarray(SEIIIR.E))
    np.testing.assert_allclose(np.asarray(expI1), np.asarray(SEIIIR.I1))
    np.testing.assert_allclose(np.asarray(expI2), np.asarray(SEIIIR.I2))
    np.testing.assert_allclose(np.asarray(expI3), np.asarray(SEIIIR.I3))
    np.testing.assert_allclose(np.asarray(expR), np.asarray(SEIIIR.R))
    np.testing.assert_allclose(1/8*np.asarray(expC), np.asarray(results.cumulative_reported_cases.values))

def test_simulate_discrete_seiiir_deterministic_larger():
    """
    Test a longer simulation with the seiiir deterministic model and verify
    the end conditions
    """
    N = 1000
    y0={'S': N, 'E': 5, 'I1': 0, 'I2': 0, 'I3':0, 'R': 0}
    sigma = 1/5.2
    gamma = 1/4
    beta = 2.2*gamma
    rho = 1
    report_delay = 7
    tf=1000
    
    results = simulate_discrete_seiiir_deterministic(y0, tf, beta, sigma, gamma, rho, N, report_delay, tx=None)
    C = results.cumulative_reported_cases.values
    T = results.SEIIIR.transmissions
    S = results.SEIIIR.S
    E = results.SEIIIR.E
    I1 = results.SEIIIR.I1
    I2 = results.SEIIIR.I2
    I3 = results.SEIIIR.I3
    R = results.SEIIIR.R
    assert sum(T) == results.cumulative_reported_cases.values[-1]
    C,T,S,E,I1,I2,I3,R = roundall(C,T,S,E,I1,I2,I3,R)
    
    assert C[-1] == R[-1]-5
    assert S[-1] == N - R[-1]+5
    assert E[-1] == 0
    assert I1[-1] == 0
    assert I2[-1] == 0
    assert I3[-1] == 0
