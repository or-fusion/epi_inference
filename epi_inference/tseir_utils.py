import numpy as np
import pyomo.environ as pe
import pyomo.dae as dae
import math

def reported_cases_from_cumulative(cumulative_reported_cases):
    cumul_rep_cases = [0] + cumulative_reported_cases
    reported_cases = [cumul_rep_cases[i] - cumul_rep_cases[i-1] for i in range(1,len(cumul_rep_cases))]
    return reported_cases

def compute_compartments_time_delay(cumulative_reported_cases, population, deltaE, deltaI, deltaP, reporting_factor):   # pragma: no cover
    print("WARNING - THIS CODE IS NOT TESTED")
    reported_cases = reported_cases_from_cumulative(cumulative_reported_cases)
    cases = [reporting_factor*c for c in reported_cases]

    S = [None]*len(cases)
    S[0] = population
    E = [0]*len(cases)
    I = [0]*len(cases)
    R = [0]*len(cases)
    for t in range(len(cases)):
        if t+1 < len(cases):
            S[t+1] = S[t] - cases[t]
            if t - deltaE - deltaI >= 0:
                R[t+1] = R[t] + cases[t-deltaE-deltaI]
        E[t] = sum(cases[t-tau] for tau in range(1,deltaE+1) if t-tau >= 0)
        I[t] = sum(cases[t-deltaE-tau] for tau in range(1,deltaI+1) if t-deltaE-tau >= 0)
    return cases, S, E, I, R

def compute_compartments_decay(cumulative_reported_cases, population, sigma, gamma_1, gamma_2, gamma_3, deltaP, reporting_factor):
    reported_cases = reported_cases_from_cumulative(cumulative_reported_cases)
    cases = [reporting_factor*c for c in reported_cases]

    S = [None]*len(cases)
    S[0] = population
    E = [0]*len(cases)
    I1 = [0]*len(cases)
    I2 = [0]*len(cases)
    I3 = [0]*len(cases)
    R = [0]*len(cases)
    for t in range(len(cases)-1):
        S[t+1] = S[t] - cases[t]
        E[t+1] = E[t] + cases[t] - sigma*E[t]
        I1[t+1] = I1[t] + sigma*E[t] - gamma_1*I1[t]
        I2[t+1] = I2[t] + gamma_1*I1[t] - gamma_2*I2[t]
        I3[t+1] = I3[t] + gamma_2*I2[t] - gamma_3*I3[t]
        R[t+1] = R[t] + gamma_3*I3[t]
    #I = [I1[t] + I2[t] + I3[t] for t in len(cases)]
    return cases, S, E, I1, I2, I3, R
