__all__ = ['run']

import sys
import datetime
import pandas as pd
from pyutilib.misc import timing
try:
    import joblib
    joblib_available = True
except:
    joblib_available = False

from ..engine.task import Task
from ..engine.task_registry import register_task
from ..engine.misc import save_metadata

from ..util import load_population, save_results
from ..collect.misc import load_collect
from ..reconstruction.deterministic import reconstruct_states_deterministic_decay
from ..reconstruction.common import reported_cases_from_cumulative


def run_county(county, df, population, CONFIG, warnings):
    #
    # Initialize results dictionary
    #
    results = {'FIPS':county}
    for key, value in CONFIG.get('factor_levels',{}).items():
        if not key in results:
            results[key] = value
    #
    # Get the cumulative cases
    #
    cumulative_reported_cases = df[county].to_list()
    if df[county][-1] == 0:
        warnings.append( "WARNING: County %s has no reported cases" % str(county))
        return results

    # reconstruct the states
    Cdates = [datetime.date.fromisoformat(day) for day in df.index.to_list()]
    reported_cases_per_day = \
            reported_cases_from_cumulative(dates=Cdates,
                                           cumulative_reported_cases=cumulative_reported_cases)

    res = reconstruct_states_deterministic_decay(dates=reported_cases_per_day.dates,
                        reported_cases_per_day=reported_cases_per_day.values,
                        population=population,
                        sigma=CONFIG['sigma'],
                        gamma=CONFIG['gamma']/3,
                        reporting_factor=CONFIG['reporting_factor'],
                        report_delay=CONFIG['deltaP'])

    # TODO - keep rdates and rcases?
    #results['rdates'] = reported_cases_per_day.dates
    #results['rcases'] = reported_cases_per_day.values
    results['dates'] = res.dates
    results['T'] = res.transmissions
    results['S'] = res.S
    results['E'] = res.E
    results['I1'] = res.I1
    results['I2'] = res.I2
    results['I3'] = res.I3
    results['R'] = res.R
    results['population'] = population

    return results


def run(CONFIG, warnings):
    #
    # Load the population data
    #
    population_df = load_population(CONFIG['population_csv']['file'], CONFIG['population_csv']['index'])
    #
    # Load the case data 
    #
    df = load_collect(CONFIG['input_csv'])
    #
    # Perform construction
    #
    results = []
    if 'county' in CONFIG:
        counties = [CONFIG['county']]
    else:
        counties = list(df.keys())

    parallel = ('parallel' in CONFIG) and (len(counties) >= CONFIG['parallel'].get('number_of_counties',10))

    if parallel and np > 1:
        if CONFIG['verbose']:
            timing.tic()
        available_counties = []
        for t in counties:
            if t not in population_df[CONFIG['population_csv']['population']]:
                warnings.append("WARNING: county %s does not have population data available" % str(t))
                continue
            available_counties.append(t)
        np = CONFIG['parallel'].get('np',2)
        if CONFIG['verbose']:
            timing.toc("Parallel Setup")
        #
        with joblib.Parallel(n_jobs=np) as parallel:
            unordered_results = parallel( joblib.delayed(run_county)(t, df, population_df[CONFIG['population_csv']['population']][t], CONFIG, warnings) for t in available_counties)
        if CONFIG['verbose']:
            timing.toc("Parallel Execution")
        #
        # Order the results
        #
        tmp = {res['FIPS'] : res for res in unordered_results}
        for t in available_counties:
            results.append( tmp[t] )
        if CONFIG['verbose']:
            timing.toc("Reorder Results")
    else:
        if CONFIG['verbose']:
            timing.tic()
        for t in counties:
            if t not in population_df[CONFIG['population_csv']['population']]:
                warnings.append("WARNING: county %s does not have population data available" % str(t))
                continue
            results.append( run_county(t, df, population_df[CONFIG['population_csv']['population']][t], CONFIG, warnings) )
        if CONFIG['verbose']:
            timing.tic("Serial Execution")
    #
    # Save results
    #
    save_results(results, CONFIG['output_json'])
    save_metadata(CONFIG, warnings)


class ParallelReconstructionDeterministicDelay(Task):

    def __init__(self):
        Task.__init__(self, "parallel_reconstruction_deterministic_delay",
            "Perform compartment reconstruction using a deterministic delay.")

    def validate(self, args):
        pass

    def run(self, data, CONFIG):
        if not joblib_available:
            raise RuntimeError("ERROR: Cannot execute the parallel_reconstruction_determinstic_delay workflow.  The 'joblib' package is missing.")
        self._warnings = []
        run(CONFIG, self._warnings)

    def warnings(self):
        return self._warnings


register_task(ParallelReconstructionDeterministicDelay())

