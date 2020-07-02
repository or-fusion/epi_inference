__all__ = ['run']

import sys
import datetime
#import pandas as pd
from pyutilib.misc import timing

import rpy2.robjects as robjects

from ..engine.task import Task
from ..engine.task_registry import register_task
from ..engine.misc import save_metadata

from ..util import load_population, save_results
from ..collect.misc import load_collect
from ..reconstruction.resampled import resampled_reconstruction, sample_county_negbin
from ..reconstruction.common import reported_cases_from_cumulative


def run_county(county, df, population, CONFIG, warnings, cache):
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
    cumulative_reported_cases = sample_county_negbin(dat=df, county=county, parameters=cache)

    # reconstruct the states
    Cdates = [datetime.date.fromisoformat(day) for day in df.index.to_list()]
    Cdates = Cdates[(len(Cdates)-len(cumulative_reported_cases)):]
    #print("HERE",df.shape[0], cumulative_reported_cases.shape[0], len(Cdates))
    reported_cases_per_day = \
            reported_cases_from_cumulative(dates=Cdates,
                                           cumulative_reported_cases=cumulative_reported_cases)

    # Setup arguments for stochastic_reconstruction
    args = {'dates':reported_cases_per_day.dates,
            'reported_cases_per_day':reported_cases_per_day.values,
            'population':population,
            'n_steps_per_day':CONFIG['n_steps_per_day']}
    for option in ['reporting_delay_mean', 'reporting_delay_dev', 'reporting_multiplier', 'fixed_incubation', 'infectious_lower', 'infectious_upper', 'seed']:
        if option in CONFIG:
            args[option] = CONFIG[option]
    res = resampled_reconstruction(**args)

    results['dates'] = res.dates
    results['transmissions'] = res.transmissions
    results['S'] = res.S
    results['E'] = res.E
    results['I1'] = res.I1
    results['I2'] = res.I2
    results['I3'] = res.I3
    results['R'] = res.R
    results['population'] = population
    results['orig_rep_cases'] = res.orig_rep_cases
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
    results = {}
    if 'county' in CONFIG:
        counties = [CONFIG['county']]
    else:
        counties = df.keys()

    if 'seed' in CONFIG:
        robjects.r('set.seed(%s)' % str(CONFIG['seed']))
    else:
        robjects.r('set.seed(123456789)')
    if CONFIG['verbose']:
        timing.tic()
    cache={}
    for t in sorted(counties):
        if t not in population_df[CONFIG['population_csv']['population']]:
            warnings.append("WARNING: county %s does not have population data available" % str(t))
            continue
        results[t] = run_county(t, df, population_df[CONFIG['population_csv']['population']][t], CONFIG, warnings, cache)
    if CONFIG['verbose']:
        timing.toc("Serial Execution")
    #
    # Save results
    #
    save_results(results, CONFIG['output_json'])
    save_metadata(CONFIG, warnings)


class ReconstructionResampled(Task):

    def __init__(self):
        Task.__init__(self, "reconstruction_resampled",
            "Perform compartment reconstruction with resampled data.")

    def validate(self, CONFIG):
        valid_options = set(['reporting_delay_mean', 'reporting_delay_dev', 'reporting_multiplier', 'fixed_incubation', 'infectious_lower', 'infectious_upper', 'seed', 'n_steps_per_day', 'population_csv', 'input_csv', 'county', 'output_json', 'verbose', 'factors', 'factor_levels', 'workflow'])
        for key in CONFIG:
            if key not in valid_options:
                raise RuntimeError("Unexpected configuration option: '%s'" % key)

    def run(self, data, CONFIG):
        self._warnings = []
        run(CONFIG, self._warnings)

    def warnings(self):
        return self._warnings


register_task(ReconstructionResampled())

