__all__ = ['run']

import sys
import datetime
import pandas as pd
from pyutilib.misc import timing

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

    # reconstruct the states
    Cdates = [datetime.date.fromisoformat(day) for day in df.index.to_list()]
    reported_cases_per_day = \
            reported_cases_from_cumulative(dates=Cdates,
                                           cumulative_reported_cases=cumulative_reported_cases)

    res = reconstruct_states_deterministic_decay(dates=reported_cases_per_day.dates,
                        reported_cases_per_day=reported_cases_per_day.values,
                        population=population,
                        sigma=CONFIG['sigma'],
                        gamma=CONFIG['gamma'],
                        reporting_factor=CONFIG['reporting_factor'],
                        report_delay=CONFIG['deltaP'],
                        county=county,
                        warnings=warnings)

    # TODO - keep rdates and rcases?
    #results['rdates'] = reported_cases_per_day.dates
    #results['rcases'] = reported_cases_per_day.values
    results['dates'] = res.dates
    results['transmissions'] = res.transmissions
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
    results = {}
    if 'county' in CONFIG:
        counties = [CONFIG['county']]
    else:
        counties = df.keys()

    if CONFIG['verbose']:
        timing.tic()
    for t in counties:
        if t not in population_df[CONFIG['population_csv']['population']]:
            warnings.append("WARNING: county %s does not have population data available" % str(t))
            continue
        results[t] = run_county(t, df, population_df[CONFIG['population_csv']['population']][t], CONFIG, warnings)
    if CONFIG['verbose']:
        timing.toc("Serial Execution")
    #
    # Save results
    #
    save_results(results, CONFIG['output_json'])
    save_metadata(CONFIG, warnings)


class ReconstructionDeterministicDelay(Task):

    def __init__(self):
        Task.__init__(self, "reconstruction_deterministic_delay",
            "Perform compartment reconstruction using a deterministic delay.")

    def validate(self, CONFIG):
        valid_options = set(['sigma', 'gamma', 'reporting_factor', 'deltaP', 'population_csv', 'input_csv', 'county', 'output_json', 'verbose', 'factors', 'factor_levels', 'workflow'])
        for key in CONFIG:
            if key not in valid_options:
                raise RuntimeError("Unexpected configuration option: '%s'" % key)

    def run(self, data, CONFIG):
        self._warnings = []
        self.validate(CONFIG)
        run(CONFIG, self._warnings)


register_task(ReconstructionDeterministicDelay())

