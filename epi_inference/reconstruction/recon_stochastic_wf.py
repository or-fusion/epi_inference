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
from ..reconstruction.stochastic import stochastic_reconstruction, stochastic_reconstruction_from_daily_transmissions
from ..reconstruction.common import reported_cases_from_cumulative


def run_county(county, reported_cases_df, population, CONFIG, warnings, alternate_transmissions_df=None):
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
    cumulative_reported_cases = reported_cases_df[county].to_list()

    #
    # Get the reported cases per day from the cumulative reported cases
    Cdates = [datetime.date.fromisoformat(day) for day in reported_cases_df.index.to_list()]
    reported_cases_per_day = \
            reported_cases_from_cumulative(dates=Cdates,
                                           cumulative_reported_cases=cumulative_reported_cases)

    #
    # Setup arguments for stochastic_reconstruction
    #
    if 'alternate_transmissions_file' in CONFIG:
        assert alternate_transmissions_df is not None
        # merge the reported cases for the county in with the transmissions

        dates  = [datetime.date.fromisoformat(day) for day in alternate_transmissions_df.index.to_list()]
        args = {'dates':dates,
                'transmissions': alternate_transmissions_df[county].to_list(),
                'population':population,
                'n_steps_per_day':CONFIG['n_steps_per_day']}
        for option in ['fixed_incubation', 'infectious_lower', 'infectious_upper', 'seed']:
            if option in CONFIG:
                args[option] = CONFIG[option]
        res = stochastic_reconstruction_from_daily_transmissions(**args)

        # add the reported cases to the results (need to match the dates - this could be done more efficiently
        date_set =  set(dates)
        # get the reported cases for each day in dates
        reported_cases_dict = {d:v for d,v in zip(reported_cases_per_day.dates, reported_cases_per_day.values) if d in date_set}
        original_reported_cases_for_dates = list()
        found_first = False
        for d in dates:
            # check that all the transmission dates are in reported cases dict
            if d not in reported_cases_dict:
                # reported cases may start later than transmissions
                # we allow this and put zeros on the front, but
                # throw an error if a date is missing after the first
                # coinciding dates
                assert found_first is False
                print(d, 'in dates, but not in reported_cases_dict')
                original_reported_cases_for_dates.append(0)
            else:
                original_reported_cases_for_dates.append(reported_cases_dict[d])
        res.orig_rep_cases = original_reported_cases_for_dates
    else:
        args = {'dates':reported_cases_per_day.dates,
                'reported_cases_per_day':reported_cases_per_day.values,
                'population':population,
                'n_steps_per_day':CONFIG['n_steps_per_day']}
        for option in ['reporting_delay_mean', 'reporting_delay_dev', 'reporting_multiplier', 'fixed_incubation', 'infectious_lower',            'infectious_upper', 'seed']:
            if option in CONFIG:
                args[option] = CONFIG[option]
        res = stochastic_reconstruction(**args)

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
    reported_df = load_collect(CONFIG['input_csv'])
    #
    # Load transmissions file if it exists
    #
    transmissions_file = CONFIG.get('alternate_transmissions_file', None)
    transmissions_df = None
    if transmissions_file is not None:
        transmissions_df = pd.read_csv(transmissions_file, index_col='Date')

    #
    # Perform construction
    #
    results = {}
    if 'county' in CONFIG:
        counties = [CONFIG['county']]
    else:
        counties = reported_df.keys()

    if CONFIG['verbose']:
        timing.tic()
    for t in counties:
        if t not in population_df[CONFIG['population_csv']['population']]:
            warnings.append("WARNING: county %s does not have population data available" % str(t))
            continue
        results[t] = run_county(t, reported_df, population_df[CONFIG['population_csv']['population']][t],
                                CONFIG, warnings, alternate_transmissions_df=transmissions_df)
    if CONFIG['verbose']:
        timing.toc("Serial Execution")
    #
    # Save results
    #
    save_results(results, CONFIG['output_json'])
    save_metadata(CONFIG, warnings)


class ReconstructionStochastic(Task):

    def __init__(self):
        Task.__init__(self, "reconstruction_stochastic",
            "Perform stochastic compartment reconstruction.")

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


register_task(ReconstructionStochastic())

