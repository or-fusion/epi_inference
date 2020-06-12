from datetime import datetime, timedelta
import numpy as np
from pyutilib.misc.misc import Bunch

"""
This module provides some common functions used for reconstruction
computations.
"""
"""
# draft of a decorator that forces named arguments
def force_keyword_args(func):
    def force_keyword_args_wrapper(func, *args, **kwargs):
        if len(args) != 0:
            raise SyntaxError('function {} accepts keyword arguments only'.format(str(func)))

        return func(**kwargs)
    return force_keyword_args_wrapper(func)
"""

def reported_cases_from_cumulative(*, dates, cumulative_reported_cases):
    """
    This function takes the cumulative reported cases (by day) and does
    the difference to return the number of reported cases within each day.
    It requires that the cumulative reported cases start at zero since
    the models in this module assume a fully suspectible population at
    the start.

    Parameters
    ----------
    dates : list of datetime objects
       This is the dates for each of the reported cases using a python datetime object.
       
    cumulative_reported_cases : list of numbers
       This is a list of the cumulative reported cases. The list must include
       all cases (i.e., start at zero) in order to be consistent with other model
       assumptions in the reconstruction module.
    """
    if cumulative_reported_cases[0] != 0:
        # This is necessary since we assume that S0, E0, and I0 are all zero.
        # Therefore, we need to start well before the beginning of the
        # reported cases
        raise ValueError('reported_cases_from_cumulative: The first'
                         ' data point in cumulative_reported_cases must'
                         ' be zero. These procedures assume the data spans'
                         ' the entire timeline.'
                         )

    if len(cumulative_reported_cases) != len(dates):
        raise ValueError("The length of the dates list is not the same as the "
                         "cumulative_reported_cases list")
    
    cumul_rep_cases = cumulative_reported_cases
    reported_cases = [j-i for i,j in zip(cumul_rep_cases, cumul_rep_cases[1:])]
    assert len(reported_cases) == len(cumul_rep_cases)-1
    reported_cases_dates = dates[1:]
    assert len(reported_cases_dates) == len(dates)-1

    # sanity check - the cumulative reported cases should be non-decreasing
    # therefore, we should not have any negative reported cases
    if not all([c >= 0 for c in reported_cases]):
        raise ValueError('reported_cases_from_cumulative: Negative reported'
                         ' cases detected for an interval. cumulative_reported_cases'
                         ' must be a list of non-decreasing numbers.')
    
    return Bunch(dates=reported_cases_dates, values=reported_cases)

def df_reported_cases_from_cumulative(df_cumulative_reported_cases):
    """
    This function takes the cumulative reported cases (by day) as a Pandas
    dataframe where each column is a different county (node) and just
    the difference to return a new dataframe with the number of reported
    cases within each day.

    It requires that the cumulative reported cases start for each county
    start at zero since the models in this module assume a fully suspectible population at
    the start.

    Parameters
    ----------
    df_cumulative_reported_cases : DataFrame
       This is a dataframe that contains cumulative reported cases. Each row
       corresponds to a different day, and each column is a different county (node).
       The index should be the 'date' column.
       The cumulative reported cases must all start at zero in order to be
       consistent with other model assumptions in the reconstruction module.
    """
    first_row = df_cumulative_reported_cases.iloc[0]
    first_row = first_row.drop(['date'])
    if (first_row != 0).any():
        # This is necessary since we assume that S0, E0, and I0 are all zero.
        # Therefore, we need to start well before the beginning of the
        # reported cases
        raise ValueError('df_reported_cases_from_cumulative: The first'
                         ' row for df_cumulative_reported_cases must'
                         ' all be zero. These procedures assume the data spans'
                         ' the entire timeline (from prior to the initial reported'
                         ' case.'
                         '\nFirst row:\n{}'.format(first_row)
                         )

    df_reported_cases = df_cumulative_reported_cases.set_index('date')
    df_reported_cases = df_reported_cases.diff()
    df_reported_cases = df_reported_cases.drop(df_reported_cases.index[0])

    # sanity check - the cumulative reported cases should be non-decreasing
    # therefore, we should not have any negative reported cases
    if (df_reported_cases < 0).any(axis=None):
        raise ValueError('df_reported_cases_from_cumulative: Negative reported'
                         ' cases detected for an interval. Columns in df_cumulative_reported_cases'
                         ' must all be non-decreasing.')

    return df_reported_cases

def np_reported_cases_from_cumulative(*, dates, counties, cumulative_reported_cases):
    """
    This function takes the cumulative reported cases (by day) and returns the
    new cases within each day.

    It requires that the cumulative reported cases start for each county
    start at zero since the models in this module assume a fully suspectible population at
    the start.

    Parameters
    ----------
    dates : numpy array of datetime objects
       The dates corresponding to the rows in the cumulative reported cases
    counties : numpy array of object (strings)
       The names of the counties (or nodes) corresponding to the columns in the
       cumulative reported cases
    cumulative_reported_cases : Numpy two-dimensional array
       This is a numpy array that contains cumulative reported cases. Each row
       corresponds to a different day, and each column is a different county (node).
       The cumulative reported cases must all start at zero in order to be
       consistent with other model assumptions in the reconstruction module.  
    """
    # check the types
    assert isinstance(dates, np.ndarray) and dates.dtype == np.object
    assert isinstance(counties, np.ndarray) and counties.dtype == np.object
    
    # check the sizes of the incoming data
    if len(dates) != cumulative_reported_cases.shape[0] or \
       len(counties) != cumulative_reported_cases.shape[1]:
        raise ValueError('Dimension error in np_reported_cases_from_cumulative.'
                         ' length of dates must be equal to the number of rows in'
                         ' cumulative_reported_cases and length of counties '
                         ' must be equal to the number of columns in'
                         ' cumulative_reported_cases.'
                         )

    if sum(cumulative_reported_cases[0,:]) > 0:
        # This is necessary since we assume that S0, E0, and I0 are all zero.
        # Therefore, we need to start well before the beginning of the
        # reported cases
        raise ValueError('np_reported_cases_from_cumulative: The first'
                         ' row for cumulative_reported_cases must'
                         ' all be zero. These procedures assume the data spans'
                         ' the entire timeline (from prior to the initial reported'
                         ' case.'
                         '\nFirst row:\n{}'.format(list((c,v) for c,v in zip(counties,cumulative_reported_cases[0,:])))
                         )

    np_reported_cases = np.diff(cumulative_reported_cases, axis=0)
    dates = dates[1:]

    # sanity check - the cumulative reported cases should be non-decreasing
    # therefore, we should not have any negative reported cases
    if np.any(np_reported_cases < 0):
        raise ValueError('df_reported_cases_from_cumulative: Negative reported'
                         ' cases detected for an interval. Columns in df_cumulative_reported_cases'
                         ' must all be non-decreasing.')

    return Bunch(dates=dates, counties=counties, values=np_reported_cases)

