__all__ = ['load_collect', 'save_collect']

import os
import csv
import pandas as pd


def load_collect(input_csv):
    """
    This function loads a CSV file containing results from the 'collect' workflow.  The
    CSV data is assumed to have a 'Date' column.

    Parameters
    ----------
    input_csv : string
        File name of the CSV file

    Returns
    -------
    Pandas dataframe containing the data from the CSV file.
    """
    try:
        assert(os.path.exists(input_csv))
    except:                                                                 # pragma: no cover
        raise RuntimeError("load_collect: ERROR. Input file "+input_csv+" does not exist.")
    return pd.read_csv(input_csv, index_col='Date')


def save_collect(output_csv, df, verbose, warnings):
    """
    This function loads a CSV file containing results from the 'collect' workflow.  The
    CSV data is assumed to have a 'Date' column.

    Parameters
    ----------
    input_csv : string
        File name of the CSV file

    Returns
    -------
    Pandas dataframe containing the data from the CSV file.
    """
    filedir = os.path.dirname(output_csv)
    if filedir and not os.path.exists(filedir):             # pragma: no cover
        os.makedirs(filedir)
    if verbose and os.path.exists(output_csv):              # pragma: no cover
        warnings.append( "WARNING: over-writing file "+output_csv )
    print("Writing file: "+output_csv)
    df.to_csv(output_csv, quoting=csv.QUOTE_NONNUMERIC)

