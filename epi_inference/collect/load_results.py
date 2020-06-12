
import sys
import os
import pandas as pd
import glob
import importlib


"""
TODO - delete?


prefix="around_md_sim_"

def load_results(expdir=None, expnum=1, trial=0, county="24031", threshold=0, start=0, stop=None, days_before_first=None, days_after_first=None, verbose=False):
    if verbose:
        print('load_results:', expdir, 'exp'+str(expnum), prefix+str(trial)+"_scn_"+county+".csv")
    filename=os.path.join(expdir, 'exp'+str(expnum), prefix+str(trial)+"_scn_"+county+".csv")
    if not os.path.exists(filename):
        print("MISSING FILE: "+filename)
        print("  Experiment "+str(expnum))
        print("  Trial "+str(trial))
        print("  County "+county)
        return [], [], None

    df = pd.read_csv(filename) 
    tmp = df.loc[ df['cumI'] >= threshold, ['cumI'] ]
    cumI = tmp['cumI'].tolist()
    index = list(range(len(cumI)))

    #
    # TODO - raise warning if both start/stop and days_after_first are specified
    #
    try:
        istart = cumI.index(next(filter(lambda x: x>0, cumI)))
        if not days_after_first is None:
            stop = min(len(cumI), istart + days_after_first)
        if not days_before_first is None:
            start = max(0, istart - days_before_first)
    except StopIteration:
        print("WARNING: no infections reported in the data.  Ignoring days_after_first and days_before_first options.")
    if stop is None or stop >= len(cumI):
        if start > 0:
            cumI = cumI[start:]
            index = index[start:]
    else:
        if start > 0:
            cumI = cumI[start:stop]
            index = index[start:stop]
        else:
            cumI = cumI[:stop]
            index = index[:stop]

    filename=os.path.join(expdir, 'exp'+str(expnum), "README.txt")
    with open(filename, "r") as INPUT:
        INPUT.readline()
        INPUT.readline()
        INPUT.readline()
        line = INPUT.readline()
        beta = line.strip().split(" ")[2]
        line = INPUT.readline()
        line.strip()
        gamma = line.strip().split(" ")[1]
        line = INPUT.readline()
        line.strip()
        sigma = line.strip().split(" ")[1]
    return index, cumI, {'beta':beta, 'gamma':gamma, 'sigma':sigma}
"""

#
# Load data files into a pandas dataframe
#
def load_data(files, index_col=0, names=None, dtype=None, header=0, labels=None):
    if labels is None:
        labels = list(range(len(files)))
    if len(files) == 1:
        DF = pd.read_csv(files[0], header=header, dtype={"FIPS":'str'}, encoding="ISO-8859-1")
        DF = DF.set_index(index_col)
        df = DF[names]
        df.index = df.index.rename('Date')
        if dtype == 'casedata':
            label = str(DF["FIPS"][0])
            if len(label) < 5:
                label = "0"*(5-len(label))+label
            df = df.rename(columns={names[0]:label})
        else:
            label = str(labels[0])
        return df.rename(columns={names[0]:label})
    else:
        dfs = []
        i = 0
        for fname in files:
            DF = pd.read_csv(fname, header=header, dtype={"FIPS":'str'}, encoding="ISO-8859-1")
            DF = DF.set_index(index_col)
            df = DF[names]
            if dtype == "casedata":
                label = str(DF["FIPS"][0])
                if len(label) < 5:
                    label = "0"*(5-len(label))+label
            else:
                label = str(labels[i])
            df = df.rename(columns={names[0]:label})
            df.index = df.index.rename('Date')
            dfs.append(df)
            i = i+1
        dfall = pd.concat(dfs, axis=1)
        dfall.index = dfall.index.rename('Date')
        return dfall
         

#
# Load data for one or more trials into a dataframe
#
def load_df_expdata(expdir=None, expnum=1, trial=None, county="24031", start=0, stop=None, days_before_first=None, days_after_first=None):
    #
    # load data
    #
    countycsvfiles = list(glob.glob(os.path.join(expdir, 'exp'+str(expnum), "*_"+county+"*.csv")))
    countycsvfiles.sort()
    if len(countycsvfiles) == 0:                    # pragma: no cover
        print("ERROR: no experimental files available for county '%s' in directory '%s'" % (county, os.path.join(expdir, 'exp'+str(expnum))))
        return None, None
    if trial is None:
        df = load_data(countycsvfiles, dtype='seir', index_col='time', names=['cumI'])
    else:
        try:
            files = [countycsvfiles[trial]]
        except Exception as exc:                    # pragma: no cover
            print("ERROR: bad trial id '%s'" % str(trial))
            print(exc)
            return None, None
        df = load_data(files, dtype='seir', index_col='time', names=['cumI'], labels=[trial])
    #print(df.head())
    #
    # filter data
    #
    cumI = df.sum(axis=1)
    #print(cumI)
    cumI = cumI.tolist()

    #
    # TODO - raise warning if both start/stop and days_after_first are specified
    #
    try:
        istart = cumI.index(next(filter(lambda x: x>0, cumI)))
        if not days_after_first is None:
            stop = min(len(cumI), istart + days_after_first+1)
        if not days_before_first is None:
            start = max(0, istart - days_before_first)
    except StopIteration:                           # pragma: no cover
        print("WARNING: no infections reported in the data.  Ignoring days_after_first and days_before_first options.")
    if stop is None or stop >= len(cumI):
        if start > 0:
            df = df.iloc[start:,]
    else:
        if start > 0:
            df = df.iloc[start:stop,]
        else:
            df = df.iloc[:stop,]

    filename=os.path.join(expdir, 'exp'+str(expnum), "README.txt")
    if os.path.exists(filename):
        with open(filename, "r") as INPUT:
            INPUT.readline()
            INPUT.readline()
            INPUT.readline()
            line = INPUT.readline()
            beta = line.strip().split(" ")[2]
            line = INPUT.readline()
            line.strip()
            gamma = line.strip().split(" ")[1]
            line = INPUT.readline()
            line.strip()
            sigma = line.strip().split(" ")[1]
            R0 = str(float(beta)/float(gamma))
    else:
        modulepath=os.path.join(expdir, 'exp'+str(expnum))
        sys.path.insert(0, modulepath)
        importlib.invalidate_caches()
        module = importlib.import_module('info')
        module = importlib.reload(module)
        sys.path.pop(0)
        beta = str(module.beta)
        gamma = str(module.gamma)
        R0 = str(module.R0)
        sigma = str(module.sigma)
    return df, {'beta':float(beta), 'R0':float(R0), 'gamma':float(gamma), 'sigma':float(sigma)}


#
# Load data for one or more trials into a dataframe
#
def load_df_casedata(files, datadir=None, datacol=None, start=0, stop=None, days_before_first=None, days_after_first=None):
    if datacol is None:
        datacol='Confirmed'
    #
    # load data
    #
    files.sort()
    if len(files) == 0:                         # pragma: no cover
        print("ERROR: no case data files were specified")
        return None
    files_with_path = []
    for f in files:
        fname = os.path.join(datadir, f)
        if not os.path.exists(fname):           # pragma: no cover
            print("ERROR: missing file '%s' in directory '%s'" % (f, datadir))
            return None
        files_with_path.append(fname)
    df = load_data(files_with_path, dtype='casedata', index_col="Date", names=[datacol])
    #
    # filter data
    #
    cumI = df.sum(axis=1)
    #print(cumI)
    cumI = cumI.tolist()

    #
    # TODO - raise warning if both start/stop and days_after_first are specified
    #
    try:
        istart = cumI.index(next(filter(lambda x: x>0, cumI)))
        if not days_after_first is None:
            stop = min(len(cumI), istart + days_after_first+1)
        if not days_before_first is None:
            start = max(0, istart - days_before_first)
    except StopIteration:                                   # pragma: no cover
        print("WARNING: no infections reported in the data.  Ignoring days_after_first and days_before_first options.")
    if stop is None or stop >= len(cumI):
        if start > 0:
            df = df.iloc[start:,]
    else:
        if start > 0:
            df = df.iloc[start:stop,]
        else:
            df = df.iloc[:stop,]
    return df

