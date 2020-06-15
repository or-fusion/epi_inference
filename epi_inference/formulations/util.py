__all__ = ['get_windows']

import datetime
from pyutilib.misc import Options

def get_windows(dates, last_day_of_window=5, window_days=7, select_window=None):
    ans = Options()
    ans.TIMES = [i for i in range(len(dates))]
    ans.WINDOW_TIMES_LIST = []
    ans.WINDOWS = []
    ans.WINDOW_TIMES = {}

    for i in ans.TIMES:
        if datetime.date.fromisoformat(dates[i]).weekday() != last_day_of_window:
            continue
        #if i < window_days-1:
        if i < window_days:
            continue

        ans.WINDOW_TIMES[i] = []
        for j in range(i+1-window_days, i+1):
            ans.WINDOW_TIMES_LIST.append((i,j)) 
            ans.WINDOW_TIMES[i].append(j) 
        ans.WINDOWS.append(i)

    if len(ans.WINDOWS) == 0:
        return ans

    if select_window is not None:
        if select_window == 'last':
            select_val = ans.WINDOWS[-1]
            ans.WINDOWS = [select_val]
            ans.WINDOW_TIMES_LIST = [(i,j) for i,j in ans.WINDOW_TIMES_LIST if i==select_val]
            ans.WINDOW_TIMES = {select_val:ans.WINDOW_TIMES[select_val]}
        else:
            select_val=None
            for i in ans.WINDOWS:
                if datetime.date.fromisoformat(dates[i]) == datetime.date.fromisoformat(select_window):
                    select_val=i
                    break
            if select_val == None:
                ans.WINDOWS = []
                ans.WINDOW_TIMES_LIST = []
                ans.WINDOW_TIMES = {}
            else:
                ans.WINDOWS = [select_val]
                ans.WINDOW_TIMES_LIST = [(i,j) for i,j in ans.WINDOW_TIMES_LIST if i==select_val]
                ans.WINDOW_TIMES = {select_val:ans.WINDOW_TIMES[select_val]}

    return ans

