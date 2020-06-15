import pytest
import os
import os.path
import shutil
import datetime

#from pyomo.common import fileutils as fileutils
#from pyutilib.misc import Options as Options
from epi_inference.formulations.util import get_windows

class TestWindow():

    def test1(self):
        # Default values
        last = datetime.date.fromisoformat('2020-06-14')
        dates=[str(last+datetime.timedelta(days=i-19)) for i in range(20)]
        ans = get_windows(dates)

        assert(ans.TIMES == list(range(20)))
        assert(ans.WINDOWS == [11,18])
        wt = {11: [5, 6, 7, 8, 9, 10, 11], 18: [12, 13, 14, 15, 16, 17, 18]}
        assert(ans.WINDOW_TIMES == wt)
        assert(ans.WINDOW_TIMES_LIST == [(i,j) for i in wt for j in wt[i]])

    def test2(self):
        # Selected a date that wasn't a valid window
        last = datetime.date.fromisoformat('2020-06-14')
        last = datetime.date.fromisoformat('2020-06-14')
        dates=[str(last+datetime.timedelta(days=i-19)) for i in range(20)]
        ans = get_windows(dates, select_window='2020-06-14')

        assert(ans.TIMES == list(range(20)))
        assert(ans.WINDOWS == [])
        wt = {}
        assert(ans.WINDOW_TIMES == wt)

    def test3(self):
        # Selected a valid date 
        last = datetime.date.fromisoformat('2020-06-14')
        last = datetime.date.fromisoformat('2020-06-14')
        dates=[str(last+datetime.timedelta(days=i-19)) for i in range(20)]
        ans = get_windows(dates, select_window='2020-06-13')

        assert(ans.TIMES == list(range(20)))
        assert(ans.WINDOWS == [18])
        wt = {18: [12, 13, 14, 15, 16, 17, 18]}
        assert(ans.WINDOW_TIMES == wt)

    def Xtest4_fails(self):
        # Shouldn't this pass?
        last = datetime.date.fromisoformat('2020-06-14')
        last = datetime.date.fromisoformat('2020-06-14')
        dates=[str(last+datetime.timedelta(days=i-19)) for i in range(20)]
        ans = get_windows(dates, window_days=5)

        assert(ans.TIMES == list(range(20)))
        assert(ans.WINDOWS == [4, 11,18])
        wt = {4: [0,1,2,3,4], 11: [7, 8, 9, 10, 11], 18: [14, 15, 16, 17, 18]}
        assert(ans.WINDOW_TIMES == wt)
        assert(ans.WINDOW_TIMES_LIST == [(i,j) for i in wt for j in wt[i]])

    def test4_succeeds(self):
        # Why isn't the first window included?
        last = datetime.date.fromisoformat('2020-06-14')
        last = datetime.date.fromisoformat('2020-06-14')
        dates=[str(last+datetime.timedelta(days=i-19)) for i in range(20)]
        ans = get_windows(dates, window_days=5)

        assert(ans.TIMES == list(range(20)))
        assert(ans.WINDOWS == [11,18])
        wt = {11: [7, 8, 9, 10, 11], 18: [14, 15, 16, 17, 18]}
        assert(ans.WINDOW_TIMES == wt)
        assert(ans.WINDOW_TIMES_LIST == [(i,j) for i in wt for j in wt[i]])


    def test5(self):
        # Shouldn't this pass?
        last = datetime.date.fromisoformat('2020-06-14')
        last = datetime.date.fromisoformat('2020-06-14')
        dates=[str(last+datetime.timedelta(days=i-29)) for i in range(30)]
        ans = get_windows(dates, window_days=14)

        assert(ans.TIMES == list(range(30)))
        assert(ans.WINDOWS == [14, 21, 28])
        wt = {14: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
              21: [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21],
              28: [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28]}
        assert(ans.WINDOW_TIMES == wt)
        assert(ans.WINDOW_TIMES_LIST == [(i,j) for i in wt for j in wt[i]])

    def test6(self):
        # Default values
        last = datetime.date.fromisoformat('2020-06-14')
        dates=[str(last+datetime.timedelta(days=i-19)) for i in range(20)]
        ans = get_windows(dates, last_day_of_window=6)  # Monday - Sunday

        assert(ans.TIMES == list(range(20)))
        assert(ans.WINDOWS == [12,19])
        wt = {12: [6, 7, 8, 9, 10, 11, 12], 19: [13, 14, 15, 16, 17, 18, 19]}
        assert(ans.WINDOW_TIMES == wt)
        assert(ans.WINDOW_TIMES_LIST == [(i,j) for i in wt for j in wt[i]])

