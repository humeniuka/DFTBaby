"""
Record how much time is spent in each part of the code
using function decorators.

In order to time a function you have to attach a timer decorator to it

from Timer import GlobalTimer as T

@T.timer
def some_function(x):
    return x*x

print T

WARNING: for parallel excution only sequential functions are timed
"""
import time

class TimedFunction:
    def __init__(self, function_name):
        self.function_name = function_name
        self.nr_calls = 0
        self.total_time = 0.0
        self.relative_time = -1.0
    def add_timing(self, dt):
        self.nr_calls += 1
        self.total_time += dt
    def __str__(self):
        txt  = "Function %s\n" % self.function_name
        txt += "  number of calls        : %d\n" % self.nr_calls
        txt += "  total time spent inside: %4.7f seconds\n" % self.total_time
        if self.relative_time > 0.0:
            txt += "  fraction of total time : %4.7f seconds\n" % self.relative_time
        return txt

class Timer:
    def __init__(self):
        self.timed_functions = {}
    def timer(self, func):
        def wrapper(*args, **kwds):
            t1 = time.time()
            res = func(*args, **kwds)
            t2 = time.time()
            try:
                T = self.timed_functions[func.__name__]
            except KeyError: 
                T = TimedFunction(func.__name__)
                self.timed_functions[func.__name__] = T
            dt = t2-t1
            T.add_timing(dt)
            return res
        return wrapper
    def relative_times(self):
        total_time = 0.0
        for func_name,T in self.timed_functions.iteritems():
            total_time += T.total_time
        # set relative times
        for func_name,T in self.timed_functions.iteritems():
            T.relative_time = T.total_time/total_time
    def __str__(self):
#        self.relative_times()
        txt =  "Code Timing\n"
        txt += "===========\n"
        for func_name,T in self.timed_functions.iteritems():
            txt += str(T)
        return txt

GlobalTimer = Timer()

if __name__ == "__main__":
    T = GlobalTimer

    @T.timer
    def test(x):
        time.sleep(0.1)
        return x

    @T.timer
    def test2(x):
        return x

    class SomeClass:
        def __init__(self):
            pass
        @T.timer
        def bla(self, x):
            return x

    test(10)
    test(10)
    test2(10)

    C = SomeClass()
    C.bla(10)

    print T
