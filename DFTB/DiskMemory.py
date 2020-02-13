# NOT WORKING PROPERLY YET
"""
Large numpy-arrays are written to a temporary file to economize memory usage. This mapping occurs only
for functions which are decorated with the @Mem.memmap decorator. All calls to np.zeros or np.asarray,
that would normally create new arrays in memory, are replaced by wrappers so that instead a memory map of this array 
is stored in a temporary file on disk. After the function has finished all temporary files are deleted.

This code is a little bit tricky because it modifies the namespace of the wrapped function so that
the name 'np.zeros' points to a different function that creates a memmap array.
"""
import resource
import numpy as _np
import tempfile
import os
from os.path import join

def memory_usage():
    """
    finds the amount of memory in MB used by the python process
    """
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000.0
    print "Memory used: %s MB" % mem
    return mem

class DiskMemory:
    _tmpdir = None
    # This dictionary contains the temporary files that were created by each function, in which the numpy
    # arrays are mapped to the hard disk. The names of the functions serve as 
    tmp_files = {}
    funcs_stack = []
    np_old = {"zeros": _np.zeros, "asarray": _np.asarray}
    def setScratchDirectory(self, tmpdir):
        self._tmpdir = tmpdir
    def memmap(self, func):
        def wrapper(*args, **kwds):
            if self._tmpdir == None:
                return func(*args, **kwds)
            else:
                print "Numpy arrays inside function %s will be mapped to hard disk" % func.__name__
                memory_usage()
                self.tmp_files[func.__name__] = []
                # The np-module imported in the namespace of the function definition
                np = func.__globals__["np"]
                # override functions for creating arrays, instead of ndarrays memmap-arrays are created
                np.zeros = self.zeros_memmap
                np.asarray = self.asarray_memmap
                # The name of th function for which the mapping is performed.
                self.funcs_stack.append( func.__name__ )
                # executed function
                res = func(*args, **kwds)
                # restore original numpy functions
                np.zeros = self.np_old["zeros"]
                np.asarray = self.np_old["asarray"]
                # convert results back to ndarrays
                if type(res) == type(()):
                    # tuple of results
                    res_nd = []
                    for r in res:
                        if type(r) == np.ndarray:
                            r = np.asarray(r)
                        res_nd.append( r )
                else:
                    if type(res) == np.core.memmap:
                        res_nd = np.asarray(res)
                    else:
                        res_nd = res
                # finished with this functions remove it from stack
                self.funcs_stack.pop()
                # remove temporary files created in func
                self.remove_tmp_files(func.__name__)
                memory_usage()

                return res_nd
            
        return wrapper
    def remove_tmp_files(self, func_name):
        for fh in self.tmp_files[func_name]:
            fh.close()
        del self.tmp_files[func_name]
    # WRAPPERS for np.zeros and np.asarray
    def zeros_memmap(self, shape, dtype=float, order='C'):
        fh = tempfile.TemporaryFile(dir=self._tmpdir, prefix="dftbaby_memmap")
        #fh,fname = tempfile.mkstemp(prefix=self._tmpdir)
        # create an array on disk that behaves in the same way as a ndarray
        arr = _np.memmap(fh, dtype=dtype, order=order, mode="w+", shape=shape)

        # keep track of the temporary files that were created for this function
        self.tmp_files[self.funcs_stack[-1]].append(fh)
        return arr
    def asarray_memmap(self, a, dtype=None, order=None):
        # if a is already a memmap object, nothing should be done
        if type(a) == _np.core.memmap:
            return a
        arr_nd = self.np_old["asarray"](a, dtype=dtype, order=order)
        fh = tempfile.TemporaryFile(dir=self._tmpdir, prefix="dftbaby_memmap")
        # create an array on disk that behaves in the same way as a ndarray
        arr_map = _np.memmap(fh, dtype=arr_nd.dtype, order=order, mode="w+", shape=arr_nd.shape)
        arr_map[:] = arr_nd[:]
        
        # keep track of the temporary files that were created for this function
        self.tmp_files[self.funcs_stack[-1]].append(fh)
        return arr_map

                
GlobalDiskMemory = DiskMemory()

if __name__ == "__main__":
    import numpy as np
    import time
    
    Mem = GlobalDiskMemory
    
    @Mem.memmap
    def test_function(x,y):
        from tempfile import mkdtemp
        import os.path as path
        
        imax = 2
        arrays = []
        for i in range(0, imax):
            #arrays.append( np.memmap(fname+"%d.dat" % i, dtype=float, mode="w+", shape=(1000,1000,1000)) )
            Z = np.zeros((1000,1000,10000), dtype=float)
            Z = np.asarray( Z )
            arrays.append( Z )
            memory_usage()
        for i in range(0, imax):
            print np.sum(arrays[i][:])
            memory_usage()
        
            
    test_function(1,2)
    
