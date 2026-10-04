"""Optional installed SuiteSparse AMD reference; never bundled into a port."""
import ctypes,ctypes.util

class SuiteSparseAMD:
    def __init__(self,path=None):
        path=path or ctypes.util.find_library('amd')
        if not path:raise RuntimeError('SuiteSparse AMD library not found; pass --amd-library')
        # Julia's macOS distribution uses an @rpath dependency; preload its sibling.
        from pathlib import Path
        sibling=Path(path).parent/'libsuitesparseconfig.dylib'
        if Path(path).is_absolute() and sibling.is_file():self.config=ctypes.CDLL(str(sibling),mode=ctypes.RTLD_GLOBAL)
        self.lib=ctypes.CDLL(str(path))
        integer=ctypes.c_int32;pointer=ctypes.POINTER(integer)
        self.order_fn=self.lib.amd_order
        self.order_fn.argtypes=[integer,pointer,pointer,pointer,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double)]
        self.order_fn.restype=ctypes.c_int
        self.lib.amd_version.argtypes=[pointer];self.lib.amd_version.restype=None
        v=(integer*3)();self.lib.amd_version(v);self.version='.'.join(map(str,v))
    def order(self,a):
        n=a['rows']
        if n!=a['cols'] or n>=2**31 or len(a['indices'])>=2**31:raise ValueError('reference requires square int32 CSR')
        integer=ctypes.c_int32
        rp=(integer*(n+1))(*a['offsets']);ci=(integer*max(1,len(a['indices'])))(*a['indices']);p=(integer*max(1,n))()
        # CSR(A) is CSC(Aᵀ); AMD symmetrizes either orientation identically.
        status=self.order_fn(n,rp,ci,p,None,None)
        if status not in (0,1):raise RuntimeError('SuiteSparse AMD failed with status '+str(status))
        out=list(p)[:n]
        if sorted(out)!=list(range(n)):raise RuntimeError('reference returned an invalid permutation')
        return out
