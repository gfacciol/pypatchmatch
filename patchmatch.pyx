import numpy as np
cimport numpy as np
from libcpp cimport bool


np.import_array()


#cdef extern from "patchmatch_simple.c":
cdef extern from "patchmatch_simple.cc":
   void patchmatch(float *u1_, int u1ncol, int u1nrow, int u1nch,
         float *u2_, int u2ncol, int u2nrow, int u2nch, 
         int w, char *method, int minoff,  int maxoff, 
         float *nnf, float *out_cost_ , 
         int iterations, int randtrials, bool color_planesi, bool use_horizontal_off)

def pm(np.ndarray[np.float32_t, ndim=3, mode='c'] u1,
       np.ndarray[np.float32_t, ndim=3, mode='c'] u2,
       np.ndarray[np.float32_t, ndim=3, mode='c'] nnf,
       np.ndarray[np.float32_t, ndim=2, mode='c'] cost,
       patchsz,
       minoff, 
       maxoff,
       n_iter=5,
       n_rand=5,
       method=u'SAD',   # SAD, SSD, ZSSD, ZSAD, NCC
       use_horizontal_off=False): 
   '''
   returns the NNF field and the COST
   '''

   cdef int w1 = u1.shape[1]
   cdef int h1 = u1.shape[0]
   cdef int c1 = u1.shape[2]
   cdef int w2 = u2.shape[1]
   cdef int h2 = u2.shape[0]
   cdef int c2 = u2.shape[2]

   methodstring = bytes(method, encoding='utf8')

   patchmatch(<float*>u1.data,w1,h1,c1, <float*>u2.data,w2,h2,c2, 
         patchsz, methodstring, minoff, maxoff,
         <float*>nnf.data, <float*>cost.data, 
         n_iter, n_rand, False, use_horizontal_off)

   return nnf, cost

