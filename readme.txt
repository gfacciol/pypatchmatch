SIMPLE PATCHMATCH CYTHON WRAPPER

build module
> python setup.py build

copy module in the right place where test-patchmatch.py is
> cp build/lib.macosx-10.11-x86_64-2.7/patchmatch.so patchmatch.so

run test program
> python test-patchmatch.py

patchmatch interface inside patchmatch.pyx
   def pm(np.ndarray[np.float32_t, ndim=3, mode='c'] u1,
          np.ndarray[np.float32_t, ndim=3, mode='c'] u2,
          np.ndarray[np.float32_t, ndim=3, mode='c'] nnf,
          np.ndarray[np.float32_t, ndim=2, mode='c'] cost,
          patchsz,         # patch size
          minoff,          # minimum offset
          maxoff,          # maximum offset
          n_iter=5,        # iterations
          n_rand=5,        # random trials per iteration
          method='SAD'):   # patch distance: SAD, SSD, ZSSD, ZSAD
      '''
      returns the NNF field and the COST
      '''

usage example
        import patchmatch as pm

        # read images img0,img1 as np.arrays
        ...

        # parameters
        patchSize = 7
        maxoff, minoff = 100, 30
        
        # fix types
        img0 = img0.astype(numpy.float32)
        img1 = img1.astype(numpy.float32)

        # create output arrays
        sz   = img0.shape;
        nnf  = numpy.ndarray((sz[0],sz[1],2)).astype(numpy.float32);
        cost = numpy.ndarray((sz[0],sz[1]  )).astype(numpy.float32);

        # call patchmatch
        pm.pm(img0, img1, nnf, cost, patchSize, numpy.int32(minoff), numpy.int32(maxoff))

        # here are the offsets
        dx = nnf[:,:,0] 
        dy = nnf[:,:,1] 


