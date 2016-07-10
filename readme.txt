A simple implementation of the PatchMatch Algorithm [Barnes, 2009],
with its Cython wrapper by Gabriele Facciolo (gfacciol@gmail.com)

PatchMatch is a randomized matching algorithm that allows to efficiently 
compute correspondece (offset) maps between two images. 

build standalone program
========================
> mkdir build; cd build
> cmake ../; make

quick test
> ./patchmatch -t SSD -w 5 ../{a,b}.png -R 300 offs.tif cost.tif backproj.png

usage: ./patchmatch [-i iter(5)] [-w patchsz(7)] [-d init_off]
        [-t patchdist={SSD(default)|SAD|ZSSD|ZSAD|NCC}]
        [-r min_off(0)] [-R max_off(10)] u v offsets [cost [backflow]]


build cython module
===================
> python setup.py build

copy module to the place where test-patchmatch.py is
> ln -s build/<lib.XXXXXXXXX-2.7>/patchmatch.so patchmatch.so

run test program
> python test-patchmatch.py


call patchmatch from python
===========================
the patchmatch interface inside patchmatch.pyx is
   def pm(np.ndarray[np.float32_t, ndim=3, mode='c'] u1,
          np.ndarray[np.float32_t, ndim=3, mode='c'] u2,
          np.ndarray[np.float32_t, ndim=3, mode='c'] nnf,
          np.ndarray[np.float32_t, ndim=2, mode='c'] cost,
          patchsz,         # patch size
          minoff,          # minimum offset
          maxoff,          # maximum offset
          n_iter=5,        # iterations
          n_rand=5,        # random trials per iteration
          method='SAD'):   # patch distance: SAD, SSD, ZSSD, ZSAD, NCC
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
        
        # fix input types
        img0 = img0.astype(numpy.float32)
        img1 = img1.astype(numpy.float32)

        # create output arrays (with the correct type)
        sz   = img0.shape;
        nnf  = numpy.ndarray((sz[0],sz[1],2)).astype(numpy.float32);
        cost = numpy.ndarray((sz[0],sz[1]  )).astype(numpy.float32);

        # call patchmatch
        pm.pm(img0, img1, nnf, cost, patchSize, minoff, maxoff)

        # here are the offsets
        dx = nnf[:,:,0] 
        dy = nnf[:,:,1] 


