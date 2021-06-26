A simple implementation of the PatchMatch Algorithm [Barnes, 2009],
with its Cython and Matlab wrapper by Gabriele Facciolo (gfacciol@gmail.com)

PatchMatch is a randomized matching algorithm that allows to efficiently 
compute the correspondence (offset) map between two images. 

Build standalone program
========================
> mkdir build; cd build
> cmake ../; make

quick test
> ./patchmatch -t SSD -i 10 -w 7 ../{b,a}.png -R 26 offset.tif cost.tif backproj.png

Usage: ./patchmatch [-i iter(5)] [-w patch_size(7)] [-d init_off]
        [-t patchdist={SSD(default)|SAD|ZSSD|ZSAD|NCC}]
        [-r min_off(0)] [-R max_off(10)] u v offset [cost [backflow]]

The offset.tif output image will have two channels with float values,
for the column and row offset, respectively, from the left to the
right image.

This algorithm does not compute a subpixel offset or perform a
left-right consistency check.

Options:

-i iter
   Number of iterations. A larger offset may need more iterations to
   converge.

-w patch_size
    Use a patch (block) of this size around each pixel in the left
    image to match to the right image.

-t patchdist
    The distance function to use to measure the similarity between
    patches in the two images.

-R max_off
    Search for the offset between -max_off and max_off. The algorithm
    should not be given a value for this much larger than what is
    expected, as that may result in incorrect results in some places.

-r min_off
    Do not search for offsets between -min_off and min_off.
   
-h 
    Search only horizontal offsets, by default offsets are 2D.   


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


build matlab mex
================
> compileMex

run test
> a = double(imread('a.png'));
> b = double(imread('b.png'));
> [n,c] = patchmatchMex(a,b);

