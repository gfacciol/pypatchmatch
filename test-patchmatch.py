'''
This code runs the patch match algorithm on two images
.....................................................
'''
import numpy
import pylab
import matplotlib;
import skimage
import skimage.io

import patchmatch as pm


files=["a.png","b.png"];

class PatchMatch:
    '''
    PatchMatch class
    '''

    def loadImages(self,files):
        '''
        load both images into self.img[0,1]
        '''
        self.img = [skimage.img_as_float(skimage.io.imread(files[i])).astype(numpy.float32) for i in (0,1)]

        
    def show(self,nffs=True):
        '''
        shows times and images
        '''
        if nffs:
            for i in range(0,3):
                pylab.imshow(self.nff[:,:,i])
                pylab.show();
        f=pylab.figure()
        f.add_subplot(1,2,1);    
        pylab.imshow(self.img[0],cmap=matplotlib.cm.Greys_r);
        f.add_subplot(1,2,2);
        pylab.imshow(self.img[1],cmap=matplotlib.cm.Greys_r);
        pylab.title("Patch Match")
        pylab.show();
    

    def match(self,files,patchSize=(5,5),iterations=20,minoff=0,maxoff=50):
        '''
        run the patchMatch algorithm on the images, returning nff array
        '''
        self.loadImages(files);

        self.size=self.img[0].shape;
        self.nff =numpy.ndarray((self.size[0],self.size[1],2)).astype(numpy.float32);
        self.cost=numpy.ndarray((self.size[0],self.size[1]  )).astype(numpy.float32);

        import time
        t = time.time()
        pm.pm(self.img[0], self.img[1], self.nff, self.cost, patchSize[0], minoff, maxoff)
        print time.time() - t

        # join nnf and costs 
        self.nff=numpy.dstack([self.nff, self.cost[:,:,numpy.newaxis]])

        if 1: # if DEMO
           self.show();
        return self.nff;


if __name__ == "__main__":
    patchmatch = PatchMatch()
    print("Please wait a few seconds...")
    patchmatch.match(files, iterations=20)
    print("Done.")
