#include <assert.h>
#include <math.h>
#include <mex.h>
#include <stdio.h>


void help() {
   mexErrMsgTxt("patchmatch Matlab wrapper\n\n"
         "[nnf,cost] = patchmatch(u1, u2, [w(5)], [minoff(0)], maxoff(imsize)], [iterations(5)], [randomtrials(5)], [method('SSD')], [use_horizontal_off(false)], [nnfin])\n"
         "   \n"
         "   u1,u2        input DOUBLE images not necessarily the same size (but same #channels)\n"
         "   w            window size (default:5)\n"
         "   minoff       minimum offset (default:0)\n"
         "   maxoff       maximum offset (default:hypot(nx,ny))\n"
         "   iterations   default:5\n"
         "   randomtrials default:5\n"
         "   method       one of SSD, SAD, ZSSD, ZSAD, NCC (default:'SSD')\n"
         "   use_horizontal_off   compute horizontal 1d nnf or 2d one (default:0, false)\n"
         "   nnfin        initial nnf (default:none)\n"
         "   \n"
         "   nnf          offset map\n"
         "   cost         cost map\n");
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* TODO: check inputs and nlhs */
    if(nrhs<2) help();

    for (int i=1;i<10;i++)
       if(i!=8) // skip the char string
         if (nrhs>=i && (!mxIsDouble(prhs[i-1])))
         {
            mexPrintf("input %d must be double\n\n", i);
            help(); 
         }

    /* read input images */
    const mwSize* dims = mxGetDimensions(prhs[0]);
    int dimy=1, dimx=1, dimc=1; 
    if(mxGetNumberOfDimensions(prhs[0])==2) {
      dimy = (int)dims[0]; dimx = (int)dims[1]; 
    } else if(mxGetNumberOfDimensions(prhs[0])==3) {
      dimy = (int)dims[0]; dimx = (int)dims[1]; dimc = (int)dims[2];
    } else
      mexErrMsgTxt("something is wrong with image 1\n\n");
    double* u1 = mxGetPr(prhs[0]);


    /* read input images */
    const mwSize* dims2 = mxGetDimensions(prhs[1]);
    int dim2y=1, dim2x=1, dim2c=1; 
    if(mxGetNumberOfDimensions(prhs[1])==2) {
      dim2y = (int)dims2[0]; dim2x = (int)dims2[1]; 
    } else if(mxGetNumberOfDimensions(prhs[1])==3) {
      dim2y = (int)dims2[0]; dim2x = (int)dims2[1]; dim2c = (int)dims2[2];
    } else
      mexErrMsgTxt("something is wrong with image 2\n\n");
    double* u2 = mxGetPr(prhs[1]);


    /* read input images */
    if(dimc!=dim2c)
      mexErrMsgTxt("images must have the same number of channels\n\n");


    // read the other params
    int w=5;
    int minoff=0;
    int maxoff=hypot(dimx,dimy);
    int iterations=5; // default 5
    int randomtrials=5; // default 5
    char* method=(char *)"SSD";  // default SSD
    bool use_horizontal_off=false; // default false
    {
    int i=3;
    if (nrhs>=i) w=mxGetPr(prhs[i-1])[0]; i++;
    if (nrhs>=i) minoff=mxGetPr(prhs[i-1])[0]; i++;
    if (nrhs>=i) maxoff=mxGetPr(prhs[i-1])[0]; i++;
    if (nrhs>=i) iterations=mxGetPr(prhs[i-1])[0]; i++;
    if (nrhs>=i) randomtrials=mxGetPr(prhs[i-1])[0]; i++;
    if (nrhs>=i) method = mxArrayToString(prhs[i-1]); i++;
    if (nrhs>=i) use_horizontal_off = mxGetPr(prhs[i-1])[0]; 
    }


    // create output matrices
    const mwSize nnfdims[3] = {dimy, dimx, 2};
    mxArray* nnfm = plhs[0] = mxCreateNumericArray(3, nnfdims, mxDOUBLE_CLASS, mxREAL);
    double *nnf=mxGetPr(nnfm);
    mxArray* costm= plhs[1] = mxCreateDoubleMatrix(dimy, dimx, mxREAL);
    double *cost=mxGetPr(costm);


    // read input nnf
    if(nrhs>=10){
      const mwSize* dims = mxGetDimensions(prhs[9]);
      double* tmp = mxGetPr(prhs[9]);
      if (dims[0]==nnfdims[0]  && dims[1]==nnfdims[1]  && dims[2]==nnfdims[2]) 
         for(int i=0;i<dimy*dimx*2;i++) { nnf[i]=tmp[i]; }
      else
         mexPrintf("input nnf has the wrong size\n");
    }


    float *u1_  =(float*)malloc(dimy*dimx*dimc*sizeof(float));
    float *u2_  =(float*)malloc(dim2y*dim2x*dim2c*sizeof(float));
    float *nnf_ =(float*)malloc(dimy*dimx*2*sizeof(float));
    float *cost_=(float*)malloc(dimy*dimx*sizeof(float));
    for(int i=0;i<dimy*dimx*dimc;i++) { u1_[i]=u1[i]; }
    for(int i=0;i<dim2y*dim2x*dim2c;i++) { u2_[i]=u2[i]; }
    
    void patchmatch(float *u1_, int nc  , int nr  , int nch, 
                float *u2_, int u2nc, int u2nr, int u2nch, 
                int w, char *method, int minoff,  int maxoff, 
                float *nnf_, float *out_cost_, 
                int iterations, int randomtrials, bool planes, bool use_horizontal_off);

    /* swap dimensions because transposed */
    patchmatch(u1_, dimy, dimx, dimc, 
          u2_, dim2y, dim2x, dim2c,
          w, method, minoff, maxoff, 
          nnf_, cost_, iterations, randomtrials, true, use_horizontal_off);

    // swap dx and dy in nnf
    for(int i=0;i<dimy*dimx;i++) { nnf[i]=nnf_[i+dimy*dimx]; nnf[i+dimy*dimx]=nnf_[i];}
    for(int i=0;i<dimy*dimx;i++) { cost[i]=cost_[i]; }
    free(u1_);
    free(u2_);
    free(nnf_);
    free(cost_);
}
