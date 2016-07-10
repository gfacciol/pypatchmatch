/* Copyright 2014, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr> */

#include <assert.h>
#include <math.h>
#include "fimage.h"
#include "string.h"

// GLOBAL: RANDOM POSITIONS TESTED IN EACH ITERATION
int RANDOMTRIALS=5;



/* MACROS FOR SYMMETRIZED ACCESS */
/* dct type II symmetry at the boundary */
/* _psym  position of the symmetrized pixel */
/* _vsym  value of the symmetrized pixel */
//#define _psym(u,x,y,c)    ( _pos(u, ( (x)<0 ? -(x)-1: ( (x) >= (u)->ncol ? -(x)+2*(u)->ncol-1 : (x)  ) )  ,  ( (y)<0 ? -(y)-1: ( (y) >= (u)->nrow ? -(y)+2*(u)->nrow-1 : (y)  ) ) ,  (c) < (u)->nch ? c: -1  )  )
// SECURE against warp around
#define _pos(u,i,j,c) (  (i) + (j)*(u)->ncol + (u)->ncol*(u)->nrow*(c) )
#define _psym(u,x,y,c)    ( _pos(u, ( (x)<0 ?  \
( -(x)-1 >= (u)->ncol  ? (u)->ncol-1: -(x)-1 )  :  \
( (x) >= (u)->ncol ? (-(x)+2*(u)->ncol-1 < 0 ? 0 : -(x)+2*(u)->ncol-1)   : (x)  ) )  \
, ( (y)<0 ?  \
( -(y)-1 >= (u)->nrow  ? (u)->nrow-1: -(y)-1 )  :  \
( (y) >= (u)->nrow ? (-(y)+2*(u)->nrow-1 < 0 ? 0 : -(y)+2*(u)->nrow-1)   : (y)  ) ) \
,  (c) < (u)->nch ? c: -1  )  )
#define _vsym(u,x,y,c)    (  (u)->gray[ _psym(u,x,y,c) ]  )



float interp_nearest(const Fimage u, const float x, const float y, const int ch)
{
    int ix = round(x);
    int iy = round(y);
    return _vsym(u, ix, iy, ch);
}


float interp(const Fimage u, const float x, const float y, const int ch) {
    return interp_nearest(u,x,y,ch);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




/*
 * Computes the weighted L2 distance - mean:  (  sum k*[u1 - mu_1 - u2 + mu_2]^2  ) 
 * between the patch of the image 1 starting at pixel (ipx,ipy) 
 * and the patch of the image 2 starting at (iqx,iqy)
 * Where mu_1 and mu_2 are the averages of u1 and u2 respectively.
 * Assumes that sum (kernel) = 1
 * */
float distance_patch_weighted_l2mean(float *u1, float *u2, float *kernel, int k_w, int k_h, int nch)
{
   int N=k_w*k_h;
   // compute mean
   float mean[nch];

   for (int cc=0; cc < nch; cc++) 
   {
      float *ptr1 = u1 + cc*N;
      float *ptr2 = u2 + cc*N;
      float *ptrk = kernel;

      mean[cc] = 0;
      for (int yy = 0; yy < N; yy++)
      {
         mean[cc] += (*ptr1 - *ptr2)* *ptrk;
         ptr1++;
         ptr2++;
         ptrk++;
      }
   }

   // compute SSD-mean
   float fDist = 0.0f;
   for (int cc=0; cc < nch; cc++) 
   {
      float *ptr1 = u1 + cc*N;
      float *ptr2 = u2 + cc*N;
      float *ptrk = kernel;
      for (int yy = 0; yy < N; yy++)
      {
         float dif = *ptr1 - *ptr2 - mean[cc];
         fDist += *ptrk * dif * dif / nch;
         ptr1++;
         ptr2++;
         ptrk++;
      }
   }

   return fDist;
}


/*
 * Computes the weighted L1 distance - mean:  (  sum k*|u1 - mu_1 - u2 + mu_2|  ) 
 * between the patch of the image 1 starting at pixel (ipx,ipy) 
 * and the patch of the image 2 starting at (iqx,iqy)
 * Where mu_1 and mu_2 are the averages of u1 and u2 respectively.
 * Assumes that sum (kernel) = 1
 * */
float distance_patch_weighted_l1mean(float *u1, float *u2, float *kernel, int k_w, int k_h, int nch)
{
   int N=k_w*k_h;
   // compute mean
   float mean[nch];

   for (int cc=0; cc < nch; cc++) 
   {
      float *ptr1 = u1 + cc*N;
      float *ptr2 = u2 + cc*N;
      float *ptrk = kernel;

      mean[cc] = 0;
      for (int yy = 0; yy < N; yy++)
      {
         mean[cc] += (*ptr1 - *ptr2)* *ptrk;
         ptr1++;
         ptr2++;
         ptrk++;
      }
   }

   // compute SAD-mean
   float fDist = 0.0f;
   for (int cc=0; cc < nch; cc++) 
   {
      float *ptr1 = u1 + cc*N;
      float *ptr2 = u2 + cc*N;
      float *ptrk = kernel;
      for (int yy = 0; yy < N; yy++)
      {
         float dif = *ptr1 - *ptr2 - mean[cc];
         fDist += *ptrk * (dif>=0?dif:-dif) / nch;
         ptr1++;
         ptr2++;
         ptrk++;
      }
   }

   return fDist;
}


/*
 * Computes the weighted L2 distance between the patch of the image 1 starting at pixel (ipx,ipy)
 * and the patch of the image 2 starting at (iqx,iqy)
 * */
float distance_patch_weighted_l2(float *u1, float *u2, float *kernel, int k_w, int k_h, int nch, float curr_min_dist)
{
    int N=k_w*k_h;
    
    // compute SSD
    float fDist = 0.0f;
    for (int cc=0; cc < nch; cc++)
    {
        float *ptr1 = u1 + cc*N;
        float *ptr2 = u2 + cc*N;
        float *ptrk = kernel;
        for (int yy = 0; yy < N && fDist<=curr_min_dist ; yy++)
        {
            float dif = *ptr1 - *ptr2;
            fDist += *ptrk * dif * dif / nch;
            ptr1++;
            ptr2++;
            ptrk++;
        }
    }
    
    return fDist;
}

/*
 * Computes the weighted L1 distance between the patch of the image 1 starting at pixel (ipx,ipy)
 * and the patch of the image 2 starting at (iqx,iqy)
 * */
float distance_patch_weighted_l1(float *u1, float *u2, float *kernel, int k_w, int k_h, int nch, float curr_min_dist)
{
    int N=k_w*k_h;
    
    float fDist = 0.0f;
    
    for (int cc=0; cc < nch; cc++)
    {
        float *ptr1 = u1 +N*cc;
        float *ptr2 = u2 +N*cc;
        float *ptrk = kernel;
        
        for (int yy = 0; yy < N && fDist<=curr_min_dist ; yy++)
        {
            float dif = *ptr1 - *ptr2;
            fDist += *ptrk * (dif>=0?dif:-dif) / nch;
            ptr1 ++;
            ptr2 ++;
            ptrk ++;
        }
    }
    
    return fDist;
}



float distance_patch(float *u1, float *u2, float *kernel, int k_w, int k_h, int nch, char * method, float curr_min_dist)
{
   float S;
   if (strcmp (method,"SSD")==0)
      S = distance_patch_weighted_l2(u1, u2, kernel, k_w, k_h, nch, curr_min_dist);
   else if (strcmp (method,"SAD")==0)
      S = distance_patch_weighted_l1(u1, u2, kernel, k_w, k_h, nch, curr_min_dist);
   else if (strcmp (method,"ZSSD")==0)
      S = distance_patch_weighted_l2mean(u1, u2, kernel, k_w, k_h, nch);
   else if (strcmp (method,"ZSAD")==0)
      S = distance_patch_weighted_l1mean(u1, u2, kernel, k_w, k_h, nch);
   else
   {
      fprintf(stderr, "UNKNOWN METHOD %s\n", method); exit(-1);
   }
   return S;
}




// ker is just needed for the size of the patch
// extract a patch where the channels are planes
void extract_patch_integer(Fimage u, int x, int y, Fimage ker, float *p)
{
    int knc  = ker->ncol;
    int knr  = ker->nrow;
    int halfknc = knc/2;
    int halfknr = knr/2;
    int nch = u->nch;
    int nc  = u->ncol;
    int nr  = u->nrow;
    int a,b,c;
    
    for (c=0;c<nch;c++)
        for (b=0;b<knr;b++)
            for (a=0;a<knc;a++)
            {
                *p++ = _vsym(u,x+a-halfknc,y+b-halfknr,c);
            }
    
}

// ker is just needed for the size of the patch
// extract a patch where the channels are planes
void extract_patch_integer_noboundary(Fimage u, const int x, const int y, Fimage ker, float *p)
{
    const int knc  = ker->ncol;
    const int knr  = ker->nrow;
    const int halfknc = knc/2;
    const int halfknr = knr/2;
    const int nch = u->nch;
    const int nc  = u->ncol;
    const int nr  = u->nrow;
    int a,b,c;
    
    for (c=0;c<nch;c++) {
        float *pu = u->gray + nc*nr*c;
        
        for (b=0;b<knr;b++) {
            float *ppu = pu + nc*(y+b-halfknr);
            
            ppu += x - halfknc;
            for (a=0;a<knc;a++) {
                *p++ = *ppu++;
            }
        }
    }
    
}


// ker is just needed for the size of the patch
// extract a patch where the channels are planes
void extract_patch_interp(Fimage u, float x, float y, Fimage ker, float *p)
{
    int knc  = ker->ncol;
    int knr  = ker->nrow;
    int halfknc = knc/2;
    int halfknr = knr/2;
    int nch = u->nch;
    int nc  = u->ncol;
    int nr  = u->nrow;
    int a,b,c;
    
    for (c=0;c<nch;c++)
        for (b=0;b<knr;b++)
            for (a=0;a<knc;a++)
            {
                *p++ = interp(u, x+a-halfknc,y+b-halfknr,c);
            }
    
}


// ker is just needed for the size of the patch
// extract a patch where the channels are planes
void extract_patch_secure(Fimage u, const float x, const float y, Fimage ker, float *p)
{
    const int knc  = ker->ncol;
    const int knr  = ker->nrow;
    const int halfknc = knc/2;
    const int halfknr = knr/2;
    const int nch = u->nch;
    const int nc  = u->ncol;
    const int nr  = u->nrow;
    
    if ((x-(int)x == 0) && (y-(int)y == 0)) {
        if(x-halfknc<0 || x+halfknc>=nc || y-halfknr<0 || y+halfknr>=nr)
            return extract_patch_integer(u, (int) x, (int) y, ker, p);
        else
            return extract_patch_integer_noboundary(u, (int) x, (int) y, ker, p);
    }
    
    return extract_patch_interp(u, x, y, ker, p);
    
}





//                                                                              vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv  these are temporary
void random_search(Fimage u1, Fimage u2, Fimage kernel, int w, int nch, char* method,  float *dx, float* dy, float *out_cost, int minoff, int maxoff)
{
    int nc = u1->ncol;
    int nr = u1->nrow;
    int pdim = w*w*nch;
    int maxtrials= RANDOMTRIALS;
#pragma omp parallel for
    for (int y=0;y<nr;y++) {
        for (int x=0;x<nc;x++) {
            int pos =x+y*nc;
            float curr_dx = dx[pos];
            float curr_dy = dy[pos];
            float curr_cost = out_cost[pos];
            float p1[pdim], p2[pdim];
            extract_patch_secure(u1, x, y, kernel, p1);
            
            for (int trials=0;trials<maxtrials; trials++) {
                float offx = (int) ( ( ( (double) rand() ) / ((double) RAND_MAX + 1.0) ) *2* maxoff) - maxoff;
                float offy = (int) ( ( ( (double) rand() ) / ((double) RAND_MAX + 1.0) ) *2* maxoff) - maxoff;
                float tmpx=x+offx;
                float tmpy=y+offy;
                
                // do not accept offsets smaller than
                if( offx*offx + offy*offy < minoff) continue;
                
                // only if the random patch is inside the image
                if(tmpx>=0 && tmpx<u2->ncol && tmpy>=0 && tmpy<u2->nrow) {
                    extract_patch_secure(u2, tmpx, tmpy, kernel, p2);
                    float new_cost = distance_patch(p1, p2, kernel->gray, w,w,nch, method, curr_cost);
                    
                    if( new_cost < curr_cost ) {
                        curr_cost=new_cost;
                        curr_dx = offx;
                        curr_dy = offy;
                    }
                } //else {printf("x");}
                
                
            }
            out_cost[pos] = curr_cost;
            dx[pos] = curr_dx;
            dy[pos] = curr_dy;
            
        }
    }
}

//                                                                              vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv  these are temporary
void backward_propagation(Fimage u1, Fimage u2, Fimage kernel, int w, int nch, char* method,  float *dx, float* dy, float *out_cost, int minoff, int maxoff)
{
    int nc = u1->ncol;
    int nr = u1->nrow;
    int pdim = w*w*nch;
#pragma omp parallel for
    for (int y=nr-1;y>=0;y--)
        for (int x=nc-1;x>=0;x--) {
            int pos =x+y*nc;
            float curr_dx = dx[pos];
            float curr_dy = dy[pos];
            float curr_cost = out_cost[pos];
            
            float p1[pdim], p2[pdim];
            extract_patch_secure(u1, x, y, kernel, p1);
            
            // scan the neighbors (forward set)
            int neighsx[] = { 0, 1, 1,-1};
            int neighsy[] = { 1, 0, 1, 1};
            for (int trials=0;trials<4; trials++) {
                // position of the neighbor
                int neighx = x+ neighsx[trials];
                int neighy = y+ neighsy[trials];
                
                if( !(neighx>=0 && neighx<u1->ncol && neighy>=0 && neighy<u1->nrow)) continue;
                float offx = dx[neighx + neighy*nc];
                float offy = dy[neighx + neighy*nc];
                
                float tmpx=x+offx;
                float tmpy=y+offy;
                
                // do not accept offsets smaller than
                if( offx*offx + offy*offy < minoff) continue;
                
                // only if the tested patch is inside the image
                if(tmpx>=0 && tmpx<u2->ncol && tmpy>=0 && tmpy<u2->nrow) {
                    extract_patch_secure(u2, tmpx, tmpy, kernel, p2);
                    float new_cost = distance_patch(p1, p2, kernel->gray, w,w,nch, method, curr_cost);
                    
                    if( new_cost < curr_cost ) {
                        curr_cost=new_cost;
                        curr_dx = offx;
                        curr_dy = offy;
                    }
                }
                
            }
            out_cost[pos] = curr_cost;
            dx[pos] = curr_dx;
            dy[pos] = curr_dy;
            
        }
}



//                                                                              vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv  these are temporary
void forward_propagation(Fimage u1, Fimage u2, Fimage kernel, int w, int nch, char* method,  float *dx, float* dy, float *out_cost, int minoff, int maxoff)
{
    int nc = u1->ncol;
    int nr = u1->nrow;
    int pdim = w*w*nch;
#pragma omp parallel for
    for (int y=0;y<nr;y++)
        for (int x=0;x<nc;x++) {
            int pos =x+y*nc;
            float curr_dx = dx[pos];
            float curr_dy = dy[pos];
            float curr_cost = out_cost[pos];
            
            float p1[pdim], p2[pdim];
            extract_patch_secure(u1, x, y, kernel, p1);
            
            // scan the neighbors (forward set)
            int neighsx[] = { 0,-1, 1,-1};
            int neighsy[] = {-1, 0,-1,-1};
            for (int trials=0;trials<4; trials++) {
                // position of the neighbor
                int neighx = x+ neighsx[trials];
                int neighy = y+ neighsy[trials];
                
                if( !(neighx>=0 && neighx<u1->ncol && neighy>=0 && neighy<u1->nrow)) continue;
                float offx = dx[neighx + neighy*nc];
                float offy = dy[neighx + neighy*nc];
                
                float tmpx=x+offx;
                float tmpy=y+offy;
                
                // do not accept offsets smaller than
                if( offx*offx + offy*offy < minoff) continue;
                
                // only if the tested patch is inside the image
                if(tmpx>=0 && tmpx<u2->ncol && tmpy>=0 && tmpy<u2->nrow) {
                    extract_patch_secure(u2, tmpx, tmpy, kernel, p2);
                    
                    float new_cost = distance_patch(p1, p2, kernel->gray, w,w,nch, method, curr_cost);
                    
                    if( new_cost < curr_cost ) {
                        curr_cost=new_cost;
                        curr_dx = offx;
                        curr_dy = offy;
                    }
                }
                
            }
            out_cost[pos] = curr_cost;
            dx[pos] = curr_dx;
            dy[pos] = curr_dy;
            
        }
}


void patchmatch(float *u1_, int u1ncol, int u1nrow, int u1nch, 
                float *u2_, int u2ncol, int u2nrow, int u2nch, 
                int w, char *method, int minoff,  int maxoff, 
                float *nnf_, float *out_cost_, 
                int iterations, int randomtrials)
{
    RANDOMTRIALS=randomtrials;
    int nch = u1nch;
    int npix1 = u1ncol*u1nrow;
    int npix2 = u2ncol*u2nrow;
    int nc = u1ncol, nr = u1nrow;
    Fimage u1       = new_fimage3(u1_, u1ncol ,u1nrow ,u1nch);
    Fimage u2       = new_fimage3(u2_, u2ncol ,u2nrow ,u2nch);
    Fimage dx       = new_fimage3(nnf_      , nc , nr, 1);
    Fimage dy       = new_fimage3(nnf_+nc*nr, nc , nr, 1);
    Fimage nnf      = new_fimage3(nnf_      , nc , nr, 2); //to be used later
    Fimage out_cost = new_fimage3(out_cost_ , nc , nr, 1);

    // fix interfacing convenction: from vector pixels to color planes
    fimage_vec_to_planar(u1);
    fimage_vec_to_planar(u2);
    
    // generate a flat window
    Fimage kernel = new_fimage3(NULL, w ,w ,1);
    for(int i=0;i<w*w  ;i++)
        kernel->gray[i] = 1.0/(w*w);
    
    printf("enter %s\n", method);
    
    for(int i=0;i<npix1;i++)
        out_cost->gray[i] = INFINITY;
    
    for (int i=0;i<iterations;i++)
    {
        printf("iteration %d\n",i);
        // random search
        random_search(u1, u2, kernel, w, nch, method,  dx->gray,  dy->gray, out_cost->gray, minoff,maxoff);
        
        // forward propagation
        forward_propagation(u1, u2, kernel, w, nch, method,  dx->gray,  dy->gray, out_cost->gray,minoff, maxoff);
        
        // backward propagation
        backward_propagation(u1, u2, kernel, w, nch, method,  dx->gray,  dy->gray, out_cost->gray,minoff, maxoff);
    }

    // cleanup the interfacing mess I just did (planar to vector)
    fimage_planar_to_vec(u1);
    fimage_planar_to_vec(u2);
    fimage_planar_to_vec(nnf);
    u1->gray       = NULL; del_fimage(u1);
    u2->gray       = NULL; del_fimage(u2);
    dx->gray       = NULL; del_fimage(dx);
    dy->gray       = NULL; del_fimage(dy);
    nnf->gray      = NULL; del_fimage(nnf);
    out_cost->gray = NULL; del_fimage(out_cost);
    del_fimage(kernel);
}





