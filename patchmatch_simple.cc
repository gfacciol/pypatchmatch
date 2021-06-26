/* Copyright (C) 2016, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>*/
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include <cmath> // std::abs for float and double
#include <numeric>
#include <vector>

#include "point.h"
#ifdef _OPENMP 
#include <omp.h>
#endif //_OPENMP

//// a structure to wrap images 
#include "img.h"


inline int check_inside_image(const int x, const int y, const struct Img &u) 
{
	if(x>=0 && y>=0 && x<u.nx && y<u.ny) return 1;
	else return 0;
}

inline float valnan(const struct Img &u, const int px, const int py, const int ch=0)
{
	return check_inside_image(px, py, u) ? u(px, py, ch) : NAN;
}

inline float valneumann(const struct Img &u, const int x, const int y, const int ch=0)
{  
   int xx=x, yy=y;
   xx = x >=  0  ? xx : 0;
   xx = x < u.nx ? xx : u.nx - 1;
   yy = y >=  0  ? yy : 0;
   yy = y < u.ny ? yy : u.ny - 1;
	return u(xx,yy,ch);
}

// // not used here but generally useful 
// typedef std::vector<float> FloatVector;

int RANDOMTRIALS=5;

inline float distance_patch_SxD(const Img &p1, const Img &p2, const float &curr_cost, const int type)
{
   float dist = 0.;
   for(int i=0;i<p1.npix*p1.nch;i++){
      float d = p1[i]-p2[i];
      dist += (type == 1 ? std::abs(d) : d*d);
//      if(dist > curr_cost) return dist;
   }
   return dist;
}
inline float distance_patch_SSD(const Img &p1, const Img &p2, const float &curr_cost)
{
   return distance_patch_SxD(p1, p2, curr_cost, 2);
}
inline float distance_patch_SAD(const Img &p1, const Img &p2, const float &curr_cost)
{
   return distance_patch_SxD(p1, p2, curr_cost, 1);
}


inline float distance_patch_ZSxD(const Img &p1, const Img &p2, const float &curr_cost, const int type)
{
   float dist = 0.;
   for(int c=0;c<p1.nch;c++){
      float mu1=0, mu2=0;
      for(int i=0;i<p1.npix;i++){
         mu1 += p1[i+p1.npix*c];
         mu2 += p2[i+p1.npix*c];
      }
      mu1/=p1.npix;
      mu2/=p1.npix;
      for(int i=0;i<p1.npix;i++){
         float d = std::abs(p1[i]-mu1-(p2[i]-mu2));
         dist += (type == 1 ? std::abs(d) : d*d);
//         if(dist > curr_cost) return dist;
      }
   }
   return dist;
}
inline float distance_patch_ZSSD(const Img &p1, const Img &p2, const float &curr_cost)
{
   return distance_patch_ZSxD(p1, p2, curr_cost, 2);
}
inline float distance_patch_ZSAD(const Img &p1, const Img &p2, const float &curr_cost)
{
   return distance_patch_ZSxD(p1, p2, curr_cost, 1);
}

inline float distance_patch_NCC(const Img &p1, const Img &p2, const float &curr_cost)
{
   float dist = 0.;
   for(int c=0;c<p1.nch;c++){
      int choffset = p1.npix*c;
      float mu1=0, mu2=0;        // means
      for(int i=0;i<p1.npix;i++){
         mu1 += p1[i+choffset];
         mu2 += p2[i+choffset];
      }
      mu1/=p1.npix;
      mu2/=p1.npix;
      double sig1=0, sig2=0;     // variances
      double distl = 0.;
      for(int i=0;i<p1.npix;i++){
         float d1 = p1[i+choffset]-mu1;
         float d2 = p2[i+choffset]-mu2;
         float d = (d1-d2);
         sig1  += d1*d1; 
         sig2  += d2*d2; 
         distl += d*d;
      }
      dist += distl/(sqrt(sig1)*sqrt(sig2));
      if(dist > curr_cost) return dist;
   }
   return dist;
}

// ker is just needed for the size of the patch
// extract a patch where the channels are planes
inline void extract_patch_integer(const Img &u, int x, int y, int knc, int knr , Img &p)
{
    int halfknc = knc/2;
    int halfknr = knr/2;
    int nch = u.nch;
    int nc  = u.ncol;
    int nr  = u.nrow;
    int a,b,c,i=0;

    for (c=0;c<nch;c++)
        for (b=0;b<knr;b++)
            for (a=0;a<knc;a++)
            {
                 p[i++] = valneumann(u,x+a-halfknc,y+b-halfknr,c);
            }
}

// ker is just needed for the size of the patch
// extract a patch where the channels are planes
inline void extract_patch_integer_noboundary(const Img &u, const int x, const int y, int knc, int knr , Img &p)
{
    const int halfknc = knc/2;
    const int halfknr = knr/2;
    const int nch = u.nch;
    const int nc  = u.ncol;
    const int nr  = u.nrow;
    int a,b,c,i=0;
    
    for (c=0;c<nch;c++) {
        int pu = nc*nr*c;
        
        for (b=0;b<knr;b++) {
            int ppu = pu + nc*(y+b-halfknr);
            
            ppu += x - halfknc;
            for (a=0;a<knc;a++) {
                p[i++] = u[ppu++];
            }
        }
    }
}


// ker is just needed for the size of the patch
// extract a patch where the channels are planes
inline void extract_patch_secure(const Img &u, const float x, const float y, Img &p)
{
    const int knc = p.ncol, knr = p.nrow;
    const int halfknc = knc/2;
    const int halfknr = knr/2;
    const int nch = u.nch;
    const int nc  = u.ncol;
    const int nr  = u.nrow;
    
    if ((x-(int)x == 0) && (y-(int)y == 0)) {
        if(x-halfknc<0 || x+halfknc>=nc || y-halfknr<0 || y+halfknr>=nr)
            return extract_patch_integer(u, (int) x, (int) y, knc, knr, p);
        else
            return extract_patch_integer_noboundary(u, (int) x, (int) y, knc, knr, p);
    }
    
    // still integer
    return extract_patch_integer(u, (int) x, (int) y, knc, knr, p);
}


typedef float(*patch_distance_func)(const Img &, const Img &, const float&); //signature of all patch distance functions

template<patch_distance_func dist>
void random_search(Img &u1, Img &u2, int w, Img &off, Img &cost, int minoff, int maxoff, bool use_horizontal_off)
{
    int maxtrials = RANDOMTRIALS;
    int thmax = 1; // thread handling data
    #ifdef _OPENMP 
    thmax = omp_get_max_threads();
    #endif
    std::vector<Img> p1(thmax, Img(w,w,u1.nch)), 
                     p2(thmax, Img(w,w,u1.nch)); 
    // random seeds
    static std::vector<unsigned int> seeds(thmax);
    for(int i=0; i<thmax; i++) seeds[i]+=i;

#pragma omp parallel for shared(p1,p2) 
    for (int y=0;y<u1.ny;y++) {
        for (int x=0;x<u1.nx;x++) {
            int thid = 0;
            #ifdef _OPENMP // thread id
            thid = omp_get_thread_num();
            #endif
            Point curr(off(x,y,0), off(x,y,1));
            float curr_cost = cost(x,y);
            extract_patch_secure(u1, x, y, p1[thid]);
            
            for (int trial=0;trial<maxtrials; trial++) {
                Point off( ((int) ( ( ( (double) rand_r(&seeds[thid]) ) / ((double) RAND_MAX + 1.0) ) *2* maxoff) - maxoff),
                           ((int) ( ( ( (double) rand_r(&seeds[thid]) ) / ((double) RAND_MAX + 1.0) ) *2* maxoff) - maxoff) * (1-use_horizontal_off) );
                
                Point tmp = Point(x,y) + off;
                // skip points that fell outside the image or offsets smaller than minoff 
                if( hypot(off[0],off[1]) < minoff || 
                    ! check_inside_image(tmp[0],tmp[1], u2)  ) continue;
                
                extract_patch_secure(u2, tmp[0], tmp[1], p2[thid]);
                float new_cost = dist(p1[thid], p2[thid], curr_cost);
                
                if( new_cost < curr_cost ) {
                    curr_cost=new_cost;
                    curr = off;
                }
                
            }
            cost(x,y) = curr_cost;
            off(x,y,0) = curr[0];
            off(x,y,1) = curr[1];
        }
    }
}

template<patch_distance_func dist>
void propagation(Img &u1, Img &u2, int w, Img &off, Img &cost, int minoff, int maxoff, const int direction)
{
    int maxtrials = RANDOMTRIALS;
    int thmax = 1; // thread handling data
    #ifdef _OPENMP 
    thmax = omp_get_max_threads();
    #endif
    std::vector<Img> p1(thmax, Img(w,w,u1.nch)), 
                     p2(thmax, Img(w,w,u1.nch)); 

    // setup scan direction
    const int tbase= direction == 1? 4 : 0;
    const int fx   = direction == 1? 0 : u1.nx-1;
    const int fy   = direction == 1? 0 : u1.ny-1;
#pragma omp parallel for shared(p1,p2)
    for (int j=0;j<u1.ny;j++) {
        int y = fy + direction * j;
        for (int i=0;i<u1.nx;i++) {
            int x = fx + direction * i;
            int thid = 0;
            #ifdef _OPENMP // thread id
            thid = omp_get_thread_num();
            #endif
            Point curr(off(x,y,0), off(x,y,1));
            float curr_cost = cost(x,y);

            extract_patch_secure(u1, x, y, p1[thid]);
            
            // scan the neighbors (backward set)   (forward set)
            static const Point neighs[] = {Point(0,1),  Point(1,0),  Point(1,1),   Point(-1,1), 
                                           Point(0,-1), Point(-1,0), Point(-1,-1), Point(1,-1)};
            for (int trial=0;trial<4; trial++) {
                // position of the neighbor
                Point neigh = Point(x,y) + neighs[trial+tbase];
                
                if( !check_inside_image(neigh[0],neigh[1], u1)) continue;
                Point noff( off(neigh[0], neigh[1], 0), off(neigh[0], neigh[1], 1) );
                
                Point tmp = Point(x,y) + noff;
                // skip points that fell outside the image or offsets smaller than minoff 
                if( hypot(noff[0],noff[1]) < minoff || 
                    ! check_inside_image(tmp[0],tmp[1], u2)  ) continue;
                
                extract_patch_secure(u2, tmp[0], tmp[1], p2[thid]);

                float new_cost = dist(p1[thid], p2[thid], curr_cost);
                if( new_cost < curr_cost ) {
                    curr_cost=new_cost;
                    curr = noff;
                }
                
            }
            cost(x,y) = curr_cost;
            off(x,y,0) = curr[0];
            off(x,y,1) = curr[1];
        }
    }
}


template<patch_distance_func dist>
void patchmatch(Img &u1, Img &u2, int w, Img &off, Img &cost,int minoff, int maxoff,  int iterations, int randomtrials, bool use_horizontal_off)
{
    RANDOMTRIALS=randomtrials;

    cost.setvalues(INFINITY);
    
    //srand(0);
    for (int i=0;i<iterations;i++)
    {
        printf("iteration %d\n",i);
        // random search
        random_search<dist>(u1, u2, w, off, cost, minoff,maxoff, use_horizontal_off);
        // forward propagation
        propagation<dist>(u1, u2, w, off, cost, minoff,maxoff, 1); 
        // backward propagation
        propagation<dist>(u1, u2, w, off, cost, minoff,maxoff, -1); 
    }
}


// interface used by python and matlab
void patchmatch(float *u1_, int nc  , int nr  , int nch, 
                float *u2_, int u2nc, int u2nr, int u2nch, 
                int w, char *method, int minoff,  int maxoff, 
                float *nnf_, float *out_cost_, 
                int iterations, int randomtrials, bool channels_as_planes, 
                bool use_horizontal_off)
{
    // fix interfacing convenction: from vector pixels to color planes
    Img u1(u1_  ,      nc  ,   nr, nch  , channels_as_planes);
    Img u2(u2_  ,      u2nc, u2nr, u2nch, channels_as_planes);
    Img nnf(nnf_,      nc  ,   nr, 2    , channels_as_planes); //to be used later
    Img cost(nc, nr, 1);

    printf("enter %s\n", method);
   if (strcmp (method,"SSD")==0)
      patchmatch<distance_patch_SSD>(u1, u2, w, nnf, cost, minoff, maxoff, iterations, randomtrials, use_horizontal_off);
   else if (strcmp (method,"SAD")==0)
      patchmatch<distance_patch_SAD>(u1, u2, w, nnf, cost, minoff, maxoff, iterations, randomtrials, use_horizontal_off);
   else if (strcmp (method,"ZSSD")==0)
      patchmatch<distance_patch_ZSSD>(u1, u2, w, nnf, cost, minoff, maxoff, iterations, randomtrials, use_horizontal_off);
   else if (strcmp (method,"ZSAD")==0)
      patchmatch<distance_patch_ZSAD>(u1, u2, w, nnf, cost, minoff, maxoff, iterations, randomtrials, use_horizontal_off);
   else if (strcmp (method,"NCC")==0)
      patchmatch<distance_patch_NCC>(u1, u2, w, nnf, cost, minoff, maxoff, iterations, randomtrials, use_horizontal_off);

    // cleanup the interfacing mess I just did (planar to vector)
    if (channels_as_planes) 
      for(int i=0;i<nc*nr*2;i++) { nnf_[i] = nnf[i]; }
    else
      for(int i=0;i<nc*nr;i++) { nnf_[2*i] = nnf[i];  nnf_[2*i+1] = nnf[i+nc*nr];}
    for(int i=0;i<nc*nr;i++)   out_cost_[i] = cost[i];
}



#ifndef DONT_USE_MAIN


extern "C" {
#include "iio.h"
}

/************ IMG IO  **************/

struct Img iio_read_vector_split(char *nm)
{
	struct Img out;
	float *tmpout = iio_read_image_float_split(nm, &out.nx, &out.ny, &out.nch);
	out.data.assign(tmpout,tmpout + out.nx * out.ny * out.nch);
	out.npix = out.nx * out.ny;
	free (tmpout);
	return out;
}



void iio_write_vector_split(char *nm, struct Img &out)
{
	// .front() -> .data() in C++11
	iio_save_image_float_split(nm, &(out.data.front()), out.nx, out.ny, out.nch);
}


void remove_nonfinite_values_Img(struct Img &u, float newval) 
{
   for(int i=0;i<u.npix*u.nch;i++) 
      if (!std::isfinite(u[i])) u[i] = newval; 
}





// c: pointer to original argc
// v: pointer to original argv
// o: option name after hyphen
// d: default value (if NULL, the option takes no argument)
static char *pick_option(int *c, char ***v, char *o, char *d)
{
   int argc = *c;
   char **argv = *v;
   int id = d ? 1 : 0;
   for (int i = 0; i < argc - id; i++)
      if (argv[i][0] == '-' && 0 == strcmp(argv[i] + 1, o)) {
	 char *r = argv[i + id] + 1 - id;
	 *c -= id + 1;
	 for (int j = i; j < argc - id; j++)
	    (*v)[j] = (*v)[j + id + 1];
	 return r;
      }
   return d;
}



/*PATCHMATCH*/



int main(int argc, char* argv[]) 
{

	//read the parameters
   char *in_disp_file = pick_option(&argc, &argv, (char*) "d", (char*) "");
   int dmin = atoi(pick_option(&argc, &argv, (char*) "r", (char*) "0"));
   int dmax = atoi(pick_option(&argc, &argv, (char*) "R", (char*) "30"));
   char* method = pick_option(&argc, &argv, (char*) "t", (char*) "SSD");   //{census|ad|sd|ncc|btad|btsd}
   int iterations = atoi(pick_option(&argc, &argv, (char*) "i", (char*) "5"));
   int w = atoi(pick_option(&argc, &argv, (char*) "w", (char*) "7"));
   bool use_horizontal_off = pick_option(&argc, &argv, (char*) "h", NULL);


	/* patameter parsing - parameters*/
	if(argc<4)
	{
		fprintf (stderr, "too few parameters\n");
		fprintf (stderr, "   usage: %s  [-i iter(5)] [-w window(7)] [-r min_offset(0)] [-R max_offset(10)] [-h] [-d init_disp] u v out [cost [backflow]]\n",argv[0]);
		fprintf (stderr, "        [-t  distance(ad)]: distance = {SSD|SAD|ZSSD|ZSAD|NCC} \n");
		return 1;
	}
	
	int i = 1;
	char* f_u     = (argc>i) ? argv[i] : NULL;      i++;
	char* f_v     = (argc>i) ? argv[i] : NULL;      i++;
	char* f_out   = (argc>i) ? argv[i] : NULL;      i++;
	char* f_cost  = (argc>i) ? argv[i] : NULL;      i++;
	char* f_back  = (argc>i) ? argv[i] : NULL;      i++;
	
	printf("%d %d\n", dmin, dmax);
	
	
	// read input
	struct Img u = iio_read_vector_split(f_u);
	struct Img v = iio_read_vector_split(f_v);

   //remove_nonfinite_values_Img(u, 0);
   //remove_nonfinite_values_Img(v, 0);

	// call
	struct Img ocost(u.nx, u.ny);
   struct Img odisp(u.nx, u.ny, 2);
   odisp.setvalues(dmin);

   if(strcmp (in_disp_file,"")!=0 ){
   	odisp = iio_read_vector_split(in_disp_file);
   }


   if (strcmp (method,"SSD")==0)
      patchmatch<distance_patch_SSD>(u, v, w, odisp, ocost,dmin, dmax,  iterations, 5, use_horizontal_off);
   else if (strcmp (method,"SAD")==0)
      patchmatch<distance_patch_SAD>(u, v, w, odisp, ocost,dmin, dmax,  iterations, 5, use_horizontal_off);
   else if (strcmp (method,"ZSSD")==0)
      patchmatch<distance_patch_ZSSD>(u, v, w, odisp, ocost,dmin, dmax,  iterations, 5, use_horizontal_off);
   else if (strcmp (method,"ZSAD")==0)
      patchmatch<distance_patch_ZSAD>(u, v, w, odisp, ocost,dmin, dmax,  iterations, 5, use_horizontal_off);
   else if (strcmp (method,"NCC")==0)
      patchmatch<distance_patch_NCC>(u, v, w, odisp, ocost,dmin, dmax,  iterations, 5, use_horizontal_off);

	// save the disparity
	
	// generate the backprojected image
	struct Img syn = Img(u.nx, u.ny, u.nch);
	for(int x=0;x<u.nx;x++)
		for(int y=0;y<u.ny;y++){
			Point q = Point(odisp(x,y,0),odisp(x,y,1));
			for(int c=0;c<u.nch;c++)
				if( check_inside_image(x+q[0], y+q[1], v) ) 
					syn(x,y,c) = v(x+q[0],y+q[1],c);
				else 
					syn(x,y,c) = u(x,y,c);
		}
	
	iio_write_vector_split(f_out, odisp);
	if(f_cost) iio_write_vector_split(f_cost, ocost);
	if(f_back) iio_write_vector_split(f_back, syn);
	
	return 0;
}

#endif
