/* Copyright 2014, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr> */

#ifndef FIMAGE_H
#define FIMAGE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* declarations */
typedef struct fimage {
	int offrow; // Useless
	int offcol; // Useless
	int nrow;
	int ncol;
   int nch;
   int planar; // Set to 1 if the image is stored as rrrrggggbbbb,
               // 0 means that the image is in vector format rgbrgbrgbrgb
               // The image can be converted using fimage_vec_to_planar
   void *meta; // Can be used to link with a deta 
               // structure containing auxiliaru info
               // about this data
               // this pointer is NOT freed by del_fimage
	float *gray;
} *Fimage;

static Fimage new_fimage(void);

static Fimage new_fimage2(int nx, int ny);

static Fimage new_fimage3(float*im, int nx, int ny, int nc);

static Fimage clone_fimage(Fimage i);

static void del_fimage(Fimage i) ;

static Fimage copy_fimage(Fimage dest, Fimage src);

// in place transform from vec to planar form of the array
static void fimage_vec_to_planar( Fimage u) 
{
   int nc  = u->ncol, nr  =u->nrow, nch  =u->nch;
   int x,y,c;
   // copy to a temporary place
   float* tmp = (float*) malloc(nc*nr*nch*sizeof*tmp);
   for(x=0;x<nc*nr*nch;x++) 
      tmp[x] = u->gray[x];
   // rearrange
   for(x=0;x<nc;x++)
   for(y=0;y<nr;y++)
   for(c=0;c<nch;c++)
      u->gray[(x+y*nc) + nc*nr*c] = tmp[(x+y*nc)*nch +c];

   u->planar=1;
   free(tmp);
}

// in place transform from planar to vec form of the array
static void fimage_planar_to_vec( Fimage u) 
{
   int nc  = u->ncol, nr  =u->nrow, nch  =u->nch;
   int x,y,c;
   // copy to a temporary place
   float* tmp = (float*) malloc(nc*nr*nch*sizeof*tmp);
   for(x=0;x<nc*nr*nch;x++) 
      tmp[x] = u->gray[x];
   // rearrange
   for(x=0;x<nc;x++)
   for(y=0;y<nr;y++)
   for(c=0;c<nch;c++)
      u->gray[(x+y*nc)*nch +c] = tmp[(x+y*nc) + nc*nr*c];

   u->planar=0;
   free(tmp);
}

static void fimage_to_vec( Fimage u) 
{
   if(u->planar == 1) fimage_planar_to_vec(u);
}

static void fimage_to_planar( Fimage u) 
{
   if(u->planar == 0) fimage_vec_to_planar(u);
}




/* definitions */
static Fimage new_fimage(void)
{
	Fimage image;

	if ( !(image = (Fimage) calloc(1,sizeof(struct fimage))) ) 
	{
		fprintf(stderr, "[new_fimage] Not enough memory\n");
		exit(1);
		return(NULL);
	}

	image->nrow = image->ncol = image->nch = 0; 
	image->offrow = image->offcol = 0;
	image->gray = NULL;
	image->gray = NULL;
	image->meta = NULL;
	image->planar = 0;
	return(image);
}


static Fimage new_fimage2(int nx, int ny){
	Fimage t= new_fimage();
	t->offrow = t->offcol = 0;
	t->ncol = nx;
	t->nrow = ny;
   t->nch  = 1;
	t->meta = NULL;
	t->planar = 1;
	if ( !(t->gray = (float*) calloc(t->ncol*t->nrow,sizeof(float))))
	{
		fprintf(stderr, "[new_fimage2] Not enough memory\n");
		exit(1);
		return(NULL);
	}
	return (t);
}

static Fimage new_fimage3(float *im, int nx, int ny, int nc){
	Fimage t= new_fimage();
	t->offrow = t->offcol = 0;
	t->ncol = nx;
	t->nrow = ny;
   t->nch  = nc;
	t->meta = NULL;
	t->planar = 1;
   if (im) t->gray = im;
   else {
      if ( !( t->gray = (float*) calloc(t->ncol*t->nrow*t->nch, sizeof(float))))
      {
         fprintf(stderr, "[new_fimage2] Not enough memory\n");
         exit(1);
         return(NULL);
      }
   }
	return (t);
}

static Fimage clone_fimage(Fimage i){
   Fimage t = new_fimage3(NULL, i->ncol, i->nrow, i->nch);
   copy_fimage(t, i);
	return (t);
}


static void del_fimage(Fimage i) {
	if (i->gray) free(i->gray);
	free(i);
}

static Fimage copy_fimage(Fimage dest, Fimage src){
	if ((dest->ncol == src->ncol) && (dest->nrow == src->nrow)){
		memcpy(dest->gray,src->gray,src->ncol*src->nrow*src->nch*sizeof(float)); 
	   dest->meta = src->meta;
	   dest->planar = src->planar;
	} else {
		fprintf(stderr, "[copy_fimage] Image sizes does not match\n");
	}
	return dest;
}

/* access to the pixel in the image u at the position [i (column) ,j (row) ] */
#define _(u,i,j) ((v)->f[ (u)->nch * ( (j) * (u)->ncol + (i)) ] )
/* index of the pixel in the image u, at position [i (column) ,j (row) ] */
#define ii(v,i,j) ( (u)->nch * ((j) * (v)->ncol + (i))  )

// TODO ADAPT THE MACROS TO CONSIDER PLANAR 
/* access to the pixel in the image u at the position [i (column) ,j (row), c (channel) ] */
#define _p(u,i,j,c) (  (u)->gray[ (u)->planar ?  (j) * (u)->ncol + (i) +(u)->ncol* (u)->nrow*(c) : (u)->nch * ( (j) * (u)->ncol + (i)) +(c) ] )
/* index of the pixel in the image u, at position [i (column) ,j (row) ] */
#define iip(u,i,j,c) ( (u)->planar ?  (j) * (u)->ncol + (i) +(u)->ncol* (u)->nrow*(c) : (u)->nch * ((j) * (u)->ncol + (i))   )

#endif
