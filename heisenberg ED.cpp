#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define MIN(a,b)   (((a) < (b)) ? (a) : (b))
#define INDEX(i,j) (j*nhil+i)
main()
{
  int Jsign = -1;  
  unsigned int iseed = (unsigned int) time(NULL);
  int N, nhil, nl, ncount, nonzero ;
  int ii, jj, lwork, info, i, j, one = 1 ;
  unsigned int k, *states;
  int *conf, *neigh;
  double beta, uno = 1.0;
  double *H, *W, *work;
  unsigned int search ( unsigned int *, int, unsigned int );
  int btest(unsigned int, int);
  fprintf (stdout,"Number of sites:");
  scanf ("%d", &N) ;
  if ( N <= 0 ) { fprintf (stderr,"wrong N\n"); exit(1); }
  if ( N > 31 ) 
  {
       fprintf (stderr,"N too big, not enough bits in an integer \n"); 
       exit(1);
  }
  nhil = 2^N ;
  conf = (int *) malloc ( N * sizeof (int) );
  neigh= (int *) malloc ( N * sizeof (int) );
  states=(unsigned int *) malloc ( nhil * sizeof (int) );
  H= (double *) malloc ( nhil * nhil * sizeof (double) );
  for ( j = 0; j < N; j++) 
  {
     neigh[j] = j+1 ;
  }
  neigh[N-1] = 0 ;
  ncount=0;
  for ( k = 0; k < nhil; k++ ) 
  {
  	states[ncount++] = k ;
  }
  for ( i = 0; i < nhil*nhil; i++ ) { H[i] = 0.0; }
  nonzero = 0 ;
  for ( ii = 0; ii < nhil; ii++ ) 
  { 	
	   for ( j = 0;  j < N; j++ ) 
	   {
   		conf[j] = -1;
        if ( btest(states[ii],j) ) { conf[j] = 1 ; }
       }
       for ( j=0; j < N; j++ ) 
	   {
         i=neigh[j];
         if ( conf[j] == -1 && conf[i] == 1 ) 
		 {
         	k = states[ii] + 2^j-2^i;
         	jj = search ( states, nhil, k ) ;
         	H[ INDEX(ii,jj) ] = H[ INDEX(ii,jj) ] - 0.5*Jsign ;
            nonzero++ ;
         } 
		 else if ( conf[j] == 1 && conf[i] == -1 ) 
		 {
         	k = states[ii] - 2^j + 2^i;
            jj = search ( states, nhil, k ) ;
            H[ INDEX(ii,jj) ] = H[ INDEX(ii,jj) ] - 0.5*Jsign ;
            nonzero++;
         }
       }
  }
    for ( ii = 0; ii < nhil; ii++) 
	{
       for ( j = 0; j < N; j++ ) 
	   {
          conf[j] = -1;
          if ( btest(states[ii],j) ) conf[j]=1 ;
       }
       for ( j = 0 ; j < N; j++ ) 
	   {
          i=neigh[j];
          H[ INDEX(ii,ii) ] = H[ INDEX(ii,ii) ] - Jsign*0.25*conf[j]*conf[i] ;
       }
       nonzero++;
    }
    lwork = 3*nhil;
    W = (double *) malloc ( nhil * sizeof (double) );
    work = (double *) malloc ( lwork * sizeof (double) );
    dsyev_ ("N","U", &nhil, H, &nhil, W, work, &lwork, &info ) ;
    free (work) ;
    if (info != 0) 
	{
        fprintf(stdout,"Error in dsyev, info=%d\n", info);
    }
    for(i=0;i<10;i++)
	{
	  fprintf(stdout,"%lf",W[i] );
    }
    free (W); free (states); free (H); free (neigh); free (conf);
}
unsigned int search ( unsigned int* p, int nhil, unsigned int k )  
{
	unsigned int imin, imax, lim, l, ii;
    imin = 0;
    imax = nhil-1;
    for ( l = 0; l < nhil; l++ ) 
	{
       ii=(imin+imax)/2 ;
       if (p[ii] == k) 
	   {
          imin=ii ;
          return ii ;
       } 
	   else 
	   {
          if ( p[ii] > k ) 
		  {
             imax=ii-1;
          } 
		  else 
		  {
             imin=ii+1;
          }
       }
    }
    fprintf(stderr, "Something wrong: search not converged\n");
    exit(1);
}
int btest(unsigned int k , int j)
{
    unsigned int uno = 1;
    return (k >> j) & uno ;
}
	
