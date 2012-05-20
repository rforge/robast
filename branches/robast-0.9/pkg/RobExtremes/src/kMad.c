#include <stdlib.h>
#include <stdio.h>

int compare_doubles(const void *a,const void *b)
{
  double *da = (double*)a;
  double *db = (double*)b;
   if (*da > *db)
      return 1;
   else if (*da < *db)
      return -1;
   else
      return 0;
}

void kMad(double *x, int *lx, int *kp, double *d, double *eps)
{
double m;
int n2, r1, r2, i1, i2,j, k = (*kp);

// special cases
    if ((lx[0]) < 1) d[0] = 0;
    if ((lx[0]) == 1) d[0] = x[0];

    if(lx[0] > 1){
       qsort(x, *lx,sizeof(double), compare_doubles);

       if ((lx[0]) % 2 == 0){
          r1 = (lx[0])/2-1;
          r2 = (lx[0])/2;
          m = (x[r1]+x[r2])/2;
		  n2 = lx[0]/2;
       }else{
          r1 = ((lx[0])+1)/2-1;
	      r2 = r1;
          m = x[r1];
		  n2 = (lx[0]+1)/2;
       }
       i1 = r1;
       while( ( x[i1] > m-eps[0] ) && (i1 > 0) ) i1--;
	   i2 = r2;
       while( ( x[i2] < m+eps[0] ) && ( i2 < lx[0]-2) ) i2++;
	

       j = i2 - i1 + 1;
       d[0] = (x[i2]-x[i1])/2;

       while(j < n2)
         {int i1l=0, i2l=0, jl=0, i1r=0, i2r=0, jr=0;
	      double xll=0, xrl=0, xlr=0, xrr=0;
          // l = left check values, r = right check values
        	// check left
	     if(i1>0){
	        i1l = i1 - 1;
	        xll = x[i1l];
	        xrl = m + k * (m-xll);
	        i2l = i2;
	        while((x[i2l] <= xrl +eps[0]) &&(i2l < lx[0]-2)) i2l ++;
	        jl = i2l - i1l + 1;
	     }
	       // check right
	     if(jl< n2){
		    i1 = i1l;
	        i2 = i2l;
	        j = jl;
	        d[0] = m-xll;
	     }else{
		    if(i2< lx[0]-2){
	           i2r = i2 + 1;
	           xrr = x[i2r];
	           xlr = m - (xrr-m) / k;
	           i1r = i1;
	           while((x[i1r] >= xlr -eps[0]) &&(i1r > 0)) i1r --;
	           jr = i2r - i1r + 1;
	        }
	      // decide which is shorter
	        if(jr > jl){
     	       i1 = i1l;
	           i2 = i2l;
	           j = jl;
	           d[0] = m-xll;
	        }else{
               i1 = i1r;
	           i2 = i2r;
	           j = jr;
	           d[0] = m-xlr;
	        }
         }
	   }
	}
}

// int main (){
//   int k = 10;
//   int n = 10;
//   int i;
//   double a[n];
//   double d;
//
//   for (i=0;i<n;i++){
//   a[i] = i*8;
//   }
//
//   kMad(a,&n,&k,&d);
//   printf("hello: %f \n", d);
//   return 0;
//  }
//
