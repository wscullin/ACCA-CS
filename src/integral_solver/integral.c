/*3456789012345678901234567890123456789012345678901234567890123456789012345678*/

#define GAE 114474615732576576.000000


/*******************************************************************
 2015-02-10 intcurve

 Calculate an approximation of the 
 integral of x^8-x^6+x^4-x^2+1 from min to max 

 This is a single file program without a makefile. It is built:
 
 mpicc -lm -g intcurve.c -o intcurve

*******************************************************************/


/*******************************************************************
 
 included headers

*******************************************************************/
#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <omp.h>

typedef struct result_record result_record;

struct result_record
{
char * solver;
double area;
int size;
};


/*******************************************************************
pdiv:
    In: x1 - a double for the first point
        x2 - a double for the second point
        div - a divisor
    Out: return - a double for x1+x2/divisor

*******************************************************************/
double pdiv(double x1, double x2, double div)
{

    return ((x1+x2)/div);
}


/*******************************************************************
line:

   Print a line of 78 = signs to set off program output
   
   In: none
   
   Out: return - none
        stdout - to screen

*******************************************************************/
void line()
{
    int m;
    
    for (m=0; m<=78; m++) printf("=");
    printf("\n");
}


/*******************************************************************
startup_message:

   Print a message at program startup
   
   In: none
   
   Out: return - none
        stdout - to screen

*******************************************************************/
void startup_message()
{
        line();
        printf("\n\tStarting calculation\n\n");
        line();
}


/*******************************************************************
final_output:

   Print the final output
   
   In: the calulated area
   
   Out: return - none
        stdout - to screen

*******************************************************************/
void final_output(result_record * results)
{

    int i;
    double error;
    double offset;

    line();
    for (i=0;i<results[0].size;i++)
    {
        
        error=results[i].area-GAE;
        offset=((GAE-results[i].area)/GAE)*100.0;
        
        printf("\n");
        printf("\tFinal area (%8s):\t%30.6lf\n",results[i].solver,results[i].area);
        printf("\tError (s-A):\t\t%30.6f\n",error);
        printf("\tError (a/A)%%:\t\t%30.10lG\n",offset);
    
    }
    printf("\n");
    line();

}


/*******************************************************************
f:

   Return the value of the polynomial equation:
   f(x)=x^8-x^6+x^4-x^2+1
   
   In: x - a double precision value to evaluate the equation for
   
   Out: return - a double for the value of the equation at x

*******************************************************************/
double f( double x )
{
    return (pow(x,8)-pow(x,6)+pow(x,4)-pow(x,2)+1);
}


/*******************************************************************
pavg:

   Return the average value of an equation at a point
   
   In: f - a double precision function
       x1 - a double for the first point
       x2 - a double for the second point
   
   Out: return - a double for average of the points

*******************************************************************/
double pavg( double(*f)(const double), double x1, double x2 )
{
    return pdiv(f(x1),f(x2),2);
}


/*******************************************************************
trule:

   Return the average value of an equation at a point
   
   In: f - a double precision function
       x1 - a double for the first point
       x2 - a double for the second point
   
   Out: return - a double for average of the points

*******************************************************************/
double trule(double x1, double x2 )
{
    return pavg(&f,x1,x2)*(x2-x1);
}


/*******************************************************************
default_solve:

   Approximate solver
   
   In: xmin - a double for lowest point in range
       xmax - a double for highest point in range
       rank - integer for the rank or thread
       nrank - total number of ranks or threads
       samples - total number of samples
   
   Out: return - a double for average of the points

*******************************************************************/
double default_solve(double xmin, double xmax, int rank, int nranks, long int samples)
{
    long int localsamples=0;
    double ss=0.0;
    double lxmin=0.0;
    double lfxs=0.0;
    double range=0.0;

    
    /* Setup for local evaluation */
    localsamples=samples/nranks;
    range=pdiv(xmax,-xmin,nranks);
    ss=range/localsamples;
    lxmin=(range*rank)+xmin;
     
    /* Caculate rank-local values */
    #pragma omp parallel
    /* loop counter*/
    long int j;
    for (j=0; j<localsamples; j++  )
    {
        lfxs=pavg(&f,lxmin,lxmin+ss)+lfxs;
        lxmin=lxmin+ss;
    }
    
    return (lfxs*((xmax-xmin)/samples));
}


/*******************************************************************
simpson:

   Solve via simpson's method... donuts!
   
   In: xmin - a double for lowest point in range
       xmax - a double for highest point in range
       rank - integer for the rank or thread
       nrank - total number of ranks or threads
       samples - total number of samples
   
   Out: return - a double for average of the points

*******************************************************************/
double simpson( double xmin, double xmax, int rank, int nranks, long int samples)
{
    long int localsamples=0;
    double lxmin=0.00;
    double lxmax=0.00;
    double xdiff=(xmax-xmin);
    double xoffset=0.0;
    double lsect=0.0;
    
   
    
    /* make values evenly divisible */
    do
    {
        localsamples=samples/nranks;
        samples++;
    }
    while((localsamples%2)!=0);
    
    xoffset=xdiff/samples;
    lxmax=(xoffset*localsamples)*(rank+1)+xmin;
    lxmin=(xoffset*localsamples)*(rank)+xmin;
    xoffset=(lxmax-lxmin)/localsamples;
    
    lsect=f(lxmin)+f(lxmax);
    
    #pragma omp parallel
    {
    long int i;
    for(i=1;i<localsamples;i=i+2) 
        lsect=( 4 * f(lxmin+(i*xoffset)))+lsect;

   
    for(i=2;i<localsamples-1;i=i+2) 
        lsect=( 2 * f(lxmin+(i*xoffset)))+lsect;
    }

    lsect=(xoffset/3)*lsect;
    
    return lsect;
}


/*******************************************************************
main:

   Main routine
   
   In: argv and argc - both unused
   
   Out: return - an integer value for the return code

*******************************************************************/
int main( int argc, char *argv[] )
{
    
    /* range across which the integral will be evaluated */
    double xmax=100.3333333;
    double xmin=-10.2666667;
 
    /* set the number of samples used in the approximation */
    long int samples=pow(10,9); 
    
    /* values used in the solve */
    double fdfs=0;
    double ldfs=0;
    double fsimpson=0;
    double lsimpson=0;

    /* loop counters */ 
    int i,j,k;
    
    /* things done to wreck havoc */
    result_record results[2];

    /* variables for process MPI information */ 
    int nranks=1;
    int rank=0;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nranks);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    MPI_Status status;
    MPI_Request request[nranks*2];
    
    /* Status to let us know things have started */
    if (rank == 0)
    {
        startup_message();
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* Setup for local evaluation */
    for (i=0; i<nranks; i++)
    {
        if ( rank == i)
        {
            ldfs=default_solve(xmin,xmax,rank,nranks,samples);
            lsimpson=simpson(xmin,xmax,rank,nranks,samples);
        }           
    }
      
    /* Output local values for each rank */
    for (k=0; k<nranks; k++)
    {
        if ( rank == k)
        {   
            printf("rank %3d deflt: % 26.6f simpson: % 26.6f\n",rank,ldfs,lsimpson);
        }           
    }
    
    MPI_Reduce(&ldfs,&fdfs,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    
    /* Sum for default solver */
    if (rank == 0)
    {

            results[0].solver="default";
            results[0].area=fdfs;
            results[0].size=2;
    }


    MPI_Reduce(&lsimpson,&fsimpson,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    /* Sum for simpson */
    if (rank == 0)
    {
        results[1].solver="simpson";
        results[1].area=fsimpson;
        results[1].size=2;
    }

     
    /* Output final results */    
    if (rank == 0)
    {
        final_output( &results);
    }


    MPI_Finalize();

    return 0;
}
