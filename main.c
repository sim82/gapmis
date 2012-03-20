/**
    libgapmis: a library for pairwise sequence aligment with a single gap.
    Copyright (C) 2012 Nikos Alachiotis, Simon Berger, Tomas Flouri, and
    Solon P. Pissis. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <errno.h>
#include <sys/time.h>
#include "gapmis.h"
#include "errors.h"

#define MAIN_MAX 	             500
#define MAIN_ALLOC_SIZE              100


/////////////////////////////////////////// Functions only used for performance measurements ///////////////////////////////////

double gettime( void )
{
    struct timeval ttime;
    gettimeofday(&ttime , 0);
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
};

double get_cells(int n, int m, int b)
{
    double nn = (double) n;
    double mm = (double) m;
    double bb = (double) b;

    if( bb >= nn)
        bb = nn - 1.0;

    return 2*( mm + 1 )*bb + ( mm + 1 ) - bb*( bb + 1 )/2;
}

double total_cells(const char ** txts, const char ** pats, int max_gap)
{
    double cells_total=0.0;
    unsigned int stride=0, strideT=0;
    const char        ** Tmp;
    int i,j;
    for ( Tmp = pats; *Tmp; ++ Tmp, ++ stride );
      for ( Tmp = txts; *Tmp; ++ Tmp, ++ strideT );
   
    for(i=0;i<stride;i++)
        for(j=0;j<strideT;j++)
            cells_total += get_cells(strlen(txts[j]),strlen(pats[i]),max_gap);

    return cells_total;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main ( int argc, char * argv [] )
 {
   struct gapmis_params         in;
   struct gapmis_align*          out;
   
   FILE                       * fd_query;   /* File descriptors for the reads*/
   FILE                       * fd_target;  /* File descriptors for the reads*/
   
   char                         read[MAIN_MAX];     /* Buffer for storing the read */
   char                         readId[MAIN_MAX];   /* Buffer for storing the Id of the read */
   char const                ** querys    = NULL;
   char const                ** querysId  = NULL;
   char const                ** targets   = NULL;
   char const                ** targetsId = NULL;
   int                          max_alloc_target;
   int                          cur_alloc_target;
   int                          max_alloc_query;
   int                          cur_alloc_query;
   int                          i;
   double			start, end;
   double			cells;

   /* HERE WE READ THE DATA */
   max_alloc_target = max_alloc_query = 0;
   cur_alloc_target = cur_alloc_query = 0;

   /* open file 1 (query sequences) */
   if ( ! ( fd_query = fopen ( "queries.fa", "r") ) ) 
    {
      fprintf ( stderr, "Cannot open file queries.fa\n" );
      return ( 0 );
    }

   /* read query sequences */
   while ( ! feof ( fd_query ) )
    {
      if ( fgetc ( fd_query ) && fgets ( readId, MAIN_MAX, fd_query ) && fgets ( read, MAIN_MAX, fd_query ) )
       {
          if ( cur_alloc_query >= max_alloc_query )
           {
             querys   = ( char const ** ) realloc ( querys,   ( max_alloc_query + MAIN_ALLOC_SIZE ) * sizeof ( char * ) ); 
             querysId = ( char const ** ) realloc ( querysId, ( max_alloc_query + MAIN_ALLOC_SIZE ) * sizeof ( char * ) );
             max_alloc_query += MAIN_ALLOC_SIZE;
           }

          read[ strlen ( read ) - 1] = 0;
          readId[ strlen ( readId ) - 1] = 0;
          querys[ cur_alloc_query ]   = strdup ( read );
          querysId[ cur_alloc_query ] = strdup ( readId );

          ++ cur_alloc_query;
       }
    }
   fclose ( fd_query );

   /* open file 2 (target sequences) */
   if ( ! ( fd_target = fopen ( "targets.fa", "r") ) ) 
    {
      fprintf ( stderr, "Cannot open file targets.fa\n" );
      return ( 0 );
    }

   /* read target sequences */
   while ( ! feof ( fd_target ) )
    {
      if ( fgetc ( fd_target ) && fgets ( readId, MAIN_MAX, fd_target ) && fgets ( read, MAIN_MAX, fd_target ) )
       {
         if ( cur_alloc_target >= max_alloc_target )
          {
            targets   = ( char const ** ) realloc ( targets,   ( max_alloc_target + MAIN_ALLOC_SIZE ) * sizeof ( char * ) );
            targetsId = ( char const ** ) realloc ( targetsId, ( max_alloc_target + MAIN_ALLOC_SIZE ) * sizeof ( char * ) );
            max_alloc_target += MAIN_ALLOC_SIZE;
          }
         
         read[ strlen ( read ) - 1] = 0;
         readId[ strlen ( readId ) - 1] = 0;
         targets[ cur_alloc_target ]   = strdup ( read );
         targetsId[ cur_alloc_target ] = strdup ( readId );

         ++ cur_alloc_target;
       }
    }
   fclose ( fd_target );

   /* adjust query and target sizes */
   targets   = ( char const ** ) realloc ( targets, ( cur_alloc_target + 1 ) * sizeof ( char * ) );
   targetsId = ( char const ** ) realloc ( targetsId, ( cur_alloc_target + 1 ) * sizeof ( char * ) );
   targets[ cur_alloc_target ]   = NULL;
   targetsId[ cur_alloc_target ] = NULL;

   querys   = ( char const ** ) realloc ( querys, ( cur_alloc_query + 1 ) * sizeof ( char * ) );
   querysId = ( char const ** ) realloc ( querysId, ( cur_alloc_query + 1 ) * sizeof ( char * ) );
   querys[ cur_alloc_query ]   = NULL;
   querysId[ cur_alloc_query ] = NULL;

   /* allocate the space for the output : | querys | * | targets | */
   out = ( struct gapmis_align * ) calloc ( cur_alloc_target * cur_alloc_query, sizeof ( struct gapmis_align ) );

   /* assign the input parameters */
   in . max_gap        = 30; 
   in . scoring_matrix = 0;
   in . gap_open_pen   = 10;
   in . gap_extend_pen = 0.5;

///////////////////////////////////////////////////////// SEQUENTIAL CALLS ////////////////////////////////////////////////////////////////////////
#if 0
/*
   double			scr;

   if ( gapmis_one_to_one ( querys[0], targets[0], &in, out ) )
     {
       results_one_to_one ( "ONE_TO_ONE", querys[0], querysId[0], targets[0], targetsId[0], &in, &out[0] );
     }

   if ( gapmis_one_to_many ( querys[0], targets, &in, out ) )
    {
      results_one_to_many ( "ONE_TO_MANY", querys[0], querysId[0], targets, targetsId, &in, out );
    }
   
   if ( gapmis_many_to_many ( querys, targets, &in, out ) )
    {
      results_many_to_many ( "MANY_TO_MANY", querys, querysId, targets, targetsId, &in, out );
    }

   if ( ! ( gapmis_one_to_one_scr ( querys[0], targets[0], &in, &scr ) ) )
     {
       printf ("\nErrno:%d", errno);
     }

   if ( ! ( gapmis_one_to_many_opt ( querys[0], targets, &in, out ) ) )
     {
       printf ("\nErrno:%d", errno);
     } */

   start = gettime();
   
   if ( ! ( gapmis_many_to_many_opt ( querys, targets, &in, out ) ) )
     {
       printf ("\nErrno:%d", errno);
     }

   end = gettime();

     /*
   for( i = 0; i < cur_alloc_query; ++i ) {
       printf( "res: %d %f %d\n", i, out[i].max_score, out[i].gap_pos );
       
   }
   */
     
#endif
///////////////////////////////////////////////////////// GPU CALLS /////////////////////////////////////////////////////////////////////////
#ifdef _USE_GPU
/*
   if ( ! ( gapmis_one_to_many_opt_gpu ( querys[0], targets, &in, out ) ) )
     {
       printf ("\nErrno:%d", errno);
     }
*/
   start = gettime();

   if ( ! ( gapmis_many_to_many_opt_gpu ( querys, targets, &in, out ) ) )
     {
       printf ("\nErrno:%d", errno);
     }

   end = gettime();
/* */
#endif
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////// SSE CALLS /////////////////////////////////////////////////////////////////////////
#ifdef _USE_SSE
   gapmis_sse_hint_num_threads(4);
/*
   if ( ! ( gapmis_one_to_many_opt_sse ( querys[0], targets, &in, out ) ) )
     {
       printf ("\nErrno:%d", errno);
     }
*/
   start = gettime();

   if ( ! ( gapmis_many_to_many_opt_sse ( querys, targets, &in, out ) ) )
     {
       printf ("\nErrno:%d", errno);
     }

   end = gettime();
     /*
   for( i = 0; i < cur_alloc_query; ++i ) {
       printf( "res: %d %f %d\n", i, out[i].max_score, out[i].gap_pos );
       
   }
     */
/* */
#endif
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   cells = total_cells( targets, querys, in . max_gap );
   printf("\n Elapsed time :%lf", ( end - start ));
   printf("\n CUPS         :%lf\n", cells / ( end - start ));

   /* Deallocation */
   for ( i = 0; i < cur_alloc_query; ++ i )
    {
      free ( ( void * ) querys[i] );
      free ( ( void * ) querysId[i] );
    }
   free ( querysId );
   free ( querys );

   for ( i = 0; i < cur_alloc_target; ++ i )
    {
      free ( ( void * ) targets[i] );
      free ( ( void * ) targetsId[i] );
    }
   free ( targetsId );
   free ( targets );

   free ( out );

   return ( 0 );
 }
