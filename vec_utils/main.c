#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include "gapmis_vec.h"
#include "../errors.h"

#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)

#define MAX                500
#define LINE_LNG            50

/* Prints the header in the output file */
static void print_header ( FILE * out, const char * filename, const struct gapmis_params* in )
 {
   time_t               t;
   time ( &t );

   fprintf ( out, "####################################\n" );
   fprintf ( out, "# Program: GapMis\n" );
   fprintf ( out, "# Rundate: %s", ctime ( &t ) );
   fprintf ( out, "# Report file: %s\n", filename );
   fprintf ( out, "# Matrix: %s\n", ( in -> scoring_matrix ? "BLOSUM62" : "EDNAFULL" ) );
   fprintf ( out, "# Gap penalty: %.3f\n", in -> gap_open_pen );
   fprintf ( out, "# Extend penalty: %.3f\n", in -> gap_extend_pen );
   fprintf ( out, "####################################\n\n" );
 }

/* Creates seq_gap and mark_mis, and computes min_mis */
static unsigned int print_alignment ( const char * seqa, unsigned int seqa_len, const char * seqb, unsigned int seqb_len, char* seq_gap, char* mark_mis, const struct gapmis_params* in, struct gapmis_align* out )
 {
   unsigned int i, j;

   if ( out -> min_gap > 0 )
    {

      for ( i = 0; i < out -> gap_pos; ++ i )
       {
         seq_gap[i] = seqa[i];
         if ( seqa[i] != seqb[i] )  
          {
            mark_mis[i] = '.';
          }
         else               
           mark_mis[i] = '|';
       }

      for ( j = 0; j < out -> min_gap; ++ j )
       {
         seq_gap[ j + i ] = '-'; 
         mark_mis[ j + i ] = ' ';
       }

      for ( ; i < seqb_len - out -> min_gap && i < seqa_len ; ++ i )
       {
         seq_gap[j + i] = seqa[i];
         if ( seqa[i] != seqb[i + out -> min_gap] ) 
          {
            mark_mis[ j + i ] = '.';
          }
         else
           mark_mis[j + i] = '|';
       }
      
      for ( ; i < seqa_len; ++ i )
       {
         seq_gap[j + i] = seqa[i];
         mark_mis[ j + i ] = '|';
       }
    }
   else
    {
      for ( i = 0; i < seqa_len; ++ i )
       {
         seq_gap[i] = seqa[i];
         if ( seqa[i] != seqb[i] )
          {
            mark_mis[i] = '.';
          }
         else           
           mark_mis[i] = '|';
       }
    }   
   
   return ( 1 );
 }

static void print_line ( const char * s, int start, int stop, int * nr_gaps, int end, int diff, FILE * output, const char* header )
 {
   int                  k;

   if ( start == stop ) return;

   if ( diff )
    {
      fprintf ( output, "%25s", "" );
    }
   else
    {
     if ( header )
       fprintf ( output, "%-13.13s %10d ", header, start + 1 - *nr_gaps );
     else
       fprintf ( output, "%-13.13s %10d ", "", start + 1 - *nr_gaps );
    }

   for ( ; start < stop; ++ start )
    {
      fputc ( s[start], output );
      if ( s[start] == '-' && ! diff ) ++ ( *nr_gaps );
    }

   if ( stop != end )
    {
      for ( k = stop; k < end; ++ k )
       {
         fputc ( ' ', output );
       }
    }
   if ( ! diff )  fprintf ( output, " %-10d", start - *nr_gaps );
   fprintf ( output, "\n" );
 }

/* Wrap two sequences s1 and s2 including the differences (diff) so that the line width is at most len */
static void wrap ( const char * s1, const char * s1_header, const char * s2, const char * s2_header, const char * diff, int len, FILE * output )
 {
   int                  m, n, i, j;
   int                  nr_gaps_a;
   int                  nr_gaps_b;
   int                  nr_lines;

   if ( ! len ) 
     return;

   m = strlen ( s1 );
   n = strlen ( s2 );

   if ( ! n && ! m ) 
     return;

   i         = 0;
   j         = 0;
   nr_gaps_a = 0;
   nr_gaps_b = 0;

   //nr_lines = m / len;
   nr_lines = ( n > m ? m : n ) / len;
   for ( i = 0; i < nr_lines; ++ i )
    {
      /* Sequence s1 */
      print_line ( s1 , i * len, ( i + 1 ) * len, &nr_gaps_a, ( i + 1 ) * len, 0, output, s1_header );

      /* Difference */
      print_line ( diff, i * len, ( i + 1 ) * len, NULL, ( i + 1 ) * len, 1, output, NULL );

      /* Sequence s2 */
      print_line ( s2 , i * len, ( i + 1 ) * len, &nr_gaps_b, ( i + 1 ) * len, 0, output, s2_header );
      fprintf ( output, "\n" );
    }

   /* Last line of first sequence and difference */
   j = i * len;
   if ( j < m || j < n ) 
    {
      print_line ( s1, i * len, min ( m, n ), &nr_gaps_a, ( i + 1 ) * len, 0, output, s1_header );
      print_line ( diff, i * len, ( m < n ) ? m : n, NULL, ( i + 1 ) * len, 1, output, NULL );
      print_line ( s2, i * len, min ( n, m), &nr_gaps_b, ( i + 1 ) * len, 0, output, s2_header );
    }

 }


/* Creates the output file with the many to many alignments */
unsigned int results_many_to_many ( const char * filename, const char const ** p, const char const ** p_header, const char const ** t, const char const ** t_header, const struct gapmis_params* in, struct gapmis_align* out )
 {

   FILE          * output;
   char          * seq_gap;            //the sequence with the inserted gap 
   char          * mark_mis;           //a string with the mismatches marked as '|' (and the matches as ' ')
   unsigned int    n;
   unsigned int    m;
   const char   ** Tmp;
   const char   ** Tmp_header;

 
   if ( ! ( output = fopen ( filename, "w" ) ) )
    {
      errno = IO; 
      return ( 0 );
    }

   for ( ; *p && *p_header; ++ p, ++ p_header )
    {
      Tmp = t;
      Tmp_header = t_header;
      for ( ; *Tmp && *Tmp_header; ++ Tmp, ++ Tmp_header, ++ out )
       {
         n = strlen ( *Tmp );
         m = strlen ( *p );
         if ( m > n )
          {
            errno = LENGTH; //Error: the length of p should be less or equal to the length of t!!!
            return ( 0 );
          }
   
         /* Dynamic memory allocation for seq_gap */
         if ( ! ( seq_gap = ( char * ) calloc ( n + 1 + ( out -> min_gap ) + 1, sizeof ( char ) ) ) )
          {
            errno = MALLOC; //Error: seq_gap could not be allocated!!!
            return ( 0 );
          } 
   
         if ( ! ( mark_mis = ( char* ) calloc ( n + ( out -> min_gap ) + 1, sizeof( char ) ) ) )
          {
            errno = MALLOC; //Error: mark_mis could not be allocated!!!
            return ( 0 );
          } 
         print_header ( output, filename, in );
   
         if ( out -> where == 1 ) //gap is in the text
          {
            print_alignment ( *Tmp, n, *p, m, seq_gap, mark_mis, in, out );
            wrap ( *p, *p_header, seq_gap, *Tmp_header, mark_mis, LINE_LNG, output ); 
          }
         else                     //gap is in the pattern
          {
            print_alignment ( *p, m, *Tmp, n, seq_gap, mark_mis, in, out );
            wrap ( seq_gap, *p_header, *Tmp, *Tmp_header, mark_mis, LINE_LNG, output ); 
          }
   
         fprintf ( output, "\n" );
         fprintf ( output, "Alignment score: %lf\n", out -> max_score );
         fprintf ( output, "Number of mismatches: %d\n", out -> num_mis );
         fprintf ( output, "Length of gap: %d\n", out -> min_gap );
   
         if( out -> min_gap > 0 )
          {
            fprintf ( output, "The gap was inserted in %.13s after position %d\n", ( out -> where == 1 ) ? *Tmp_header : *p_header, out -> gap_pos );
          } 
   
         fprintf ( output, "\n\n" );
         free ( mark_mis );
         free ( seq_gap ); 
    }
  }
  if ( fclose ( output ) ) 
   {
     errno = IO; 
     return ( 0 );
   }
     
   return ( 1 );    
 }


/* used only in main function */
static void print ( const char * label, struct gapmis_align * out )
{
    printf ( "%smax_score: %f\nmin_gap: %d\nwhere: %d\ngap_pos: %d\nnum_mis: %d\n", label, out -> max_score, out -> min_gap, out -> where, out -> gap_pos, out -> num_mis );
}

int main ( int argc, char * argv [] )
 {
   struct gapmis_params         in;
   struct gapmis_align*          out;
   
   FILE                       * fd_query;   /* File descriptors for the reads*/
   FILE                       * fd_target;  /* File descriptors for the reads*/
   
   char                         read[MAX];     /* Buffer for storing the read */
   char                         readId[MAX];   /* Buffer for storing the Id of the read */
   char const                ** querys    = NULL;
   char const                ** querysId  = NULL;
   char const                ** targets   = NULL;
   char const                ** targetsId = NULL;
   int                          max_alloc_target;
   int                          cur_alloc_target;
   int                          max_alloc_query;
   int                          cur_alloc_query;
   int                          i;
   
   /* HERE WE READ THE DATA */
   printf("\nReading the data...\n");

   max_alloc_target = max_alloc_query = 0;
   cur_alloc_target = cur_alloc_query = 0;

   const size_t ALLOC_SIZE = 100;
   
   /* open file 1 (query sequences) */
   if ( ! ( fd_query = fopen ( "queries.fa", "r") ) ) 
    {
      fprintf ( stderr, "Cannot open file queries.fa\n" );
      return ( 0 );
    }

   /* read query sequences */
   while ( ! feof ( fd_query ) )
    {
      if ( fgetc ( fd_query ) && fgets ( readId, MAX, fd_query ) && fgets ( read, MAX, fd_query ) )
       {
          if ( cur_alloc_query >= max_alloc_query )
           {
             querys   = ( char const ** ) realloc ( querys,   ( max_alloc_query + ALLOC_SIZE ) * sizeof ( char * ) ); 
             querysId = ( char const ** ) realloc ( querysId, ( max_alloc_query + ALLOC_SIZE ) * sizeof ( char * ) );
             max_alloc_query += ALLOC_SIZE;
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
      if ( fgetc ( fd_target ) && fgets ( readId, MAX, fd_target ) && fgets ( read, MAX, fd_target ) )
       {
         if ( cur_alloc_target >= max_alloc_target )
          {
            targets   = ( char const ** ) realloc ( targets,   ( max_alloc_target + ALLOC_SIZE ) * sizeof ( char * ) );
            targetsId = ( char const ** ) realloc ( targetsId, ( max_alloc_target + ALLOC_SIZE ) * sizeof ( char * ) );
            max_alloc_target += ALLOC_SIZE;
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
   in . max_gap        = 27; 
   in . scoring_matrix = 0;
   in . gap_open_pen   = 10;
   in . gap_extend_pen = 0.5;

   /* ONE_TO_ONE test */
//   if ( gapmis_one_to_one ( querys[0], targets[0], &in, out ) )
//     {
       /* ONE_TO_ONE result */
//       results_one_to_one ( "ONE_TO_ONE", querys[0], querysId[0], targets[0], targetsId[0], &in, &out[0] );
//     }

   /* ONE_TO_MANY test */
//   if ( gapmis_one_to_many ( querys[0], targets, &in, out ) )
//    {
      /* ONE_TO_MANY result */
//      results_one_to_many ( "ONE_TO_MANY", querys[0], querysId[0], targets, targetsId, &in, out );
//    }
   
   printf("\nAlignment...\n");
   /* MANY_TO_MANY test */
   if ( gapmis_many_to_many ( querys, targets, &in, out ) )
    {
      printf("\nWriting the results...\n");
      /* MANY_TO_MANY results */
      results_many_to_many ( "MANY_TO_MANY", querys, querysId, targets, targetsId, &in, out );
    }

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
