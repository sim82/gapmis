#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <errno.h>
#include <time.h>
#include "gapmis.h"
#include "errors.h"
#include "EDNAFULL.h"
#include "EBLOSUM62.h"

/* 
 * forward declarations of internal helper functions
 */
static unsigned int dp_algorithm ( int ** G, unsigned int** H, const char* t, unsigned int n, const char * p, unsigned int m, const struct gapmis_params * in );
static unsigned int dp_algorithm_scr ( int ** G, const char* t, unsigned int n, const char * p, unsigned int m, const struct gapmis_params * in );
static int nuc_delta ( char a, char b );
static int pro_delta ( char a, char b );
static unsigned int nuc_char_to_index ( char a );
static unsigned int pro_char_to_index ( char a );
static unsigned int i_limits( unsigned int n, unsigned int m, unsigned int * up, unsigned int * down, unsigned int max_gap );
static unsigned int opt_solution ( int ** G, unsigned int n, unsigned int m, const struct gapmis_params * in, struct gapmis_align * out, unsigned int * start );
static unsigned int opt_solution_scr ( int ** G, unsigned int n, unsigned int m, const struct gapmis_params * in, double * scr );
static unsigned int backtracing ( unsigned int ** H, unsigned int m, unsigned int n, unsigned int start, struct gapmis_align * out );
static double total_scoring( unsigned int gap, double matrix_score, double gap_open_penalty, double gap_extend_penalty );
static unsigned int num_mismatch ( const char * seqa, unsigned int seqa_len, const char * seqb, unsigned int seqb_len, struct gapmis_align * out );


/* Computes the optimal semi-global alignment with the maximum score between each pattern p and all texts t */
unsigned int gapmis_many_to_many_opt ( const char const ** p, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out )
 {
   for ( ; *p; ++ p, ++out )
    {
      if ( ! ( gapmis_one_to_many_opt ( *p, t, in, out ) ) )
        return ( 0 );  
    }

   return ( 1 );
 }

/* Computes the optimal semi-global alignment with the maximum score between a pattern p and all texts t */
unsigned int gapmis_one_to_many_opt ( const char * p, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out )
 {
   double       tmp_scr     = -DBL_MAX;
   double       scr         = 0;
   unsigned int i           = 0;
   unsigned int max_t       = 0;
   const char   ** Tmp      = t;	
	
   for ( ; *Tmp; ++Tmp, ++ i )	//computing the alignment with the maximum score
    {
      if ( ! ( gapmis_one_to_one_scr ( p, *Tmp, in, &scr ) ) )
        return ( 0 );
      if ( scr > tmp_scr )
       {
         max_t = i;
         tmp_scr = scr;
       } 
    } 
   
   if ( ! ( gapmis_one_to_one ( p, t[ max_t ], in, out ) ) ) //computing the rest details of the alignment with the maximum score
     return ( 0 );

   return ( 1 );
 }

/* Computes only the maximum score of the optimal semi-global alignment between pattern p and text t */
unsigned int gapmis_one_to_one_scr ( const char * p, const char * t, const struct gapmis_params * in, double * scr )
 {
   int               ** G; 		//dynamic programming matrix
   unsigned int         n;
   unsigned int         m;
   unsigned int         i;

   /* Checks the input parameters */
   n = strlen ( t );
   m = strlen ( p );
   if ( m > n )
    {
      errno = LENGTH; //Error: the length of p should be less or equal to the length of t!!!
      return ( 0 );
    }
   
   if ( in -> scoring_matrix > 1 )
    {
      errno = MATRIX; //Error: the value of the scoring matrix parameter should be either 0 (nucleotide sequences) or 1 (protein sequences)!!!
      return ( 0 );
    }

   if ( in -> max_gap >= n )
    {
      errno = MAXGAP; //Error: the value of the max gap parameter should be less than the length of t!!!
      return ( 0 );
    }

   /* 2d dynamic memory allocation for matrices G and H*/
   if ( ! ( G = ( int ** ) malloc ( ( n + 1 ) * sizeof ( int * ) ) ) )
    {
      errno = MALLOC; //Error: DP matrix could not be allocated!!!
      return ( 0 );
    } 
   
   if ( ! ( G[0] = ( int * ) calloc ( ( n + 1 ) * ( m + 1 ), sizeof ( int ) ) ) )
    {
      errno = MALLOC; //Error: DP matrix could not be allocated!!!
      return ( 0 );
    } 
   
   for ( i = 1; i < n + 1; ++ i )
     G[i] = ( void * ) G[0] + i * ( m + 1 ) * sizeof ( int );
     
   /* dynamic programming algorithm */
   if ( ! ( dp_algorithm_scr( G, t, n, p, m, in ) ) )
    {
      //Error: dp_algorithm_scr() failed due to bad character!!!
      return ( 0 );	
    }
   
   /* computes the optimal alignment based on the matrix score and the affine gap penalty function */
   opt_solution_scr ( G, n, m, in, scr );
   
   free ( G[0] );
   free ( G );	
   return ( 1 );
   
 } 

#if 1
 
/* Computes the optimal semi-global alignment between each pattern p and all texts t */
unsigned int gapmis_many_to_many ( const char const ** p, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out )
 {
   
   unsigned int         stride = 0;
   const char        ** Tmp;
   for ( Tmp = t; *Tmp; ++ Tmp, ++ stride );  //counting the total number of texts

   for ( ; *p; ++ p, out += stride )
    {
      if ( ! ( gapmis_one_to_many ( *p, t, in, out ) ) )
        return ( 0 );
    }

   return ( 1 );
 }

/* Computes the optimal semi-global alignment between a pattern p and all texts t */
unsigned int gapmis_one_to_many ( const char * p, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out )
 {
   for ( ; *t; ++ t, ++ out )
    {
      if ( ! ( gapmis_one_to_one ( p, *t, in, out ) ) )
        return ( 0 );
    }

   return ( 1 );
 }

/* Computes the optimal semi-global alignment between pattern p and text t */
unsigned int gapmis_one_to_one ( const char * p, const char * t, const struct gapmis_params * in, struct gapmis_align * out )
 {
   int               ** G; 		//dynamic programming matrix
   unsigned int      ** H; 		//backtracing matrix
   unsigned int         n;
   unsigned int         m;
   unsigned int         i;
   unsigned int         start = 0;

   /* Checks the input parameters */
   n = strlen ( t );
   m = strlen ( p );
   if ( m > n )
    {
      errno = LENGTH; //Error: the length of p should be less or equal to the length of t!!!
      return ( 0 );
    }
   
   if ( in -> scoring_matrix > 1 )
    {
      errno = MATRIX; //Error: the value of the scoring matrix parameter should be either 0 (nucleotide sequences) or 1 (protein sequences)!!!
      return ( 0 );
    }

   if ( in -> max_gap >= n )
    {
      errno = MAXGAP; //Error: the value of the max gap parameter should be less than the length of t!!!
      return ( 0 );
    }

   /* 2d dynamic memory allocation for matrices G and H*/
   if ( ! ( G = ( int ** ) malloc ( ( n + 1 ) * sizeof ( int * ) ) ) )
    {
      errno = MALLOC; //Error: DP matrix could not be allocated!!!
      return ( 0 );
    } 
   
   if ( ! ( G[0] = ( int * ) calloc ( ( n + 1 ) * ( m + 1 ), sizeof ( int ) ) ) )
    {
      errno = MALLOC; //Error: DP matrix could not be allocated!!!
      return ( 0 );
    } 
   
   for ( i = 1; i < n + 1; ++ i )
     G[i] = ( void * ) G[0] + i * ( m + 1 ) * sizeof ( int );
   
   if ( ! ( H = ( unsigned int ** ) malloc ( ( n + 1 ) * sizeof ( unsigned int * ) ) ) )
    {
      errno = MALLOC; //Error: DP matrix could not be allocated!!!
      return ( 0 );
    } 
   
   if ( ! ( H[0] = ( unsigned int * ) calloc ( ( n + 1 ) * ( m + 1 ) , sizeof ( unsigned int ) ) ) )
    {
      errno = MALLOC; //Error: DP matrix could not be allocated!!!
      return ( 0 );
    }
   
   for ( i = 1 ; i < n + 1 ; ++ i )
     H[i] = ( void * ) H[0] + i * ( m + 1 ) * sizeof ( unsigned int );
   
   /* dynamic programming algorithm */
   if ( ! ( dp_algorithm( G, H, t, n, p, m, in ) ) )
    {
      //Error: dp_algorithm_scr() failed due to bad character!!!
      return ( 0 );	
    }
   
   /* computes the optimal alignment based on the matrix score and the affine gap penalty function */
   opt_solution ( G, n, m, in, out, &start );
   
   /* computes the position of the gap */
   if ( out -> min_gap > 0 ) 
        backtracing ( H, m, n, start, out );

   /* computes the number of mismatches in the alignment */
   if ( out -> where == 1 )     //gap is the text
     num_mismatch ( t, n, p, m, out );
   else                         //gap is in the pattern or there is no gap
     num_mismatch ( p, m, t, n, out );
   
   free ( G[0] );
   free ( H[0] );
   free ( G );
   free ( H );

   return ( 1 );
   
 } 
#endif
/* The dynamic programming algorithm for calculating matrices G and H */
static unsigned int dp_algorithm ( int ** G, unsigned int ** H, const char * t, unsigned int n, const char * p, unsigned int m, const struct gapmis_params * in )
 {
   int                  gap;
   int                  mis;
   unsigned int         i;
   unsigned int         j;
   int                  matching_score;
   unsigned int 	j_min;
   unsigned int 	j_max;
   unsigned int 	valM;
   unsigned int 	i_max;

   for( i = 0; i < n + 1 ; i++ )      H[i][0] = i;
   for( j = 0; j < m + 1 ; j++ )      H[0][j] = j;

   i_max = min ( n, m + in -> max_gap );

   for( i = 1; i < i_max + 1; i++)
     {
       j_min = max ( 1, (int) ( i - in -> max_gap ));
       j_max = min ( m, (int) ( i + in -> max_gap ));
       for( j = j_min; j <= j_max; j++ )
        {
           matching_score = ( in -> scoring_matrix ? ( int ) pro_delta( t[i - 1], p[j - 1] ) : ( int ) nuc_delta( t[i - 1], p[j - 1] ) );
           if ( matching_score == BADCHAR )
             {  
                errno = BADCHAR; //Error: unrecognizable character!!!
                return ( 0 );
             }	

           mis = G[i - 1][j - 1] + matching_score;
	   gap = G[j][j];
	   valM = i - j;

	   if( j > i )	
	     {
	       gap = G[i][i];
               valM = j - i;
	     }

           if( gap > mis )	H[i][j] = valM;
	   if( i == j )		gap = mis - 1;

           G[i][j] = max ( mis, gap );
         }
      }
   return ( 1 );
 }

/* The dynamic programming algorithm for calculating only matrix G */
static unsigned int dp_algorithm_scr ( int ** G, const char * t, unsigned int n, const char * p, unsigned int m, const struct gapmis_params * in )
 {
   int                  gap;
   int                  mis;
   unsigned int         i;
   unsigned int         j;
   int                  matching_score;
   unsigned int 	j_min;
   unsigned int 	j_max;
   unsigned int 	i_max;

   i_max = min ( n, m + in -> max_gap );

   for( i = 1; i < i_max + 1; i++)
     {
       j_min = max ( 1, (int) ( i - in -> max_gap ));
       j_max = min ( m, (int) ( i + in -> max_gap ));
       for( j = j_min; j <= j_max; j++ )
        {
           matching_score = ( in -> scoring_matrix ? ( int ) pro_delta( t[i - 1], p[j - 1] ) : ( int ) nuc_delta( t[i - 1], p[j - 1] ) );
           if ( matching_score == BADCHAR )
            {  
               errno = BADCHAR; //Error: unrecognizable character!!!
               return ( 0 );
            }	

           mis = G[i - 1][j - 1] + matching_score;
	   gap = G[j][j];
	   if( j > i )		gap = G[i][i];
	   if( i == j )		gap = mis - 1;
           G[i][j] = max ( mis, gap );
         }
      }
   return ( 1 );
 }

/* Returns the score for matching character a and b based on EDNAFULL matrix */
static int nuc_delta ( char a, char b )
 {
   unsigned int index_a = nuc_char_to_index ( a );
   unsigned int index_b = nuc_char_to_index ( b );

   if ( ( index_a < NUC_SCORING_MATRIX_SIZE ) && ( index_b < NUC_SCORING_MATRIX_SIZE ) )
     return ( EDNAFULL_matrix[ index_a ][ index_b ] );
   else //Error: unrecognizable character!!!
     return ( BADCHAR );
 }

/* Returns the score for matching character a and b based on EBLOSUM62 matrix */
static int pro_delta ( char a, char b )
 {
   unsigned int index_a = pro_char_to_index ( a );
   unsigned int index_b = pro_char_to_index ( b );

   if ( ( index_a < PRO_SCORING_MATRIX_SIZE ) && ( index_b < PRO_SCORING_MATRIX_SIZE ) )
     return ( EBLOSUM62_matrix[ index_a ][ index_b ] );
   else //Error: unrecognizable character!!!
     return ( BADCHAR );
 }

/* Returns the index of char a in EDNAFULL matrix */
static unsigned int nuc_char_to_index ( char a )
 {
   unsigned int         index; 

   switch ( a )
    {
      case 'A':
        index = 0; break;

      case 'T':
        index = 1; break;

      case 'G':
        index = 2; break;

      case 'C':
        index = 3; break;

      case 'S':
        index = 4; break;

      case 'W':
        index = 5; break;

      case 'R':
        index = 6; break;

      case 'Y':
        index = 7; break;

      case 'K':
        index = 8; break;

      case 'M':
        index = 9; break;

      case 'B':
        index = 10; break;

      case 'V':
        index = 11; break;

      case 'H':
        index = 12; break;

      case 'D':
        index = 13; break;

      case 'N':
        index = 14; break;

      default:
        index = BADCHAR; break;
    }
   
   return ( index );
 }

/* Returns the index of char a in EBLOSUM62 matrix */
static unsigned int pro_char_to_index ( char a )
 {
   unsigned int         index; 

   switch ( a )
    {
      case 'A':
        index = 0; break;

      case 'R':
        index = 1; break;

      case 'N':
        index = 2; break;

      case 'D':
        index = 3; break;

      case 'C':
        index = 4; break;

      case 'Q':
        index = 5; break;

      case 'E':
        index = 6; break;

      case 'G':
        index = 7; break;

      case 'H':
        index = 8; break;

      case 'I':
        index = 9; break;

      case 'L':
        index = 10; break;

      case 'K':
        index = 11; break;

      case 'M':
        index = 12; break;

      case 'F':
        index = 13; break;

      case 'P':
        index = 14; break;

      case 'S':
        index = 15; break;

      case 'T':
        index = 16; break;

      case 'W':
        index = 17; break;

      case 'Y':
        index = 18; break;

      case 'V':
        index = 19; break;

      case 'B':
        index = 20; break;

      case 'Z':
        index = 21; break;

      case 'X':
        index = 22; break;

      case '*':
        index = 23; break;

      default:
        index = BADCHAR; break;
    }
   return ( index );
 }

/* Computes the optimal alignment using matrix G */
static unsigned int opt_solution ( int ** G, unsigned int n, unsigned int m, const struct gapmis_params * in, struct gapmis_align * out, unsigned int * start )
 {
   unsigned int         i;
   double               score = -DBL_MAX;
   unsigned int         up    = 0;
   unsigned int         down  = 0;
   double               temp_score;
   
   i_limits ( n, m, &up, &down, in -> max_gap );			// computes the i coordinates for matrix G for the last column

   for ( i = up ; i <= down ; i++ )
    {
      temp_score = 0.0;
      if ( i < m )
       {
         if ( m - i <= in -> max_gap )
          {
            temp_score = total_scoring ( m - i, G[i][m], - in -> gap_open_pen, - in -> gap_extend_pen );
            if ( temp_score > score )
             {
               score            = temp_score;
               out -> max_score = score; 
               out -> min_gap   = m - i;
               out -> where     = 1;		//where: gap is in the text and start backtracing from the last column
               ( *start )       = i;		//backtrace from cell G[start,m]
             }
          }
       }
      else if ( i > m )
       {
         if ( i - m <= in -> max_gap )
          {
            temp_score = total_scoring( i - m, G[i][m], - in -> gap_open_pen, - in -> gap_extend_pen );
            if (  temp_score > score )
             {
               score            = temp_score;
               out -> max_score = score; 
               out -> min_gap   = i - m;
               out -> where     = 2;		//where: gap is in the pattern and start backtracing from last column
               ( *start )       = i;		//backtrace from cell G[start,m]
             }
          }
       }
      else if ( i == m )
       {
         temp_score = total_scoring( 0, G[i][m], - in -> gap_open_pen, - in -> gap_extend_pen );
         if (  temp_score > score ) // mgap = 0 
          {
            score            = temp_score;
            out -> max_score = score; 
            out -> min_gap   = 0;
            out -> where     = 0;		//there is no gap
            ( *start )       = m;		//no need to backtrace
          }
       }
    }

   return ( 1 );
 }

/* Computes only the maximum score of the optimal alignment using matrix G */
static unsigned int opt_solution_scr ( int ** G, unsigned int n, unsigned int m, const struct gapmis_params * in, double * scr )
 {
   unsigned int         i;
   double               score = -DBL_MAX;
   unsigned int         up    = 0;
   unsigned int         down  = 0;
   double               temp_score;
   
   i_limits ( n, m, &up, &down, in -> max_gap );			

   for ( i = up ; i <= down ; i++ )
    {
      temp_score = 0.0;
      if ( i < m )
       {
         if ( m - i <= in -> max_gap )
          {
            temp_score = total_scoring ( m - i, G[i][m], - in -> gap_open_pen, - in -> gap_extend_pen );
          }
       }
      else if ( i > m )
       {
         if ( i - m <= in -> max_gap )
          {
            temp_score = total_scoring( i - m, G[i][m], - in -> gap_open_pen, - in -> gap_extend_pen );        
          }
       }
      else if ( i == m )
       {
         temp_score = total_scoring( 0, G[i][m], - in -> gap_open_pen, - in -> gap_extend_pen );
       }
      if (  temp_score > score )
       {
         score     = temp_score;
         ( * scr ) = score;
       }
    }

   return ( 1 );
 }

/* Computes the limits of the i-th coordinate for the matrix G in constant time */
static unsigned int i_limits ( unsigned int n, unsigned int m, unsigned int * up, unsigned int * down, unsigned int max_gap )
 {
   (* up )   = ( ( int ) m - ( int ) max_gap < 0 ) ?  0 : m - max_gap;
   (* down ) = ( m + max_gap > n )                 ?  n : m + max_gap;
   return ( 0 );
 }


/*
Gives the total score of an alignment in constant time
Note: double matrix_score is only the value of G[i][m], i.e. the score of an alignment WITHOUT the affine gap penalty
*/
static double total_scoring( unsigned int gap, double matrix_score, double gap_open_penalty, double gap_extend_penalty )
 {
   return ( matrix_score + ( ( gap > 0 ) ? ( gap - 1 ) * gap_extend_penalty + gap_open_penalty : 0 ) );
 }

/* Computes the position of the gap */
static unsigned int backtracing ( unsigned int ** H, unsigned int m, unsigned int n, unsigned int start, struct gapmis_align * out )
 {
   unsigned int         i, j;

   out -> gap_pos = 0;

   if ( out -> where == 1 || out -> where == 2 )
    {
      i = start; j = m; 	//we start backtracing from the last column
    }
   else
    {
      i = n; j = start;	//we start backtracing from the last row
    }
   while ( i >= 0 && j >= 0 )
    {
      if ( H[i][j] == 0 )
       {
         --i; --j;
       }
      else				
       {
         out -> gap_pos = ( i > j ) ? j : i;
         break;	
       }
    }

   return ( 1 );
 }

/* Computes the number of mismatches between seqa and seqb*/
static unsigned int num_mismatch ( const char * seqa, unsigned int seqa_len, const char * seqb, unsigned int seqb_len, struct gapmis_align * out )
 {
   unsigned int i;
   unsigned int min_mis = 0;

   if ( out -> min_gap > 0 )
    {
      for ( i = 0; i < out -> gap_pos; ++ i )
        if ( seqa[i] != seqb[i] ) ++ min_mis;

      for ( ; i < seqb_len - out -> min_gap && i < seqa_len ; ++ i )
        if ( seqa[i] != seqb[i + out -> min_gap] ) ++ min_mis;
    }
   else
    {
      for ( i = 0; i < seqa_len; ++ i )
        if ( seqa[i] != seqb[i] ) ++ min_mis;
    }	
   out -> num_mis = min_mis;

   return ( 1 );
 }

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

/* Creates the output file with the one to one alignment */
unsigned int results_one_to_one ( const char * filename, const char * p, const char * p_header, const char * t, const char * t_header, const struct gapmis_params* in, struct gapmis_align* out )
 {

   FILE          * output;
   char          * seq_gap;            //the sequence with the inserted gap 
   char          * mark_mis;           //a string with the mismatches marked as '|' (and the matches as ' ')
   unsigned int    n;
   unsigned int    m;
 
   n = strlen ( t );
   m = strlen ( p );
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
   
   if ( ! ( output = fopen ( filename, "w" ) ) )
    {
      errno = IO; 
      return ( 0 );
    }
   
   print_header ( output, filename, in );
   
   if ( out -> where == 1 ) //gap is in the text
    {
      print_alignment ( t, n, p, m, seq_gap, mark_mis, in, out );
      wrap ( p, p_header, seq_gap, t_header, mark_mis, LINE_LNG, output ); 
    }
   else                     //gap is in the pattern
    {
      print_alignment ( p, m, t, n, seq_gap, mark_mis, in, out );
      wrap ( seq_gap, p_header, t, t_header, mark_mis, LINE_LNG, output ); 
    }
   
   fprintf ( output, "\n" );
   fprintf ( output, "Alignment score: %lf\n", out -> max_score );
   fprintf ( output, "Number of mismatches: %d\n", out -> num_mis );
   fprintf ( output, "Length of gap: %d\n", out -> min_gap );
   
   if( out -> min_gap > 0 )
    {
      fprintf ( output, "The gap was inserted in %.13s after position %d\n", ( out -> where == 1 ) ? t_header : p_header, out -> gap_pos );
    } 
   
   fprintf ( output, "\n\n" );
   
   if ( fclose ( output ) ) 
    {
      errno = IO; 
      return ( 0 );
    }
   
   free ( mark_mis );
   free ( seq_gap );
   
   return ( 1 );	
 }

/* Creates the output file with the one to many alignments */
unsigned int results_one_to_many ( const char * filename, const char * p, const char const * p_header, const char const ** t, const char const ** t_header, const struct gapmis_params* in, struct gapmis_align* out )
 {

   FILE          * output;
   char          * seq_gap;            //the sequence with the inserted gap 
   char          * mark_mis;           //a string with the mismatches marked as '|' (and the matches as ' ')
   unsigned int    n;
   unsigned int    m;
  
   m = strlen ( p );
   if ( ! ( output = fopen ( filename, "w" ) ) )
    {
      errno = IO; 
      return ( 0 );
    }
   
   for ( ; *t && *t_header; ++ t, ++ t_header, ++ out )
    {
       n = strlen ( *t );
       
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
         print_alignment ( *t, n, p, m, seq_gap, mark_mis, in, out );
         wrap ( p, p_header, seq_gap, *t_header, mark_mis, LINE_LNG, output ); 
       }
      else                     //gap is in the pattern
       {
         print_alignment ( p, m, *t, n, seq_gap, mark_mis, in, out );
         wrap ( seq_gap, p_header, *t, *t_header, mark_mis, LINE_LNG, output ); 
       }
   
      fprintf ( output, "\n" );
      fprintf ( output, "Alignment score: %lf\n", out -> max_score );
      fprintf ( output, "Number of mismatches: %d\n", out -> num_mis );
      fprintf ( output, "Length of gap: %d\n", out -> min_gap );
   
      if( out -> min_gap > 0 )
       {
         fprintf ( output, "The gap was inserted in %.13s after position %d\n", ( out -> where == 1 ) ? *t_header : p_header, out -> gap_pos );
       } 
   
      fprintf ( output, "\n\n" );
      free ( mark_mis );
      free ( seq_gap ); 
    }

   if ( fclose ( output ) ) 
    {
      errno = IO; 
      return ( 0 );
    }
     
   return ( 1 );	
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


#ifdef _USE_GPU
/* 
 * forward declarations of internal helper functions
 */
static cl_platform_id get_gpu_id(int * error);
static cl_device_id get_dev_id(cl_platform_id gpu_id, int * error);
static cl_context create_context(cl_device_id dev_id, int * error);
static cl_command_queue create_cmd_queue (cl_device_id dev_id, cl_context context, int * error);
static cl_kernel load_kernel (char * name, char * kernel_name, cl_device_id dev_id, cl_context context, int * error);
static cl_mem malloc_device (cl_context context, size_t size, int * error);
static void init_device_mem_int (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, int * mem, size_t size, int * error);
static void init_device_mem_uint (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, unsigned int * mem, size_t size, int * error);
static void init_device_mem_float (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, float * mem, size_t size, int * error);
static void set_kernel_arguments (cl_kernel kernel, cl_command_queue cmd_queue, cl_mem cl_mem0, cl_mem cl_mem1, cl_mem cl_mem2, cl_mem cl_mem3, cl_mem cl_mem4, cl_mem cl_mem5, cl_mem cl_mem6, cl_mem cl_mem7);
static void read_device_mem_float (cl_command_queue cmd_queue, size_t size, float * mem, cl_mem dev_mem, int * error);
static void fill_txtsLenVec (unsigned int totalTxts, const char const ** in, int * outVec);
static void fill_argsVec (unsigned int totalPats, unsigned int totalTxts, const char const ** in, unsigned int max_gap, unsigned int pBlockSize, unsigned int maxPatLen, unsigned int maxTxtLen, int * outVec);
static void fill_patsVec (unsigned int total, unsigned int blockSize, const char const ** in, unsigned int * outVec, int matrix);
static void fill_txtsVec (unsigned int total, unsigned int blockSize, const char const ** in, unsigned int * outVec, int matrix);
static unsigned int get_pblock_size (unsigned int input, unsigned int mult);
static unsigned int get_max_length (unsigned int total, const char ** input);
static unsigned int get_min_length (unsigned int total, const char ** input);
static unsigned int get_number_of_sequences (const char ** input);
static unsigned int get_number_of_groups (int elements, int groupSize);
static unsigned int kernel_launch (cl_kernel kernel, cl_context context, cl_command_queue cmd_queue, const char const ** p, const char const ** t, const struct gapmis_params * in, float * scores);
static unsigned int kernel_launch_l (cl_kernel kernel, cl_context context, cl_command_queue cmd_queue, const char const ** p, const char const ** t, const struct gapmis_params * in, float * scores, struct gapmis_align * out);
static void update_group_match (float * groupScores, int * groupMatch, float * groupMatchScores, unsigned int patGroupSize, unsigned int txtGroupSize, int pats, int txts, int patGroupIndex, int txtGroupIndex);
static void set_invalid(int * groupMatch, int groupSize);
static void set_minimum(float * groupMatchScores, int groupSize);
static void set_null (const char ** input, int size);
static void initialize_pointers (const char * groupPatterns[], int groupIndex, int groupSize, const char const ** source, int sourceSize);

unsigned int gapmis_one_to_many_opt_gpu ( const char * p1, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out )
{
	const char * p[] = { p1, NULL};

	if ( in -> scoring_matrix > 1 )
	{
		errno = MATRIX;
		return ( 0 );
	}

	unsigned int pats = get_number_of_sequences (p);
	unsigned int txts = get_number_of_sequences (t);
	unsigned int	maxPatLen = get_max_length (pats, p);
	unsigned int	minTxtLen = get_min_length (txts, t);

	if(maxPatLen > minTxtLen)
	{
		errno = LENGTH;
      		return ( 0 );
	}

	if ( in -> max_gap >= minTxtLen )
	{
		errno = MAXGAP; 
		return ( 0 );
	}


	int err = -1;

	cl_platform_id gpu_id = get_gpu_id(&err);
	
	if(err)
	{	
	 	errno = NOGPU;
      		return ( 0 );
	}

	cl_device_id dev_id = get_dev_id(gpu_id, &err);

	if(err)
	{	
	 	errno = NOGPU;
      		return ( 0 );
	}

	cl_context context = create_context(dev_id, &err);
	if(err)
	{	
	 	errno = GPUERROR;
      		return ( 0 );
	}

	cl_command_queue cmd_queue = create_cmd_queue (dev_id, context, &err);
	if(err)
	{	
	 	errno = GPUERROR;
      		return ( 0 );
	}

	cl_kernel kernel;

	if(in->scoring_matrix==0)
		kernel = load_kernel ("kernel_dna.cl", "gapmis_kernel", dev_id, context, &err);
	else
		kernel = load_kernel ("kernel_pro.cl", "gapmis_kernel", dev_id, context, &err);

	if(err)
	{	
	 	errno = KERNEL;
      		return ( 0 );
	}

	const unsigned int patGroupSize = 1;
	const unsigned int txtGroupSize = 768;
	unsigned int i, j;	
	unsigned int patGroups = get_number_of_groups (pats, patGroupSize);
	unsigned int txtGroups = get_number_of_groups (txts, txtGroupSize);	

	const char * groupPatterns[patGroupSize+1];
	set_null (groupPatterns, patGroupSize+1);

	const char * groupTexts[txtGroupSize+1];
	set_null (groupTexts, txtGroupSize+1);

	float * groupScores;
        groupScores = calloc (patGroupSize*txtGroupSize, sizeof(float) );

	int groupMatch [patGroupSize];
	float groupMatchScores [patGroupSize];
	set_invalid(groupMatch,patGroupSize);
	set_minimum(groupMatchScores,patGroupSize);	

	for(i=0;i<patGroups;i++)
	{
		set_null (groupPatterns, patGroupSize+1);
		initialize_pointers (groupPatterns,i,patGroupSize,p,pats);
		set_invalid(groupMatch,patGroupSize);
		set_minimum(groupMatchScores,patGroupSize);
		
		for(j=0;j<txtGroups;j++)
		{			
			set_null (groupTexts, txtGroupSize+1);
			initialize_pointers (groupTexts,j,txtGroupSize,t,txts);

			if(kernel_launch (kernel, context, cmd_queue, groupPatterns, groupTexts, in, groupScores))
				return (0);			

			update_group_match (groupScores,groupMatch,groupMatchScores,patGroupSize,txtGroupSize, pats, txts, i, j);
		
		}

		for(j=0;j<patGroupSize;j++)
		{
			if(i*patGroupSize+j<pats)
			{
				groupPatterns[0] = p[i*patGroupSize+j];
				groupPatterns[1] = NULL;

				groupTexts[0] = t[groupMatch[j]];
				groupTexts[1] = NULL;
				
				if(kernel_launch_l (kernel, context, cmd_queue, groupPatterns, groupTexts, in, groupScores,&out[i*patGroupSize+j]))
					return (0);				
			}
		}
	}

        free ( groupScores );
        clReleaseContext ( context );
	clReleaseCommandQueue ( cmd_queue );
        clReleaseKernel(kernel);

	return ( 1 );
 }

unsigned int gapmis_many_to_many_opt_gpu ( const char const ** p, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out )
{

	if ( in -> scoring_matrix > 1 )
	{
		errno = MATRIX;
		return ( 0 );
	}

	unsigned int pats = get_number_of_sequences (p);
	unsigned int txts = get_number_of_sequences (t);
	unsigned int	maxPatLen = get_max_length (pats, p);
	unsigned int	minTxtLen = get_min_length (txts, t);

	if(maxPatLen > minTxtLen)
	{
		errno = LENGTH;
      		return ( 0 );
	}

	if ( in -> max_gap >= minTxtLen )
	{
		errno = MAXGAP; 
		return ( 0 );
	}


	int err = -1;

	cl_platform_id  gpu_id = get_gpu_id(&err);
	
	if(err)
	{	
	 	errno = NOGPU;
      		return ( 0 );
	}

	cl_device_id  dev_id = get_dev_id(gpu_id, &err);

	if(err)
	{	
	 	errno = NOGPU;
      		return ( 0 );
	}

	cl_context context = create_context(dev_id, &err);
	if(err)
	{	
	 	errno = GPUERROR;
      		return ( 0 );
	}

	cl_command_queue cmd_queue = create_cmd_queue (dev_id, context, &err);
	if(err)
	{	
	 	errno = GPUERROR;
      		return ( 0 );
	}

	cl_kernel kernel;

	if(in->scoring_matrix==0)
		kernel = load_kernel ("kernel_dna.cl", "gapmis_kernel", dev_id, context, &err);
	else
		kernel = load_kernel ("kernel_pro.cl", "gapmis_kernel", dev_id, context, &err);

	if(err)
	{	
	 	errno = KERNEL;
      		return ( 0 );
	}

	const unsigned int patGroupSize = 1024;
	const unsigned int txtGroupSize = 32;
	unsigned int i, j;	
	unsigned int patGroups = get_number_of_groups (pats, patGroupSize);
	unsigned int txtGroups = get_number_of_groups (txts, txtGroupSize);	

	const char * groupPatterns[patGroupSize+1];
	set_null (groupPatterns, patGroupSize+1);

	const char * groupTexts[txtGroupSize+1];
	set_null (groupTexts, txtGroupSize+1);

	float * groupScores;
        groupScores = calloc (patGroupSize*txtGroupSize, sizeof(float) );

	int groupMatch [patGroupSize];
	float groupMatchScores [patGroupSize];
	set_invalid(groupMatch,patGroupSize);
	set_minimum(groupMatchScores,patGroupSize);	
      
	for(i=0;i<patGroups;i++)
	{
		set_null (groupPatterns, patGroupSize+1);
		initialize_pointers (groupPatterns,i,patGroupSize,p,pats);
		set_invalid(groupMatch,patGroupSize);
		set_minimum(groupMatchScores,patGroupSize);
		
		for(j=0;j<txtGroups;j++)
		{			
			set_null (groupTexts, txtGroupSize+1);
			initialize_pointers (groupTexts,j,txtGroupSize,t,txts);

			if(kernel_launch (kernel, context, cmd_queue, groupPatterns, groupTexts, in, groupScores))
				return (0);			
			
			update_group_match (groupScores,groupMatch,groupMatchScores,patGroupSize,txtGroupSize, pats, txts, i, j);
		
		}

		for(j=0;j<patGroupSize;j++)
		{
			if(i*patGroupSize+j<pats)
			{
				groupPatterns[0] = p[i*patGroupSize+j];
				groupPatterns[1] = NULL;

				groupTexts[0] = t[groupMatch[j]];
				groupTexts[1] = NULL;
				
				if(kernel_launch_l (kernel, context, cmd_queue, groupPatterns, groupTexts, in, groupScores,&out[i*patGroupSize+j]))
					return (0);				
			}
		}
	}

        free ( groupScores );
        clReleaseContext ( context );
	clReleaseCommandQueue ( cmd_queue );
        clReleaseKernel(kernel);
	
        return ( 1 );
 }


static void kernel_one_to_one_dna (unsigned int groupID, unsigned int localID, unsigned int * patsVec, unsigned int * txtsVec, int * argsVec, int * txtsLenVec, float * pensVec, int * hproVec, int * dproVec, float * scrsVec, unsigned int ** H, struct gapmis_align * out, int * start)
{
	int i, j, cur_diag_nxt, mis, gap,
	    max_gap = argsVec[2], 
	    blockSize = argsVec[3],
	    maxPatLen = argsVec[4],
	    dproVecGsize = argsVec[5],
	    hproVecGsize = argsVec[6],
	    m = argsVec[groupID + 7],
	    n = txtsLenVec[localID],
	    doffset = 1,gapmismax;

	unsigned int j_min, j_max, abs_ij, patChar, txtChar;

	float temp_score, score = -1000000.0;		
	
	unsigned int min_gap, where;

	for( i = 0; i < n + 1 ; i++ )      H[i][0] = i;
	for( j = 0; j < m + 1 ; j++ )      H[0][j] = j;

	if(max ( 1,  m - max_gap )==1)
	{
		score = ( m - 1 ) * pensVec[1] + pensVec[0];		 
               	min_gap   = m;
               	where     = 1;		
              	(*start)       = 0;	
	}	

	for( i = 1; i <= m; i++ )
	{
		patChar = patsVec [groupID * maxPatLen + i - 1];

		j_min = max ( 1,  i - max_gap );

		j_max = min ( n,  i + max_gap );

		if(i<= max_gap+1)
		{
			cur_diag_nxt = 0;
		}
		else
		{
			cur_diag_nxt = hproVec[hproVecGsize*groupID + doffset*blockSize + localID];		
			doffset++;		
		}

		for( j = j_min; j <= j_max; j++)
		{
			txtChar = txtsVec[(j-1)*blockSize + localID]; 

			mis = cur_diag_nxt + EDNAFULL_matrix [txtChar][patChar];
			
			if(j<i)
			{	
				abs_ij = i-j;
				gap = dproVec[dproVecGsize*groupID + j*blockSize + localID];
			}
			else
			{	
				abs_ij = j-i;
				gap = dproVec[dproVecGsize*groupID + i*blockSize + localID]; 
			}

			if(i==j)
			{
				gap = mis - 1;
				dproVec[dproVecGsize*groupID + j*blockSize + localID]= mis;
			}
	
			if( gap > mis )	H[j][i] = abs_ij;

			cur_diag_nxt = hproVec[hproVecGsize*groupID + j*blockSize + localID];
	
			gapmismax = max ( mis, gap );
		
			hproVec[hproVecGsize*groupID + j*blockSize + localID] = gapmismax;

			if(i==m)
			{				
				if(abs_ij<=max_gap)
				{
					temp_score = (float)gapmismax;			
				
					if(abs_ij>0)
						temp_score += ( abs_ij - 1 ) * pensVec[1] + pensVec[0];			

					if(temp_score>score)
					{
						score = temp_score;

						if(i>j)
						{							
               						 min_gap   = m-j;
               						 where     = 1;		
              						 (*start)  = j;		
						}

						if(i==j)
						{
							
            						min_gap   = 0;
            						where     = 0;		
            						(*start)  = m;		
						}

						if(i<j)
						{						  
						       min_gap   = j-m;
						       where     = 2;		
						       ( *start )   = j;		
						}								
					}
				}
			}								
		}
	}

	scrsVec[groupID*blockSize+localID]=score;
	out->max_score = score;
	out->min_gap = min_gap;
	out->where = where;	
}

static void kernel_one_to_one_pro (unsigned int groupID, unsigned int localID, unsigned int * patsVec, unsigned int * txtsVec, int * argsVec, int * txtsLenVec, float * pensVec, int * hproVec, int * dproVec, float * scrsVec, unsigned int ** H, struct gapmis_align * out, int * start)
{
	int i, j, cur_diag_nxt, mis, gap,
	    max_gap = argsVec[2], 
	    blockSize = argsVec[3],
	    maxPatLen = argsVec[4],
	    dproVecGsize = argsVec[5],
	    hproVecGsize = argsVec[6],
	    m = argsVec[groupID + 7],
	    n = txtsLenVec[localID],
	    doffset = 1,gapmismax;

	unsigned int j_min, j_max, abs_ij, patChar, txtChar;

	float temp_score, score = -1000000.0;		
	
	unsigned int min_gap, where;

	for( i = 0; i < n + 1 ; i++ )      H[i][0] = i;
	for( j = 0; j < m + 1 ; j++ )      H[0][j] = j;

	if(max ( 1,  m - max_gap )==1)
	{
		score = ( m - 1 ) * pensVec[1] + pensVec[0];		 
               	min_gap   = m;
               	where     = 1;		
              	(*start)       = 0;	
	}	

	for( i = 1; i <= m; i++ )
	{
		patChar = patsVec [groupID * maxPatLen + i - 1];

		j_min = max ( 1,  i - max_gap );

		j_max = min ( n,  i + max_gap );

		if(i<= max_gap+1)
		{
			cur_diag_nxt = 0;
		}
		else
		{
			cur_diag_nxt = hproVec[hproVecGsize*groupID + doffset*blockSize + localID];		
			doffset++;		
		}

		for( j = j_min; j <= j_max; j++)
		{
			txtChar = txtsVec[(j-1)*blockSize + localID]; 

			mis = cur_diag_nxt + EBLOSUM62_matrix [txtChar][patChar];
			
			if(j<i)
			{	
				abs_ij = i-j;
				gap = dproVec[dproVecGsize*groupID + j*blockSize + localID];
			}
			else
			{	
				abs_ij = j-i;
				gap = dproVec[dproVecGsize*groupID + i*blockSize + localID]; 
			}

			if(i==j)
			{
				gap = mis - 1;
				dproVec[dproVecGsize*groupID + j*blockSize + localID]= mis;
			}
	
			if( gap > mis )	H[j][i] = abs_ij;

			cur_diag_nxt = hproVec[hproVecGsize*groupID + j*blockSize + localID];
	
			gapmismax = max ( mis, gap );
		
			hproVec[hproVecGsize*groupID + j*blockSize + localID] = gapmismax;

			if(i==m)
			{				
				if(abs_ij<=max_gap)
				{
					temp_score = (float)gapmismax;				
		
				
					if(abs_ij>0)
						temp_score += ( abs_ij - 1 ) * pensVec[1] + pensVec[0];			

					if(temp_score>score)
					{
						score = temp_score;

						if(i>j)
						{							
               						 min_gap   = m-j;
               						 where     = 1;		
              						 (*start)  = j;		
						}

						if(i==j)
						{							
            						min_gap   = 0;
            						where     = 0;		
            						(*start)  = m;		
						}

						if(i<j)
						{						  
						       min_gap   = j-m;
						       where     = 2;		
						       ( *start ) = j;		
						}								
					}
				}
			}								
		}
	}

	scrsVec[groupID*blockSize+localID]=score;
	out->max_score = score;
	out->min_gap = min_gap;
	out->where = where;
}

static unsigned int kernel_launch (cl_kernel kernel, cl_context context, cl_command_queue cmd_queue, const char const ** p, const char const ** t, const struct gapmis_params * in, float * scores)
{
        int error=1;
	unsigned int	pats = get_number_of_sequences (p);
	unsigned int	txts = get_number_of_sequences (t);

	unsigned int	maxTxtLen = get_max_length (txts, t);
	unsigned int	maxPatLen = get_max_length (pats, p);
	unsigned int	pBlockSize = get_pblock_size (txts,32); 
	unsigned int 	hproVecLen = pats * pBlockSize * (maxTxtLen + 1);
	unsigned int 	dproVecLen = pats * pBlockSize * (maxPatLen + 1);
	
	unsigned int * txtsVec = calloc(maxTxtLen*pBlockSize, sizeof(unsigned int));
	unsigned int * patsVec = calloc(maxPatLen*pats, sizeof(unsigned int));
	         int * argsVec = malloc(sizeof(int)*(pats+7));		
		 int * txtsLenVec = calloc(pBlockSize,sizeof(int));
	       float * pensVec = malloc(sizeof(float)*2);
	         int * hproVec = calloc(hproVecLen,sizeof(int));
	         int * dproVec = calloc(dproVecLen,sizeof(int));

	cl_int err;	

	if(patsVec==NULL   || txtsVec==NULL      || argsVec==NULL || 
           pensVec == NULL || txtsLenVec == NULL || hproVec==NULL || 
	   dproVec==NULL)
	{	
	 	errno = MALLOC;
      		return ( 1 );
	}

	fill_txtsVec (txts, pBlockSize, t, txtsVec, in->scoring_matrix);
	fill_patsVec (pats, maxPatLen, p, patsVec, in->scoring_matrix);
	fill_argsVec (pats, txts, p, in->max_gap, pBlockSize, maxPatLen, maxTxtLen, argsVec);
	fill_txtsLenVec (txts, t, txtsLenVec);

	pensVec[0] = - in -> gap_open_pen;
	pensVec[1] = - in -> gap_extend_pen;	
	

	cl_mem txtsVec_device = malloc_device (context, (maxTxtLen*pBlockSize)*sizeof(unsigned int), &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 1 );	
	}
	
	init_device_mem_uint (context, cmd_queue, txtsVec_device, txtsVec, maxTxtLen*pBlockSize, &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 1 );	
	}
	
	cl_mem patsVec_device = malloc_device (context, (maxPatLen*pats)*sizeof(unsigned int), &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 1 );	
	}

	init_device_mem_uint (context, cmd_queue, patsVec_device, patsVec,maxPatLen*pats, &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 1 );	
	}

	cl_mem argsVec_device = malloc_device (context, (pats+7)*sizeof(int), &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 1 );	
	}

	init_device_mem_int (context, cmd_queue, argsVec_device, argsVec, pats+7, &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 1 );	
	}

	cl_mem txtsLenVec_device = malloc_device (context, pBlockSize*sizeof(int), &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 1 );	
	}

	init_device_mem_int (context, cmd_queue, txtsLenVec_device, txtsLenVec, pBlockSize, &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 1 );	
	}

	cl_mem pensVec_device = malloc_device (context, 2*sizeof(float), &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 1 );	
	}

	init_device_mem_float (context, cmd_queue, pensVec_device, pensVec, 2, &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 1 );	
	}

	cl_mem hproVec_device = malloc_device (context, hproVecLen*sizeof(int), &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 1 );	
	}

	init_device_mem_int (context, cmd_queue, hproVec_device, hproVec, hproVecLen, &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 1 );	
	}

	cl_mem dproVec_device = malloc_device (context, dproVecLen*sizeof(int), &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 1 );	
	}

	init_device_mem_int (context, cmd_queue, dproVec_device, dproVec, dproVecLen, &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 1 );	
	}

	cl_mem scrsVec_device = malloc_device (context, (pats*pBlockSize)*sizeof(float), &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 1 );	
	}

	err = clFinish(cmd_queue);
	if(error)
	{
	 	errno = GPUERROR;
      		return ( 1 );	
	}

	set_kernel_arguments (kernel, cmd_queue, patsVec_device, txtsVec_device, argsVec_device, txtsLenVec_device, pensVec_device, hproVec_device, dproVec_device, scrsVec_device);

	err = clFinish(cmd_queue);
	if(error)
	{
	 	errno = GPUERROR;
      		return ( 1 );	
	}	

	size_t WorkSizeGlobal[] = {pBlockSize * pats};
	size_t WorkSizeLocal[] = {pBlockSize};

	err = clEnqueueNDRangeKernel(cmd_queue, kernel, 1, NULL, WorkSizeGlobal, WorkSizeLocal, 0, NULL, NULL);
	if(error)
	{
	 	errno = KERNEL;
      		return ( 1 );	
	}

	err=clFinish(cmd_queue);
	if(error)
	{
	 	errno = GPUERROR;
      		return ( 1 );	
	}	

	read_device_mem_float (cmd_queue, pats*pBlockSize, scores, scrsVec_device, &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 1 );	
	}

	free (txtsVec);
	free (patsVec);
	free (argsVec);
	free (txtsLenVec);
	free (pensVec);
	free (hproVec);
	free (dproVec);

	clReleaseMemObject(patsVec_device);
	clReleaseMemObject(txtsVec_device);
	clReleaseMemObject(argsVec_device);
	clReleaseMemObject(txtsLenVec_device);
	clReleaseMemObject(pensVec_device);
	clReleaseMemObject(hproVec_device);
	clReleaseMemObject(dproVec_device);
	clReleaseMemObject(scrsVec_device);

	return (0);
}

static unsigned int kernel_launch_l (cl_kernel kernel, cl_context context, cl_command_queue cmd_queue, const char const ** p, const char const ** t, const struct gapmis_params * in, float * scores, struct gapmis_align * out)
{
	unsigned int 	i;
		 int	start;
	unsigned int	pats = 1;
	unsigned int	txts = 1;

	unsigned int	maxTxtLen = strlen(t[0]);
	unsigned int	maxPatLen = strlen(p[0]);
	unsigned int	pBlockSize = 1;
	unsigned int 	hproVecLen = pats * pBlockSize * (maxTxtLen + 1);
	unsigned int 	dproVecLen = pats * pBlockSize * (maxPatLen + 1);
	
	unsigned int * txtsVec = calloc(maxTxtLen*pBlockSize, sizeof(unsigned int));
	unsigned int * patsVec = calloc(maxPatLen*pats, sizeof(unsigned int));
	         int * argsVec = malloc(sizeof(int)*(pats+7));		
		 int * txtsLenVec = calloc(pBlockSize,sizeof(int));
	       float * pensVec = malloc(sizeof(float)*2);
	         int * hproVec = calloc(hproVecLen,sizeof(int));
	         int * dproVec = calloc(dproVecLen,sizeof(int));
       unsigned int ** H;

	H = malloc((maxTxtLen+1)*sizeof(unsigned int*));
	H[0] = calloc ((maxTxtLen+1)*(maxPatLen+1), sizeof(unsigned int));
	for ( i = 1 ; i < maxTxtLen + 1 ; ++ i )
     		H[i] = (void*)H[0] + i*(maxPatLen+1)* sizeof(unsigned int);


	if(patsVec==NULL   || txtsVec==NULL      || argsVec==NULL || 
           pensVec == NULL || txtsLenVec == NULL || hproVec==NULL || 
	   dproVec==NULL || H==NULL)
	{
	 	errno = MALLOC;
      		return ( 1 );	
	}

	fill_txtsVec (txts, pBlockSize, t, txtsVec,in->scoring_matrix);
	fill_patsVec (pats, maxPatLen, p, patsVec,in->scoring_matrix);
	fill_argsVec (pats, txts, p, in->max_gap, pBlockSize, maxPatLen, maxTxtLen, argsVec);
	fill_txtsLenVec (txts, t, txtsLenVec);

	pensVec[0] = - in -> gap_open_pen;
	pensVec[1] = - in -> gap_extend_pen;


	if(in->scoring_matrix==0)
		kernel_one_to_one_dna (0, 0, patsVec, txtsVec, argsVec, txtsLenVec, pensVec, hproVec, dproVec, scores, H, out, &start);
	else
		kernel_one_to_one_pro (0, 0, patsVec, txtsVec, argsVec, txtsLenVec, pensVec, hproVec, dproVec, scores, H, out, &start);  

   	if ( out -> min_gap > 0 ) 
        	backtracing ( H, maxPatLen, maxTxtLen, start, out );


   	if ( out -> where == 1 ) 
     		num_mismatch ( t[0], maxTxtLen, p[0], maxPatLen, out );
   	else                   
     		num_mismatch ( p[0], maxPatLen, t[0], maxTxtLen, out );

	free (txtsVec);
	free (patsVec);
	free (argsVec);
	free (txtsLenVec);
	free (pensVec);
	free (hproVec);
	free (dproVec);
        free ( H[0] );
	free (H);
	
	return 0;
}

static cl_platform_id get_gpu_id(int * error)
{
	cl_platform_id gpu_id = NULL;

	cl_uint platforms_total=0;
	
	cl_int err;

	err = clGetPlatformIDs (0, NULL, &platforms_total);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}
	if(platforms_total<=0)
	{
		*error = 1;
		return NULL;
	}
	cl_platform_id * gpu_id_vec = malloc(sizeof(cl_platform_id) * platforms_total);

	err = clGetPlatformIDs (platforms_total, gpu_id_vec, NULL);

   // printf( "platforms: %d\n", platforms_total );
    int i;
    int use_platform = -1;
    
    
    // choose the nvidia platform
    for( i = 0; i < platforms_total; ++i ) {
        char str[256];
        
        clGetPlatformInfo( gpu_id_vec[i], CL_PLATFORM_VENDOR, sizeof(str), str, NULL );
        
     //   printf( "vendor: %s\n", str );
        
        if( strstr( str, "NVIDIA" ) != NULL ) {
            use_platform = i;
            break;
        }
    }
    
	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return NULL;
	}
	
	if( use_platform == -1 ) {
        *error = 1;
        return NULL;
    }
	gpu_id = gpu_id_vec[use_platform];

	free( gpu_id_vec );
	*error = 0;

	return gpu_id;
}

static cl_device_id get_dev_id(cl_platform_id  gpu_id, int * error)
{
	cl_device_id dev_id = NULL;

	cl_uint devices_total;

	cl_int err;
	
	err = clGetDeviceIDs(gpu_id, CL_DEVICE_TYPE_GPU, 0, NULL, &devices_total);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	if(devices_total<=0)
	{
		*error = 1;
		return NULL;
	}

	cl_device_id * dev_id_vec = malloc(sizeof(cl_device_id) * devices_total);

	err = clGetDeviceIDs(gpu_id, CL_DEVICE_TYPE_ALL, devices_total, dev_id_vec, NULL);

    
    
    
	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	dev_id = dev_id_vec[0];

	free( dev_id_vec );
	*error = 0;
	return dev_id;
}

static cl_context create_context(cl_device_id dev_id, int * error)
{


	cl_context context;

	cl_int err;
    
	context = clCreateContext (0,1, &dev_id, NULL,NULL, &err);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}
	
	*error = 0;
	return context;
}

static cl_command_queue create_cmd_queue (cl_device_id dev_id, cl_context context, int * error)
{
	cl_int err;

	cl_command_queue cmd_queue;

	cmd_queue = clCreateCommandQueue(context, dev_id, 0, &err);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	return cmd_queue;
}

static cl_kernel load_kernel (char * name, char * kernel_name, cl_device_id dev_id, cl_context context, int * error)
{
	cl_kernel kernel;

	FILE * fp = fopen(name, "r");

	if(fp==NULL)
	{
		*error = 1;
		return NULL;
	}

	char * source;
	size_t size;

	fseek(fp, 0, SEEK_END);
	size = ftell(fp);
	fseek(fp, 0, SEEK_SET);
		
	source = malloc((size+1)*sizeof(char));
	
	size = fread(source, 1, size, fp);

	fclose(fp);
	
	source[size] = '\0';

	cl_int err;

	cl_program program = clCreateProgramWithSource (context, 1, (const char **) &source, &size, &err);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	clBuildProgram(program, 0, NULL, NULL, NULL, NULL);

	cl_build_status status;

	clGetProgramBuildInfo(program, dev_id, CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &status, NULL);

	if(status!=CL_BUILD_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	kernel = clCreateKernel(program, kernel_name , &err);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	*error = 0;
        free ( source );
	clReleaseProgram(program);
	return kernel;
}

static cl_mem malloc_device (cl_context context, size_t size, int * error)
{
	cl_mem mem = NULL;

	cl_int err;

	mem = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &err);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return NULL;
	}

	*error=0;
	return mem;	
}

static void init_device_mem_int (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, int * mem, size_t size, int * error)
{
	cl_int err;

	err = clEnqueueWriteBuffer(cmd_queue, dev_mem, CL_FALSE, 0, size * sizeof(int), mem, 0, NULL, NULL);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	err = clFinish(cmd_queue);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	*error = 0;
	return;
}

static void init_device_mem_uint (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, unsigned int * mem, size_t size, int * error)
{
	cl_int err;

	err = clEnqueueWriteBuffer(cmd_queue, dev_mem, CL_FALSE, 0, size * sizeof(unsigned int), mem, 0, NULL, NULL);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	err = clFinish(cmd_queue);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	*error = 0;
	return;
}

static void init_device_mem_float (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, float * mem, size_t size, int * error)
{
	cl_int err;

	err = clEnqueueWriteBuffer(cmd_queue, dev_mem, CL_FALSE, 0, size * sizeof(float), mem, 0, NULL, NULL);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	err = clFinish(cmd_queue);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	*error = 0;
	return;
}

static void set_kernel_arguments (cl_kernel kernel, cl_command_queue cmd_queue, cl_mem cl_mem0, cl_mem cl_mem1, cl_mem cl_mem2, cl_mem cl_mem3, cl_mem cl_mem4, cl_mem cl_mem5, cl_mem cl_mem6, cl_mem cl_mem7)
{
	clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &cl_mem0);
	clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &cl_mem1);
	clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &cl_mem2);
	clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &cl_mem3);
	clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &cl_mem4);
	clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &cl_mem5);
	clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &cl_mem6);
	clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &cl_mem7);
	clFinish(cmd_queue);
}

static void read_device_mem_float (cl_command_queue cmd_queue, size_t size, float * mem, cl_mem dev_mem, int * error)
{
	cl_int err;

	err = clEnqueueReadBuffer(cmd_queue, dev_mem, CL_FALSE, 0, size * sizeof(float), mem, 0, NULL, NULL);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	err = clFinish(cmd_queue);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}
	
	*error = 0;
	return;
}

static void update_group_match (float * groupScores, int * groupMatch, float * groupMatchScores, unsigned int patGroupSize, unsigned int txtGroupSize, int pats, int txts, int patGroupIndex, int txtGroupIndex)
{
	int i,j;

	float max_score;
	int position=-1;

	for(i=0;i<patGroupSize;i++)
	{
		if(patGroupIndex*patGroupSize+i<pats)
		{		
			max_score=groupMatchScores[i];
			position=groupMatch[i];
	
			for(j=0;j<txtGroupSize;j++)
			{		
				if( groupScores[i*txtGroupSize+j] > max_score && txtGroupIndex*txtGroupSize + j < txts)
				{
					max_score = groupScores[i*txtGroupSize+j];
					position = txtGroupIndex*txtGroupSize+j;	
				}
			}

			groupMatch[i]=position;
			groupMatchScores[i]=max_score;
		}
	}
}

static void set_invalid(int * groupMatch, int groupSize)
{
	int i;
	for(i=0;i<groupSize;i++)
		groupMatch[i]=-1;
}

static void set_minimum(float * groupMatchScores, int groupSize)
{
	int i;
	for(i=0;i<groupSize;i++)
		groupMatchScores[i]=-DBL_MAX;
}

static unsigned int get_number_of_groups (int elements, int groupSize)
{	
	unsigned int div32 = elements / groupSize;
	unsigned int mod32 = elements % groupSize;
	
	unsigned int groups = mod32!=0?(div32+1):div32;	
	
	return groups;
}

static void set_null (const char ** input, int size)
{
	int i;
	for(i=0;i<size;i++)
		input[i]=NULL;
}

static void initialize_pointers (const char * groupPatterns[], int groupIndex, int groupSize, const char const ** source, int sourceSize)
{
	int i, elements = sourceSize - groupIndex * groupSize;

	if(elements>=groupSize)
	{
		for(i=0;i<groupSize;i++)
			groupPatterns[i]=source[i + groupIndex * groupSize];
	}
	else
	{
		for(i=0;i<elements;i++)
			groupPatterns[i]=source[i + groupIndex * groupSize];

		for(i=elements;i<groupSize;i++)
			groupPatterns[i]=NULL;
	}
}

static unsigned int get_number_of_sequences (const char ** input)
{
	unsigned int 	total = 0;
	const char 	** Tmp;

	for ( Tmp = input; *Tmp; ++ Tmp, ++ total );

	return total;  
}

static unsigned int get_max_length (unsigned int total, const char ** input)
{
	unsigned int	i, curLen, maxLen=0;
   
	for (i=0;i<total;i++)
	{
		curLen = strlen(input[i]);
		if(curLen>maxLen)
			maxLen = curLen;
	} 

	return maxLen;
}

static unsigned int get_min_length (unsigned int total, const char ** input)
{
	unsigned int	i, curLen, minLen=4294967295u;
   
	for (i=0;i<total;i++)
	{
		curLen = strlen(input[i]);
		if(curLen<minLen)
			minLen = curLen;
	} 

	return minLen;
}

static unsigned int get_pblock_size (unsigned int input, unsigned int mult)
{
	unsigned int div32 = input / mult;
	unsigned int mod32 = input % mult;
	
	unsigned int result = mod32!=0?(div32+1)*mult:div32*mult;	
	
	return result;
}

static void fill_txtsVec (unsigned int total, unsigned int blockSize, const char const ** in, unsigned int * outVec, int matrix)
{
	unsigned int i,j,len;
	
	if(matrix==0)
		for(i=0;i<total;i++)
		{		
			len = strlen(in[i]);
	
			for(j=0;j<len;j++)
				outVec[j*blockSize + i] = nuc_char_to_index ( in[i][j]);
		}
	else
		for(i=0;i<total;i++)
		{		
			len = strlen(in[i]);
	
			for(j=0;j<len;j++)
				outVec[j*blockSize + i] = pro_char_to_index ( in[i][j]);
		}
}

static void fill_patsVec (unsigned int total, unsigned int blockSize, const char const ** in, unsigned int * outVec, int matrix)
{
	unsigned int i,j,len,spos;

	if(matrix==0)
		for(i=0;i<total;i++)
		{		
			len = strlen(in[i]);
			spos = i * blockSize;

			for(j=0;j<len;j++)
				outVec[spos + j] = nuc_char_to_index ( in[i][j]);
		}
	else
		for(i=0;i<total;i++)
		{		
			len = strlen(in[i]);
			spos = i * blockSize;

			for(j=0;j<len;j++)
				outVec[spos + j] = pro_char_to_index ( in[i][j]);
		}	
}

static void fill_argsVec (unsigned int totalPats, unsigned int totalTxts, const char const ** in, unsigned int max_gap, unsigned int pBlockSize, unsigned int maxPatLen, unsigned int maxTxtLen, int * outVec)
{
	unsigned int i;

	outVec[0] = totalPats;
	outVec[1] = totalTxts;
	outVec[2] = max_gap;
	outVec[3] = pBlockSize;
	outVec[4] = maxPatLen;
	outVec[5] = pBlockSize * (maxPatLen + 1);
	outVec[6] = pBlockSize * (maxTxtLen + 1);

	for(i=7;i<totalPats+7;i++)		
		outVec[i] = strlen(in[i-7]);
}

static void fill_txtsLenVec (unsigned int totalTxts, const char const ** in, int * outVec)
{
	unsigned int i;

	for(i=0;i<totalTxts;i++)		
		outVec[i] = strlen(in[i]);
}

#endif




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
   double			scr;
   
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

///////////////////////////////////////////////////////// SEQUENTIAL CALLS ////////////////////////////////////////////////////////////////////////
#if 1
   /* ONE_TO_ONE test */
   if ( gapmis_one_to_one ( querys[0], targets[0], &in, out ) )
     {
      /* ONE_TO_ONE result */
       results_one_to_one ( "ONE_TO_ONE", querys[0], querysId[0], targets[0], targetsId[0], &in, &out[0] );
     }

   /* ONE_TO_MANY test */
   if ( gapmis_one_to_many ( querys[0], targets, &in, out ) )
    {
    /* ONE_TO_MANY result */
      results_one_to_many ( "ONE_TO_MANY", querys[0], querysId[0], targets, targetsId, &in, out );
    }
   
   /* MANY_TO_MANY test */
   if ( gapmis_many_to_many ( querys, targets, &in, out ) )
    {
      /* MANY_TO_MANY results */
      results_many_to_many ( "MANY_TO_MANY", querys, querysId, targets, targetsId, &in, out );
    }
   /* ONE_TO_ONE_scr test */
   if ( ! ( gapmis_one_to_one_scr ( querys[0], targets[0], &in, &scr ) ) )
     {
       printf ("\nErrno:%d", errno);
     }

   /* ONE_TO_MANY_opt test */
   if ( ! ( gapmis_one_to_many_opt ( querys[0], targets, &in, out ) ) )
     {
       printf ("\nErrno:%d", errno);
     }
   
   /* MANY_TO_MANY_opt test */
   if ( ! ( gapmis_many_to_many_opt ( querys, targets, &in, out ) ) )
     {
       printf ("\nErrno:%d", errno);
     }
     /*
   for( i = 0; i < cur_alloc_query; ++i ) {
       printf( "res: %d %f %d\n", i, out[i].max_score, out[i].gap_pos );
       
   }
   */
     
#endif
///////////////////////////////////////////////////////// GPU CALLS /////////////////////////////////////////////////////////////////////////
#ifdef _USE_GPU
   /* ONE_TO_MANY_opt test */
   if ( ! ( gapmis_one_to_many_opt_gpu ( querys[0], targets, &in, out ) ) )
     {
       printf ("\nErrno:%d", errno);
     }
   
   /* MANY_TO_MANY_opt test */
   if ( ! ( gapmis_many_to_many_opt_gpu ( querys, targets, &in, out ) ) )
     {
       printf ("\nErrno:%d", errno);
     }
/* */
#endif
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////// GPU CALLS /////////////////////////////////////////////////////////////////////////
#ifdef _USE_SSE
   /* ONE_TO_MANY_opt test */
   if ( ! ( gapmis_one_to_many_opt_sse ( querys[0], targets, &in, out ) ) )
     {
       printf ("\nErrno:%d", errno);
     }
   
   /* MANY_TO_MANY_opt test */
   if ( ! ( gapmis_many_to_many_opt_sse ( querys, targets, &in, out ) ) )
     {
       printf ("\nErrno:%d", errno);
     }
     /*
   for( i = 0; i < cur_alloc_query; ++i ) {
       printf( "res: %d %f %d\n", i, out[i].max_score, out[i].gap_pos );
       
   }
     */
/* */
#endif
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
