#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <errno.h>
#include "gapmis.h"
#include "errors.h"
#include "EDNAFULL.h"
#include "EBLOSUM62.h"

		

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
	
   for ( ; *Tmp; ++Tmp, ++ i )	//computing the maximum score
    {
      if ( ! ( gapmis_one_to_one_scr ( p, *Tmp, in, &scr ) ) )
        return ( 0 );
      if ( scr > tmp_scr )
       {
         max_t = i;
         tmp_scr = scr;
       } 
    } 
   
   if ( ! ( gapmis_one_to_one ( p, t[ max_t ], in, out ) ) ) //computing the rest of the alignment with the maximum score
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
      return ( 1 );
    } 
   
   for ( i = 1; i < n + 1; ++ i )
     G[i] = ( void * ) G[0] + i * ( m + 1 ) * sizeof ( int );
     
   /* dynamic programming algorithm */
   if ( ! ( dp_algorithm_scr( G, t, n, p, m, in ) ) )
    {
      //Error: dp_algorithm_scr() failed due to bad character!!!
      return ( 0 );	
    }
   
   /* computes the optimal alignment based on the matrix score and the gap function */
   opt_solution_scr ( G, n, m, in, scr );
   
   free ( G[0] );
   free ( G );	
   return ( 1 );
   
 } 

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
   unsigned int         start;

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
      return ( 1 );
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
   
   /* computes the optimal alignment based on the matrix score and the gap function */
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
	   gap = G[i][i];
	   valM = j - i;

	   if( i > j )	
	     {
	       gap = G[j][j];
               valM = i - j;
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
	   gap = G[i][i];
	   if( i > j )		gap = G[j][j];
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
   else //Error
     return ( BADCHAR );
 }

/* Returns the score for matching character a and b based on EBLOSUM62 matrix */
static int pro_delta ( char a, char b )
 {
   unsigned int index_a = pro_char_to_index ( a );
   unsigned int index_b = pro_char_to_index ( b );

   if ( ( index_a < PRO_SCORING_MATRIX_SIZE ) && ( index_b < PRO_SCORING_MATRIX_SIZE ) )
     return ( EBLOSUM62_matrix[ index_a ][ index_b ] );
   else //Error
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
static unsigned int opt_solution ( int ** G, 
                                   unsigned int n, 
                                   unsigned int m, 
                                   const struct gapmis_params * in,
                                   struct gapmis_align * out,
                                   unsigned int * start 
)
 {
   unsigned int         i;
   double               score = -DBL_MAX;
   unsigned int         up    = 0;
   unsigned int         down  = 0;
   
   i_limits ( n, m, &up, &down, in -> max_gap );			// computes the i coordinates for matrix G for the last column

   for ( i = up ; i <= down ; i++ )
    {
      double temp_score = 0.0;
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
   
   i_limits ( n, m, &up, &down, in -> max_gap );			

   for ( i = up ; i <= down ; i++ )
    {
      double temp_score = 0.0;
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
Note: double matrix_score is only the value of G[i][m], i.e. the score of an alignment WITHOUT the gap penalty
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

/* used only in main function */
static void print ( const char * label, struct gapmis_align * out )
 {
   printf ( "%smax_score: %f\nmin_gap: %d\nwhere: %d\ngap_pos: %d\nnum_mis: %d\n", label, out -> max_score, out -> min_gap, out -> where, out -> gap_pos, out -> num_mis );
 }

int main ( int argc, char * argv [] )
 {
   unsigned int                 i, j;
   double		        scr;
   struct gapmis_params         in;
   struct gapmis_align          out;
   struct gapmis_align          out_arr[3];
   struct gapmis_align          out_arr2[9];

   const char * p1 = "AAACCCTGCTATAGTCAGTGTGAGACGAACGCATAAAGGAAAGATGTTACAGCCCGTCTGACCTTCAGAGGCTTTATCGCGGACGTAAACCCTATACAGGATCCCAACCTGTATTTATCTTCTCATTGGGAGGAACACGTGGGCAGACTA";
   const char * t1 = "AAACCCTGCTATAGTCAGTGTGTGACGAACGCATAAAGGAAAGATGTTACAGCCAGTCTGACGGTCAGAGGCTTTATCGCGGAGGTAAACCCTATATAGGATCCCAGCCTGTATTTATCTTCTCATTGGGAGGAACACGTGGTTACCCCCCCCCCCCCCCCCCC";

   const char * p2 = "AAACCCTGTATATCCGTGTGAGACGAACGCGTAAGGAAAGTGATACAGCCGTCGACCCTCAGAGGCTTTATCGGGACGTAAACCATATACGGATCCCAACTGTATTTATCGTCTCATGGGAGGAACTGTGGGCAGACTA";

   const char * t2 = "AAACCCTGTATATCCGTGTGAGAGCGAACGCGTAAGGAAAGTGATACAGCCGTTCGACCCTCAGAGGCTGGTCTATCGGGACGTAAACCATGAAATCAATATACGGATCCCAACCTGATATTTGATCGGTCTCATGGGAGGAACTGTGGGCAGACTA";
   const char * p3 = "AAACCCTGTATATCCTATGGACAACGGTATGACATGATCAGCGTCGCCCCAAGGTATTCGAGCGTAAACCATATACGGATCACACCTGTATTGCGTCTCTGGGAGGAACGAGGGCAACTA";

   const char * t3 = "ACAACCCTTAGACCGTGATGTATATCCTCATGGACAACGGGTTAATCGGAAACATGATCAGCGTCGGCCCCATAGGTATTTTCGAGCGTAAACCATATAACCGGATACGCACGACTACACCTGTATTGCGTCTCTGGGGAATCAGAGGAACGAGGGCAACTA";

   const char * texts[] = { t2, t1, t3, NULL };
   const char * pats[]  = { p1, p2, p3, NULL };

   in . max_gap        = strlen ( t2 ) - 1;
   in . scoring_matrix = 0;
   in . gap_open_pen   = 10;
   in . gap_extend_pen = 0.5;

   /* ONE_TO_ONE test */
   if ( ! ( gapmis_one_to_one ( p1, t2, &in, &out ) ) )
     {
       printf ( "%d\n", errno );
       return 1;
     } 
   print ( "ONE_TO_ONE\n", &out );
   printf ( "\n\n" );

   /* ONE_TO_MANY test */
   gapmis_one_to_many ( p1, texts, &in, out_arr );
   for ( i = 0; i < 3; ++ i )
     print ( "ONE_TO_MANY\n", &out_arr[i] );
   printf ( "\n\n" );
   
   /* MANY_TO_MANY test */
   gapmis_many_to_many ( pats, texts, &in, out_arr2 );
   for ( i = 0; i < 3; ++ i )
    for ( j = 0; j < 3; ++ j )
     {
       printf ( "%d ", i * 3 + j );
       print ( "MANY_TO_MANY\n", &out_arr2[i * 3 + j] );
     }
   printf ( "\n\n" );

   /* ONE_TO_ONE_SCR test */
   gapmis_one_to_one_scr ( p1, t2, &in, &scr );
   printf ( "ONE_TO_ONE_SCR\n%lf", scr );
   printf ( "\n\n" );

   /* ONE_TO_MANY_OPT test */
   gapmis_one_to_many_opt ( p1, texts, &in, &out );
   print ( "ONE_TO_MANY_OPT\n", &out );
   printf ( "\n\n" );

   /* MANY_TO_MANY_OPT test */
   gapmis_many_to_many_opt ( pats, texts, &in, out_arr );
   for ( i = 0; i < 3; ++ i )
     print ( "MANY_TO_MANY_OPT\n", &out_arr[i] );
   printf ( "\n\n" );

   return ( 0 );
 }
