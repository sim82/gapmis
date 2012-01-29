#ifndef         GAPMIS_H
#define         GAPMIS_H

#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)

struct gapmis_params
 {
   unsigned int         max_gap;
   unsigned int         scoring_matrix;
   double               gap_open_pen;
   double               gap_extend_pen;
 };

struct gapmis_align
 {
   double               max_score;
   unsigned int         min_gap;
   unsigned int         where;
   unsigned int         gap_pos;
   unsigned int         num_mis;
 };

/* Computes the optimal semi-global alignment between t and p */
unsigned int gapmis_one_to_one ( const char * p, const char * t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes the optimal semi-global alignment between a set of texts and p */
unsigned int gapmis_one_to_many ( const char * p, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes the optimal semi-global alignment between a set of factors and a set of patterns */
unsigned int gapmis_many_to_many ( const char const ** p, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out );

static unsigned int dp_algorithm ( int ** G, unsigned int** H, const char* t, unsigned int n, const char * p, unsigned int m, const struct gapmis_params * in );
static int nuc_delta ( char a, char b );
static int pro_delta ( char a, char b );
static unsigned int nuc_char_to_index ( char a );
static unsigned int pro_char_to_index ( char a );
static unsigned int i_limits( unsigned int n, unsigned int m, unsigned int * up, unsigned int * down, unsigned int max_gap );
static unsigned int opt_solution ( int ** G, unsigned int ** H, unsigned int n, unsigned int m, const struct gapmis_params * in, struct gapmis_align * out, unsigned int * start );
static unsigned int backtracing ( unsigned int ** H, unsigned int m, unsigned int n, unsigned int start, struct gapmis_align * out );
static double total_scoring( unsigned int gap, double matrix_score, double gap_open_penalty, double gap_extend_penalty );
static unsigned int num_mismatch ( const char * seqa, unsigned int seqa_len, const char * seqb, unsigned int seqb_len, struct gapmis_align * out );

#endif
