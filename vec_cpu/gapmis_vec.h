#ifndef gapmis_vec_h__
#define gapmis_vec_h__


#ifdef __cplusplus
extern "C" {
#endif

//
//typedef struct gv_aligner {
//    void *cxx_aligner;
//    int vec_width;
//} gv_aligner_t;
//
//
//typedef struct gv_result {
//    double max_score;
//    unsigned int min_gap;
//    unsigned int where;
//    unsigned int start;
//} gv_result_t;
//
//void gv_init_aligner( gv_aligner_t *c_ali, int n, int matrix, char *states, int n_states );
//void gv_delete_aligner( gv_aligner_t *c_ali );
//
//void gv_reset_text( gv_aligner_t *c_ali, char **t );
//void gv_align( gv_aligner_t *c_ali, char * qs, int qs_len, int max_gap, double gap_open, double gap_extend, gv_result_t *res );


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
unsigned int gapmis_one_to_many ( const char * p, const char ** t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes the optimal semi-global alignment between a set of factors and a set of patterns */
unsigned int gapmis_many_to_many ( const char ** p, const char ** t, const struct gapmis_params * in, struct gapmis_align * out );


#ifdef __cplusplus
}
#endif


#endif



