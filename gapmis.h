#ifndef         GAPMIS_H
#define         GAPMIS_H

#ifdef _USE_GPU
#include "CL/opencl.h"
#endif

#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)
#define NUC_SCORING_MATRIX_SIZE         15		
#define PRO_SCORING_MATRIX_SIZE         24
#define LINE_LNG 			50
#define MAX 				500
#define ALLOC_SIZE                      100

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

/* Computes only the maximum score of the optimal semi-global alignment between pattern p and text t */
unsigned int gapmis_one_to_one_scr ( const char * p, const char * t, const struct gapmis_params * in, double* scr );

/* Computes the optimal semi-global alignment with the maximum score between a pattern p and all texts t */
unsigned int gapmis_one_to_many_opt ( const char * p, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes the optimal semi-global alignment with the maximum score between each pattern p and all texts t */
unsigned int gapmis_many_to_many_opt ( const char const ** p, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes the optimal semi-global alignment between pattern p and text t */
unsigned int gapmis_one_to_one ( const char * p, const char * t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes the optimal semi-global alignment between a pattern p and all texts t */
unsigned int gapmis_one_to_many ( const char * p, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes the optimal semi-global alignment between each pattern p and all texts t */
unsigned int gapmis_many_to_many ( const char const ** p, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out );

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


#ifdef _USE_GPU
/* Computes the optimal semi-global alignment with the maximum score between a pattern p and all texts t */
unsigned int gapmis_one_to_many_opt_gpu ( const char * p, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes the optimal semi-global alignment with the maximum score between each pattern p and all texts t */
unsigned int gapmis_many_to_many_opt_gpu (const char const ** p, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out );

static cl_platform_id * get_gpu_id(int * error);
static cl_device_id * get_dev_id(cl_platform_id * gpu_id, int * error);
static cl_context create_context(cl_device_id * dev_id, int * error);
static cl_command_queue create_cmd_queue (cl_device_id * dev_id, cl_context context, int * error);
static cl_kernel load_kernel (char * name, char * kernel_name, cl_device_id * dev_id, cl_context context, int * error);
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
#endif

#endif
