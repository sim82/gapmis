
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <getopt.h>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>

//#include "functions.h"
#include "gapmis_vec.h"
#include "EDNAFULL.h"
#include "EBLOSUM62.h"
#include "vec_unit.h"
#include "aligned_buffer.h"
#include "cycle.h"

const size_t NUC_SCORING_MATRIX_SIZE = 15;
const size_t PRO_SCORING_MATRIX_SIZE = 24;
const int ERR = 24;      //error number returned if char_to_index returns an invalid index


template<typename score_t, size_t VW>
class aligner {
    typedef vector_unit<score_t, VW> vu;
    typedef typename vu::vec_t vec_t;


public:


    aligner( size_t n, unsigned int matrix, const std::vector<char> &qstates )
    : n_(n),
      matrix_(matrix),
      aprofile_( n * qstates.size() * VW ),
      out_score_(VW),
      out_min_gap_(VW),
      out_where_(VW),
      out_start_(VW),


      qstates_(qstates),
      qs_back_(256, size_t(-1)),
      ncup_(0)
    {
        //create_aprofile( t, n, aprofile_.begin(), matrix, qstates );

        //std::cout << "states: " << qstates_.size() << "\n";

        for( size_t i = 0; i < qstates_.size(); ++i ) {
            qs_back_.at( qstates_[i] ) = i;
        }

    }


    void reset_profile( const char **t ) {
        create_aprofile( t, n_, aprofile_.begin(), matrix_, qstates_ );
    }


    score_t *get( size_t i, size_t j ) {
        const size_t row_size = (n_ + 1) * VW;
        return g_( row_size * j + i * VW );
    }

    score_t *geth( size_t i, size_t j ) {
        const size_t row_size = (n_ + 1) * VW;
        return h_( row_size * j + i * VW );
    }

    uint64_t ncup() const {
        return ncup_;
    }

    void align( const char *p, const size_t m, const size_t max_gap ) {


        const size_t row_size = (n_ + 1) * VW;
        const size_t matrix_size = row_size * (m+1);

        if( g_.size() < matrix_size ) {
            g_.resize( matrix_size, 0.0 );
            h_.resize( matrix_size );
        }

        for( size_t i = 0; i < m + 1 ; i++) {
            const size_t addr = row_size * i;

            vu::store( vu::set1(i), h_(addr));
            vu::store( vu::setzero(), g_(addr));
        }
        for( size_t j = 0; j < n_ + 1 ; j++)
        {
            const size_t addr = j * VW;

            vu::store( vu::set1(j), h_(addr));
            vu::store( vu::setzero(), g_(addr));
        }
        ncup_ += n_ * m * VW;

#if 0
        for( size_t j = 1; j < m + 1; ++j ) {
            size_t p_comp = qs_back_.at( p[j] );
            assert( p_comp != size_t(-1) );

            score_t * __restrict aprof_iter = aprofile_( n_ * VW * p_comp );



            //   for( size_t i = 1; i < n_ + 1; ++i ) {
            for( size_t i = 1; i < j; ++i ) {
                vec_t matching_score = vu::load( aprof_iter + (i-1) * VW );
                const size_t diag_addr = row_size * (j-1) + (i-1) * VW;
                const size_t cur_addr = row_size * j + i * VW;




                const size_t mdiag_addr = row_size * i + i * VW;


                const vec_t g_diag = vu::load( g_(diag_addr ));
                const vec_t g_mdiag = vu::load( g_(mdiag_addr ));
                const vec_t mis = vu::add( g_diag, matching_score );
                const vec_t gap = g_mdiag;

                const vec_t cmp_mask = vu::cmp_lt( mis, gap );
                const vec_t h = vu::bit_and( cmp_mask, vu::set1( j - i ));

                vu::store( h, h_(cur_addr) );
                const vec_t g = vu::max( mis, gap );
                vu::store( g, g_(cur_addr));

                //                    if( t[i] == p[j] )
                //                    {
                //                        G[i][j] = G[i-1][j-1] + matching_score;
                //                        H[i][j] = 0;
                //                    }
                //                    else
                //                    {
                //                        mis = G[i-1][j-1] + matching_score;
                //                        gap = G[i][i];
                //
                //                        if ( gap > mis )
                //                            H[i][j] = j - i;
                //                        else
                //                            H[i][j] = 0;
                //
                //                        G[i][j] = std::max ( mis, gap );
                //                    }
            }

            if( j < n_ + 1 )
            {
                size_t i = j;

                vec_t matching_score = vu::load( aprof_iter + (i-1) * VW );
                const size_t diag_addr = row_size * (j-1) + (i-1) * VW;
                const size_t cur_addr = row_size * j + i * VW;

                vu::store( vu::add( matching_score, vu::load(g_(diag_addr))), g_(cur_addr) );
                vu::store( vu::setzero(), h_(cur_addr));
            }

            //                    G[i][j] = G[i-1][j-1] + matching_score;
            //                    H[i][j] = 0;


            for( size_t i = j + 1; i < n_ + 1; ++ i)
            {
                vec_t matching_score = vu::load( aprof_iter + (i-1) * VW );
                const size_t diag_addr = row_size * (j-1) + (i-1) * VW;
                const size_t cur_addr = row_size * j + i * VW;



                const size_t mdiag_addr = row_size * j + j * VW;


                const vec_t g_diag = vu::load( g_(diag_addr ));
                const vec_t g_mdiag = vu::load( g_(mdiag_addr ));
                const vec_t mis = vu::add( g_diag, matching_score );
                const vec_t gap = g_mdiag;

                const vec_t cmp_mask = vu::cmp_lt( mis, gap );
                const vec_t h = vu::bit_and( cmp_mask, vu::set1( i - j ));

                vu::store( h, h_(cur_addr) );
                const vec_t g = vu::max( mis, gap );
                vu::store( g, g_(cur_addr));

                //                    if( t[i] == p[j] )
                //                    {
                //                        G[i][j] = G[i-1][j-1] + matching_score;
                //                        H[i][j] = 0;
                //                    }
                //                    else
                //                    {
                //                        mis = G[i-1][j-1] + matching_score;
                //                        gap = G[j][j];
                //
                //                        if ( gap > mis )
                //                            H[i][j] = i - j;
                //                        else
                //                            H[i][j] = 0;
                //
                //                        G[i][j] = std::max ( mis, gap );
                //                    }
            }
        }


#else


        for( size_t j = 1; j < m + 1; ++j ) {
            size_t p_comp = qs_back_.at( p[j-1] );

//            std::cout << "p[j]" << int(p[j]) << "\n";
            assert( p_comp != size_t(-1) );

            score_t * __restrict aprof_iter = aprofile_( n_ * VW * p_comp );



            //for( size_t i = 1; i < n_ + 1; ++i ) {
            const size_t i_min = std::max ( 1, int( int(j) - max_gap ));
            const size_t i_max = std::min ( int(n_), int( int(j) + max_gap ));

            /*const size_t i_min = 1;
            const size_t i_max = n_*/;
            for( size_t i = i_min; i <= i_max; ++i ) {
                vec_t matching_score = vu::load( aprof_iter + (i-1) * VW );
                const size_t diag_addr = row_size * (j-1) + (i-1) * VW;
                const size_t cur_addr = row_size * j + i * VW;

                if( i < j ) {


                    const size_t mdiag_addr = row_size * i + i * VW;


                    const vec_t g_diag = vu::load( g_(diag_addr ));
                    const vec_t g_mdiag = vu::load( g_(mdiag_addr ));
                    const vec_t mis = vu::add( g_diag, matching_score );
                    const vec_t gap = g_mdiag;

                    const vec_t cmp_mask = vu::cmp_lt( mis, gap );
                    const vec_t h = vu::bit_and( cmp_mask, vu::set1( j - i ));

                    vu::store( h, h_(cur_addr) );
                    const vec_t g = vu::max( mis, gap );
                    vu::store( g, g_(cur_addr));
                }
                else if ( i > j )
                {

                    const size_t mdiag_addr = row_size * j + j * VW;


                    const vec_t g_diag = vu::load( g_(diag_addr ));
                    const vec_t g_mdiag = vu::load( g_(mdiag_addr ));
                    const vec_t mis = vu::add( g_diag, matching_score );
                    const vec_t gap = g_mdiag;

                    const vec_t cmp_mask = vu::cmp_lt( mis, gap );
                    const vec_t h = vu::bit_and( cmp_mask, vu::set1( i - j ));

                    vu::store( h, h_(cur_addr) );
                    const vec_t g = vu::max( mis, gap );
                    vu::store( g, g_(cur_addr));

                }
                else if (i == j)
                {
                    vu::store( vu::add( matching_score, vu::load(g_(diag_addr))), g_(cur_addr) );
                    vu::store( vu::setzero(), h_(cur_addr));
                }


            }

        }
#endif
    }
    inline vec_t total_scoring( unsigned int gap, vec_t matrix_score, double gap_open_penalty, double gap_extend_penalty )
    {
        vec_t inc;

        if( gap > 0 ) {
            inc = vu::set1((gap-1) * gap_extend_penalty + gap_open_penalty);
        } else {
            inc = vu::setzero();
        }


        return vu::add( matrix_score, inc );
        //return ( matrix_score + ( ( gap > 0 ) ? ( gap - 1 ) * gap_extend_penalty + gap_open_penalty : 0 ) );
    }

    /* Computes the limits of the i-th coordinate for the matrix G in constant time */
    static unsigned int i_limits ( unsigned int n, unsigned int m, unsigned int * up, unsigned int * down, unsigned int max_gap )
    {
        if ( ( int ) m - ( int ) max_gap < 0 )  (* up )    = 0;
        else                                (* up )    = m - max_gap;
        if ( m + max_gap > n )              (* down )  = n;
        else                                (* down )  = m + max_gap;

        return ( 0 );
    }
    /* Computes the limits of the j-th coordinate for matrix G and H in constant time */
    unsigned int j_limits ( unsigned int i, unsigned int m, unsigned int * left, unsigned int * right, unsigned int MAXgap )
     {
       if ( (int) i - (int) MAXgap > 0 )    (* left )   = i - MAXgap;
       else                                 (* left )   = 1;
       if ( i + MAXgap > m )                (* right )  = m;
       else                                 (* right )  = i + MAXgap;
       return ( 0 );
     }


    /*
    Computes the optimal alignment using matrix G in O(2*MAXgap+1) time
    Note:   double gap_open_penalty, double gap_extend_penalty, double gap_open_offset_penalty are arguments given by the user to represent the gap penalty.
     */
    unsigned int opt_solution ( unsigned int m, unsigned int MAXgap,
            double gap_open_penalty,
            double gap_extend_penalty

    )
    {

        // NOTE: in the vectorized version, the rows have the length of the text sequences, so
        // we search for the best score in the last _row_!!!
        // this seems to be the opposite of the sequential version.

        vec_t score = vu::set1(vu::SMALL_VALUE);
        vec_t min_gap = vu::setzero();
        vec_t where = vu::setzero();
        vec_t start = vu::setzero();

        //unsigned int i, j;

        unsigned int up = 0;
        unsigned int down = 0;
        i_limits( n_, m, &up, &down, MAXgap );                   // computes the i coordinates for matrix G for the last column



        const size_t row_size = (n_ + 1) * VW;
        for ( size_t i = up ; i <= down ; i++ )
        {

//            std::cout << "im: " << i << " " << m << "\n";
            //            double temp_score = 0.0;
            if ( i < m )
            {
                if ( m - i <= MAXgap )
                {
                    vec_t g = vu::load( g_(row_size * m + i * VW));

                    vec_t temp_score = total_scoring ( m - i, g, -gap_open_penalty, -gap_extend_penalty );

//                    score = vu::max( score, temp_score );
                    vec_t cr = vu::cmp_lt( score, temp_score );
                    vec_t ncr = vu::bit_invert(cr);

                                        //score = vu::max( score, temp_score );
                    score = vu::bit_or( vu::bit_and( score, ncr ), vu::bit_and( temp_score, cr ));
                    min_gap = vu::bit_or( vu::bit_and( min_gap, ncr ), vu::bit_and( vu::set1(m - i ), cr ));
                    where = vu::bit_or( vu::bit_and( where, ncr ), vu::bit_and( vu::set1(1), cr ));
                    start = vu::bit_or( vu::bit_and( start, ncr ), vu::bit_and( vu::set1(i), cr ));
//                    std::cout << "tmp ";
//                    vu::println( g );
//                    vu::println( temp_score );
//                    vu::println( score );

                    //                    temp_score = total_scoring ( m - i, G[i][m], gap_open_penalty, gap_extend_penalty );
                    //                    if ( temp_score > score )
                    //                    {
                    //                        score = temp_score;
                    //                        ( *MAXscore ) = score;
                    //                        ( *MINgap ) = m - i;
                    //                        ( *where ) = 1;         //where: gap is in the text and start backtracing from the last column
                    //                        ( *start ) = i;         //backtrace from cell G[start,m]
                    //                    }
                }
            }
            else if ( i > m )
            {
                if ( i - m <= MAXgap )
                {
                    vec_t g = vu::load( g_(row_size * m + i * VW));
                    vec_t temp_score = total_scoring ( i - m, g, -gap_open_penalty, -gap_extend_penalty );

                    vec_t cr = vu::cmp_lt( score, temp_score );
                    vec_t ncr = vu::bit_invert(cr);

                    //score = vu::max( score, temp_score );
                    score = vu::bit_or( vu::bit_and( score, ncr ), vu::bit_and( temp_score, cr ));
                    min_gap = vu::bit_or( vu::bit_and( min_gap, ncr ), vu::bit_and( vu::set1(i-m), cr ));
                    where = vu::bit_or( vu::bit_and( where, ncr ), vu::bit_and( vu::set1(2), cr ));
                    start = vu::bit_or( vu::bit_and( start, ncr ), vu::bit_and( vu::set1(i), cr ));

//                    std::cout << "tmp2 ";
//                    vu::println( g );
//                    vu::println( temp_score );
//                    vu::println( score );

                    //temp_score = total_scoring( i - m, G[i][m], gap_open_penalty, gap_extend_penalty );
                    //                    if (  temp_score > score )
                    //                    {
                    //                        score = temp_score;
                    //                        ( *MAXscore ) = score;
                    //                        ( *MINgap ) = i - m;
                    //                        ( *where ) = 2;         //where: gap is in the pattern and start backtracing from last column
                    //                        ( *start ) = i;         //backtrace from cell G[start,m]
                    //                    }
                }
            }
            else if ( i == m )
            {
                vec_t g = vu::load( g_(row_size * m + i * VW));
                vec_t temp_score = total_scoring ( 0, g, -gap_open_penalty, -gap_extend_penalty );

//                score = vu::max( score, temp_score );
                vec_t cr = vu::cmp_lt( score, temp_score );
                vec_t ncr = vu::bit_invert(cr);
                                    //score = vu::max( score, temp_score );


                score = vu::bit_or( vu::bit_and( score, ncr ), vu::bit_and( temp_score, cr ));
                min_gap = vu::bit_or( vu::bit_and( min_gap, ncr ), vu::bit_and( vu::set1(0), cr ));
                where = vu::bit_or( vu::bit_and( where, ncr ), vu::bit_and( vu::set1(0), cr ));
                start = vu::bit_or( vu::bit_and( start, ncr ), vu::bit_and( vu::set1(m), cr ));
//                std::cout << "tmp3 ";
//                vu::println( g );
//                vu::println( temp_score );
//                vu::println( score );

                //temp_score = total_scoring( 0, G[i][m], gap_open_penalty, gap_extend_penalty );
                //                if (  temp_score > score ) // mgap = 0
                //                {
                //                    score = temp_score;
                //                    ( *MAXscore ) = score;
                //                    ( *MINgap ) = 0;
                //                    ( *where ) = 0;         //there is no gap
                //                    ( *start ) = m;         //no need to backtrace
                //                }
            }
        }

//                unsigned int left = 0;
//                unsigned int right = 0;
//                j_limits ( n_, m, &left, &right, MAXgap );       // computes the j coordinates for matrix G for the last row
//
//                for ( j = left ; j < right ; j++ )
//                {
//                    double temp_score = 0;
//                    if ( n_ - j <= MAXgap )
//                    {
//                        vec_t g = vu::load( g_(row_size * n_ + j * VW));
//                        vec_t temp_score = total_scoring ( n_ - j, g, gap_open_penalty, gap_extend_penalty );
//                        score = vu::max( score, temp_score );
//
//        //                temp_score = total_scoring( n - j, G[n][j], gap_open_penalty, gap_extend_penalty );
//        //                if (  temp_score > score )
//        //                {
//        //                    score = temp_score;
//        //                    ( *MAXscore ) = score;
//        //                    ( *MINgap ) = n - j;
//        //                    ( *where ) = 3;         //where: gap is in the pattern and start backtracing from last row
//        //                    ( *start ) = j;         //backtrace from cell G[n,start]
//        //                }
//                    }
//                }

        vu::store( score, out_score_(0) );
        vu::store( min_gap, out_min_gap_(0) );
        vu::store( where, out_where_(0) );
        vu::store( start, out_start_(0) );
        return 1;
    }


    //    aligned_buffer<score_t>::iterator result_begin() {
    //        return result_.begin();
    //    }
    //
    //    aligned_buffer<score_t>::iterator result_end() {
    //        return result_.end();
    //    }

    score_t get_out_score( size_t idx ) {
        assert( idx < VW );
        return out_score_[idx];
    }
    score_t get_out_min_gap( size_t idx ) {
        assert( idx < VW );
        return out_min_gap_[idx];
    }
    score_t get_out_where( size_t idx ) {
        assert( idx < VW );
        return out_where_[idx];
    }

    score_t get_out_start( size_t idx ) {
        assert( idx < VW );
        return out_start_[idx];
    }

    void backtrace( gapmis_align *out, size_t m, size_t num_valid ) {


        for( size_t v = 0; v < num_valid; ++v ) {
            int         i, j;
            out[v].gap_pos = 0;

//             std::cout << "where: " << out_where_[v] << "\n";
//             std::cout << "start: " << out_start_[v] << "\n";
            
            
            if ( out_where_[v] == 1 || out_where_[v] == 2 )
            {
                i = out_start_[v]; j = m;         //we start backtracing from the last column
            }
            else
            {
                i = n_; j = out_start_[v]; //we start backtracing from the last row
            }
            while ( i >= 0 )
            {
                const size_t row_size = (n_ + 1) * VW;
                int h = h_[row_size * j + i * VW + v];

                if ( h == 0 )
                {
                    --i; --j;
                }
                else
                {
                    out[v].gap_pos = ( i > j ) ? j : i;
                    break;
                }
            }
        }
    }

private:
    /* Returns the score for matching character a and b based on EDNAFULL matrix */
    static int nuc_delta ( char a, char b )
    {
        unsigned int index_a = nuc_char_to_index ( a );
        unsigned int index_b = nuc_char_to_index ( b );

        if ( ( index_a < NUC_SCORING_MATRIX_SIZE ) && ( index_b < NUC_SCORING_MATRIX_SIZE ) )
            return ( EDNAFULL_matrix[ index_a ][ index_b ] );
        else //Error
            return ( ERR );
    }

    /* Returns the score for matching character a and b based on EBLOSUM62 matrix */
    static int pro_delta ( char a, char b )
    {
        unsigned int index_a = pro_char_to_index ( a );
        unsigned int index_b = pro_char_to_index ( b );

        if ( ( index_a < PRO_SCORING_MATRIX_SIZE ) && ( index_b < PRO_SCORING_MATRIX_SIZE ) )
            return ( EBLOSUM62_matrix[ index_a ][ index_b ] );
        else //Error
            return ( ERR );
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
            fprintf ( stderr, "Error: unrecognizable character in one of the nucleotide sequences ('%d')!!!\n", a );
            index = ERR;
            throw std::runtime_error( "bailing out" );
            break;
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
            fprintf ( stderr, "Error: unrecognizable character in one of the protein sequences!!!\n" );
            index = ERR; break;
        }
       return ( index );
     }

    template<typename oiter>
    void create_aprofile( const char **t, size_t n, oiter ostart, unsigned int matrix, const std::vector<char> &qstates ) {

        if( matrix ) {
            for( size_t qs = 0; qs != qstates.size(); ++qs ) {
                for( size_t i = 0; i < n; ++i ) {
                    for( size_t j = 0; j < VW; ++j, ++ostart ) {
                        *ostart =  (score_t) pro_delta( t[j][i], qstates[qs] );
                    }
                }
            }
        } else {
            for( size_t qs = 0; qs != qstates.size(); ++qs ) {


                for( size_t i = 0; i < n; ++i ) {
                    for( size_t j = 0; j < VW; ++j, ++ostart ) {
                        *ostart =  (score_t) nuc_delta( t[j][i], qstates[qs] );
                    }
                }
            }
        }



    }


    const size_t n_;
    const int matrix_;
    aligned_buffer<score_t> aprofile_;
    aligned_buffer<score_t> g_;
    aligned_buffer<score_t> h_;

    aligned_buffer<score_t> out_score_;
    aligned_buffer<score_t> out_min_gap_;
    aligned_buffer<score_t> out_where_;
    aligned_buffer<score_t> out_start_;



    std::vector<char> qstates_;
    std::vector<size_t> qs_back_;

    uint64_t ncup_;
};


/* Computes the optimal semi-global alignment between t and p */
unsigned int gapmis_one_to_one ( const char * p, const char * t, const struct gapmis_params * in, struct gapmis_align * out ) {
    const char *px[2] = {p, 0};
    const char *tx[2] = {t, 0};

    return gapmis_many_to_many( px, tx, in, out );
}

/* Computes the optimal semi-global alignment between a set of texts and p */
unsigned int gapmis_one_to_many ( const char * p, const char ** t, const struct gapmis_params * in, struct gapmis_align * out ) {

    const char *px[2] = {p, 0};
    return gapmis_many_to_many( px, t, in, out );
}

/* Computes the optimal semi-global alignment between a set of factors and a set of patterns */
unsigned int gapmis_many_to_many ( const char ** p, const char ** t, const struct gapmis_params * in, struct gapmis_align * out ) {
    size_t         num_t = 0;

    for ( const char **Tmp = t; *Tmp; ++ Tmp, ++num_t );  //Counting the number of texts


    // states_c controls for which sequences characters the reference profile will be generated
    const char *states_c = in->scoring_matrix == 0 ? "ACGT" : "ARNDCQEGHILKMFPSTWYVBZX";

    const size_t n_states = strlen( states_c );
    const std::vector<char> states( states_c, states_c + n_states);

    const char ** t_iter = t;
    const size_t VW = 8;

    size_t block_start = 0;

    std::vector<size_t> p_sizes;
    std::vector<const char*> p_ptrs;
    for( const char ** p_iter = p; *p_iter != 0; ++p_iter ) {
        p_sizes.push_back(strlen(*p_iter));
        p_ptrs.push_back( *p_iter );
    }

    size_t len = strlen(*t_iter);
    
    aligner<short,VW> ali(len, in->scoring_matrix, states );
    while( true ) {
        const char *block[VW];
        std::fill( block, block + VW, (char *)0 );

        size_t num_valid = 0;



        for( size_t i = 0; i < VW; ++i ) {
            if( * t_iter != 0 ) {
                block[i] = *t_iter;
                ++num_valid;
                ++t_iter;

//                std::cout << "len: " << strlen(block[i]) << "\n";
//                if( len == size_t(-1)) {
//                    len = strlen(block[i]);
//                } else {
                
                assert( len == strlen(block[i]) );
                
//                }

            } else {
                block[i] = block[0];
                assert( block[0] != 0 );
            }
        }

        ali.reset_profile( block );

        for( size_t i = 0; i != p_ptrs.size(); ++i ) {
            ali.align( p_ptrs[i], p_sizes[i], in->max_gap );
            ali.opt_solution( p_sizes[i], in->max_gap, in->gap_open_pen, in->gap_extend_pen );


//            for( size_t x = 0; x < len+1; ++x ) {
//                for( size_t y = 0; y < p_sizes[i]+1; ++y ) {
//                    std::cout << ali.get( x, y )[0] << "\t";
//                }
//                std::cout << "\n";
//            }

            gapmis_align ali_out[VW];
            for( size_t j = 0; j < num_valid;++j ) {
                gapmis_align &x = ali_out[j];
                memset( &x, 0, sizeof( gapmis_align ));

                x.max_score = ali.get_out_score( j );
                x.min_gap = ali.get_out_min_gap( j );
                x.where = ali.get_out_where( j );

                //out[i * num_t + block_start + j] = x;

            }
            if(! false )
            {
                ali.backtrace( ali_out, p_sizes[i], num_valid );
            }
            for( size_t j = 0; j < num_valid;++j ) {
                out[i * num_t + block_start + j] = ali_out[j];
            }


        }


        block_start += VW;
        if( *t_iter == 0 ) {
            break;
        }

    }


    return 1;
}


//
//void gv_init_aligner( gv_aligner_t *c_ali, int n, int matrix, char *states, int n_states ) {
//    aligner<short,8> *ali = new aligner<short,8>( n, matrix, std::vector<char>( states, states + n_states) );
//
//    c_ali->cxx_aligner = ali;
//    c_ali->vec_width = 8;
//}
//void gv_delete_aligner( gv_aligner_t *c_ali ) {
//    assert( c_ali != 0 );
//    assert( c_ali->vec_width == 8 );
//    assert( c_ali->cxx_aligner != 0 );
//
//    aligner<short,8> *ali = reinterpret_cast<aligner<short,8> *> (c_ali->cxx_aligner);
//    delete ali;
//
//    memset( c_ali, 0, sizeof( gv_aligner_t ));
//
//}
//
//
//void gv_reset_text( gv_aligner_t *c_ali, char **t ) {
//    assert( c_ali != 0 );
//    assert( c_ali->vec_width == 8 );
//    assert( c_ali->cxx_aligner != 0 );
//
//    aligner<short,8> *ali = reinterpret_cast<aligner<short,8> *> (c_ali->cxx_aligner);
//
//    ali->reset_profile(t);
//}
//
//
//
//void gv_align( gv_aligner_t *c_ali, char * qs, int qs_len, int max_gap, double gap_open, double gap_extend, gv_result_t *res ) {
//    assert( c_ali != 0 );
//    assert( c_ali->vec_width == 8 );
//    assert( c_ali->cxx_aligner != 0 );
//
//    aligner<short,8> *ali = reinterpret_cast<aligner<short,8> *> (c_ali->cxx_aligner);
//
//    ali->align( qs, qs_len, max_gap );
//    ali->opt_solution( qs_len, max_gap, gap_open, gap_extend );
//
//    for( size_t i = 0; i < 8; ++i ) {
//        res[i].max_score = ali->get_max_score(i);
//    }
//
//
//}
