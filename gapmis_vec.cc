/*
 * Copyright (C) 2012 Simon A. Berger
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 */

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
#include <limits>
#include <cerrno>
//#include "functions.h"

#ifndef _USE_SSE
#define _USE_SSE
#endif

#define _HACK_CLEAN_NAMESPACE
#include "gapmis.h"
#include "EDNAFULL.h"
#include "EBLOSUM62.h"
#include "vec_utils/vec_unit.h"
#include "vec_utils/aligned_buffer.h"
#include "vec_utils/cycle.h"
#include "vec_utils/thread.h"

#include "errors.h"

const size_t NUC_SCORING_MATRIX_SIZE = 15;
const size_t PRO_SCORING_MATRIX_SIZE = 24;
const int ERR = 24;      //error number returned if char_to_index returns an invalid index


class bad_char_error : public std::exception {
public:
    int c_;
    
    bad_char_error( int c ) : c_(c) {}
    virtual ~bad_char_error() throw() {}
    
    
};

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



        for( size_t j = 1; j < m + 1; ++j ) {
            char pc = p[j-1];
            
            if( pc < 0 || size_t(pc) > qs_back_.size() ) {
                throw bad_char_error(pc);
            }
            size_t p_comp = qs_back_.at( p[j-1] );

            if( p_comp == size_t(-1) ) {
                throw bad_char_error(pc);
            }
            
//            std::cout << "p[j]" << int(p[j]) << "\n";
            //assert( p_comp != size_t(-1) );

            score_t * __restrict aprof_iter = aprofile_( n_ * VW * p_comp );



            //for( size_t i = 1; i < n_ + 1; ++i ) {
            const size_t i_min = std::max ( 1, int( int(j) - max_gap ));
            const size_t i_max = std::min ( int(n_), int( int(j) + max_gap ));

            /*const size_t i_min = 1;
            const size_t i_max = n_*/;
            _mm_prefetch( g_(row_size * j-1), _MM_HINT_T0);
#if 0
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
#else
            //for( size_t i = i_min; i <= i_max; ++i ) {
            size_t i = i_min;
//             size_t diag_addr = row_size * (j-1) + (i-1) * VW;
//             size_t cur_addr = row_size * j + i * VW;
            
            for( ; i < j; ++i /*, diag_addr += VW, cur_addr += VW */) {
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
            }
            for( ;i == j; ++i /*, diag_addr += VW, cur_addr += VW */) {
                vec_t matching_score = vu::load( aprof_iter + (i-1) * VW );
                const size_t diag_addr = row_size * (j-1) + (i-1) * VW;
                const size_t cur_addr = row_size * j + i * VW;
                vu::store( vu::add( matching_score, vu::load(g_(diag_addr))), g_(cur_addr) );
                vu::store( vu::setzero(), h_(cur_addr));
                
            }
            for( ; i <= i_max; ++i/*, diag_addr += VW, cur_addr += VW */) {
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
                
            }
            

            
#endif

        }
    }
#if 1
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
#endif
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
                //int h = h_[row_size * j + i * VW + v];
                int h = h_.at(row_size * j + i * VW + v);
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

        if ( ( index_a < NUC_SCORING_MATRIX_SIZE ) && ( index_b < NUC_SCORING_MATRIX_SIZE ) ) {
            return ( EDNAFULL_matrix[ index_a ][ index_b ] );
        } else {
            assert(0 && "nuc_char_to_index did not throw exception as it was supposed to do" );
        }
    }

    /* Returns the score for matching character a and b based on EBLOSUM62 matrix */
    static int pro_delta ( char a, char b )
    {
        unsigned int index_a = pro_char_to_index ( a );
        unsigned int index_b = pro_char_to_index ( b );

        if ( ( index_a < PRO_SCORING_MATRIX_SIZE ) && ( index_b < PRO_SCORING_MATRIX_SIZE ) ) {
            return ( EBLOSUM62_matrix[ index_a ][ index_b ] );
        } else {
            assert(0 && "pro_char_to_index did not throw exception as it was supposed to do" );
        }
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
            
            index = ERR;
            throw bad_char_error(a);
            //break;
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
            throw bad_char_error(a);
            //index = ERR; break;
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




/* Computes the optimal semi-global alignment between a set of texts and p */
unsigned int gapmis_one_to_many_opt_sse ( const char * p, const char ** t, const struct gapmis_params * in, struct gapmis_align * out ) {

    const char *px[2] = {p, 0};
    return gapmis_many_to_many_opt_sse( px, t, in, out );
}

static size_t s_num_threads = 1;

void gapmis_sse_hint_num_threads( size_t num ) {
    s_num_threads = num;
}

template<size_t VW>
struct block {
    const char *seqs[VW];
    size_t block_start;
    size_t num_valid;
};

// private state per worker thread
struct worker_private {
    // private state
    std::vector<double> max_scores_; //( p_ptrs.size(), -std::numeric_limits<double>::infinity());
    std::vector<size_t> max_text_idx_; //( p_ptrs.size(), size_t(-1));
    
    int errno_;
};


// worker thread main object: keeps references to shared and private state.
template<typename score_t, size_t VW>
class worker {
public:
    worker( const std::vector<block<VW> > &blocks, size_t rank, size_t stride, size_t len, const struct gapmis_params * in, const std::vector<size_t> &p_sizes, const std::vector<const char*> &p_ptrs, const std::vector<char> &states, worker_private *priv ) 
    : blocks_(blocks),
      rank_(rank), 
      stride_(stride), 
      len_(len), 
      in_(in), 
      p_sizes_(p_sizes), 
      p_ptrs_(p_ptrs), 
      states_(states),
      priv_(priv)
    {
        priv_->max_scores_.resize( p_ptrs_.size(), -std::numeric_limits<double>::infinity());
        priv_->max_text_idx_.resize( p_ptrs_.size(), size_t(-1));
        priv_->errno_ = 0;
    }
    
    void operator()() {
        priv_->errno_ = UNHANDLED_INTERNAL;
        try {
        
            aligner<score_t,VW> ali(len_, in_->scoring_matrix, states_ );
            
            typename std::vector<block<VW> >::const_iterator it = blocks_.begin() + rank_;
            for( ; it < blocks_.end(); it += stride_ ) 
            {
                block<VW> block = *it;
//                 std::cout << "block: " << block.block_start << " " << rank_ << "\n";
                ali.reset_profile( block.seqs );
                
                for( size_t i = 0; i != p_ptrs_.size(); ++i ) {
                    ali.align( p_ptrs_[i], p_sizes_[i], in_->max_gap );
                    ali.opt_solution( p_sizes_[i], in_->max_gap, in_->gap_open_pen, in_->gap_extend_pen );
                    
                    //                 std::cout << "align:\n";
                    
                    for( size_t j = 0; j < block.num_valid;++j ) {
                        double score = ali.get_out_score(j);
                        
                        if( score > priv_->max_scores_[i] ) {
                            priv_->max_scores_[i] = score;
                            priv_->max_text_idx_[i] = block.block_start + j;
                        }
                        //out[i * num_t + block_start + j] = x;
                    }
                    
                }
            }
            
            priv_->errno_ = 0; // only if this point is reached, there was no error.
            
        } catch( std::bad_alloc x ) { /* in sutter we trust... LALALA my code is now exception safe LALALA */
            priv_->errno_ = MALLOC;
        } catch( bad_char_error x ) {
            priv_->errno_ = BADCHAR;
        } catch( std::logic_error x ) {
            std::cerr << "BUG: unhandled std::logic_error: " << x.what() << std::endl;
            priv_->errno_ = UNHANDLED_INTERNAL;
        } catch( std::exception x ) {
            std::cerr << "BUG: unhandled std::exception: " << x.what() << std::endl;
            priv_->errno_ = UNHANDLED_INTERNAL;
        } catch(...) {
            std::cerr << "BUG: unhandled exception of unknown type (sorry, for the useless message...)" << std::endl;
            priv_->errno_ = UNHANDLED_INTERNAL;
        }
        
        
    }
    
    
private:
    // shared read-only state
    const std::vector<block<VW> > &blocks_;
    const size_t rank_;
    const size_t stride_;
    const size_t len_;
    const struct gapmis_params * in_;
    const std::vector<size_t> &p_sizes_;
    const std::vector<const char*> &p_ptrs_;
    const std::vector<char> &states_;
    
    worker_private *priv_;
};

/* Computes the optimal semi-global alignment between a set of factors and a set of patterns */
unsigned int gapmis_many_to_many_opt_sse ( const char ** p, const char ** t, const struct gapmis_params * in, struct gapmis_align * out ) {
    try {
        
        size_t         num_t = 0;
        
        for ( const char **Tmp = t; *Tmp; ++ Tmp, ++num_t );  //Counting the number of texts


    // states_c controls for which sequences characters the reference profile will be generated
        const char *states_c = in->scoring_matrix == 0 ? "ACGT" : "ARNDCQEGHILKMFPSTWYVBZX";
        
        const size_t n_states = strlen( states_c );
        const std::vector<char> states( states_c, states_c + n_states);
        
        const char ** t_iter = t;
        
        //     typedef short score_t;
        const size_t vec_width = 16; // width of the vector unit in bytes. Either 16 (16bytes = 128bit = sse) or 32 (=256bit = AVX)
        typedef float score_t;
        const size_t VW = vec_width / sizeof(score_t);
        
        // copy pattern pointes to vector and pre-calculate the pattern lengths
        std::vector<size_t> p_sizes;
        std::vector<const char*> p_ptrs;
        for( const char ** p_iter = p; *p_iter != 0; ++p_iter ) {
            p_sizes.push_back(strlen(*p_iter));
            p_ptrs.push_back( *p_iter );
        }
        
        
        size_t len = strlen(*t_iter);
        
        
        // generate the list of blocks
        size_t block_start = 0;
        std::vector<block<VW> > blocks;
        
        while( true ) {
            //const char *block[VW];
            
            block<VW> block;
            std::fill( block.seqs, block.seqs + VW, (char *)0 );
            
            block.num_valid = 0;
            
            for( size_t i = 0; i < VW; ++i ) {
                if( * t_iter != 0 ) {
                    block.seqs[i] = *t_iter;
                    ++block.num_valid;
                    ++t_iter;

                    
                    if( len != strlen(block.seqs[i]) ) {
                        errno = TEXTLEN;
                        return 0;
                    }
                    
                } else {
                    block.seqs[i] = block.seqs[0];
                    assert( block.seqs[0] != 0 );
                }
            }
            
            block.block_start = block_start;
            blocks.push_back( block );
            
            
            block_start += VW;
            if( *t_iter == 0 ) {
                break;
            }
        }
        
        // the s_num_threads > 48 check is only meant as temporary sanity check.
        if( s_num_threads == 0 || s_num_threads > 48 ) {
            errno = THREADCOUNT;
            return 0;
        }
        const size_t num_threads = s_num_threads;
        //s_num_threads = 2;
        
        
        // create worker threads.
        // all data supplied to the workers, except for the wpriv[] ptrs are shared/read-only.
        std::vector<worker_private> wpriv( num_threads );
        thread_group tg;
        for( size_t i = 0; i < num_threads; ++i ) {
            tg.create_thread(worker<score_t, VW>(blocks, i, num_threads, len, in, p_sizes, p_ptrs, states, &wpriv[i]));
        }
        
        // wait for threads to finish
        tg.join_all();

        // check thread error states and propagate, if necessary
        for( size_t i = 0; i < num_threads; ++i ) {
            if( wpriv[i].errno_ != 0 ) {
                errno = wpriv[i].errno_;
                std::cerr << "errno set by thread " << i << ": " << wpriv[i].errno_ << "\n";
                return 0;
            }
        }
        
        assert( !wpriv.empty() );
        // get results fromt first thread as reference
        std::vector<double> max_scores = wpriv.front().max_scores_; //( p_ptrs.size(), -std::numeric_limits<double>::infinity());
        std::vector<size_t> max_text_idx = wpriv.front().max_text_idx_; //( p_ptrs.size(), size_t(-1));
        
        // merge result scores from threads
        for( size_t j = 1; j < num_threads; ++j ) {

            for( size_t i = 0; i < max_scores.size(); ++i ) {
                
                // if two threads return an equal scores (=float equality) for a pattern, choose the pattern with the lower index.
                // this should give equivalent results to the single threaded version.
                if( wpriv[j].max_scores_[i] > max_scores[i] || (wpriv[j].max_scores_[i] == max_scores[i] && wpriv[j].max_text_idx_[i] < max_text_idx[i] )) {
                    max_scores[i] = wpriv[j].max_scores_[i];
                    max_text_idx[i] = wpriv[j].max_text_idx_[i];
                }
            }        
        }
        
        //computing the rest details of the alignment with the maximum score
        for( size_t i = 0; i != p_ptrs.size(); ++i ) {
            size_t text_idx = max_text_idx[i];
            double score = max_scores[i];
            
            assert( text_idx < num_t );
            
            gapmis_one_to_one( p[i], t[text_idx], in, &out[i] );
            
            assert( score == out[i].max_score );
            
        }
    } catch( std::bad_alloc x ) { /* in sutter we trust... LALALA my code is now exception safe LALALA */
        errno = MALLOC;
        return 0;
    } catch( bad_char_error x ) {
        errno = BADCHAR;
        return 0;
    }
    
    errno = 0;
    return 1;
}


