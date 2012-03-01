#include <stdio.h>
#include <string.h>
#include "gapmis_vec.h"



/* used only in main function */
static void print ( const char * label, struct gapmis_align * out )
{
    printf ( "%smax_score: %f\nmin_gap: %d\nwhere: %d\ngap_pos: %d\nnum_mis: %d\n", label, out -> max_score, out -> min_gap, out -> where, out -> gap_pos, out -> num_mis );
}

int main ( int argc, char * argv [] )
{
    unsigned int                 i, j;
    struct gapmis_params         in;
    struct gapmis_align          out;
    struct gapmis_align          out_arr[6];
    struct gapmis_align          out_arr2[9];




    
    const char * p1 = "AAACCCTGCTATAGTCAGTGTGAGACGAACGCATAAAGGAAAGATGTTACAGCCCGTCTGACCTTCAGAGGCTTTATCGCGGACGTAAACCCTATACAGGATCCCAACCTGTATTTATCTTCTCATTGGGAGGAACACGTGGGCAGACTA";
    const char * t1 = "AAACCCTGCTATAGTCAGTGTGTGACGAACGCATAAAGGAAAGATGTTACAGCCAGTCTGACGGTCAGAGGCTTTATCGCGGAGGTAAACCCTATATAGGATCCCAGCCTGTATTTATCTTCTCATTGGGAGGAACACGTGGTTACCCCCCCCCCCCCCCCCCC";

    const char * p2 = "AAACCCTGTATATCCGTGTGAGACGAACGCGTAAGGAAAGTGATACAGCCGTCGACCCTCAGAGGCTTTATCGGGACGTAAACCATATACGGATCCCAACTGTATTTATCGTCTCATGGGAGGAACTGTGGGCAGACTA";

    const char * t2 = "AAACCCTGTATATCCGTGTGAGAGCGAACGCGTAAGGAAAGTGATACAGCCGTTCGACCCTCAGAGGCTGGTCTATCGGGACGTAAACCATGAAATCAATATACGGATCCCAACCTGATATTTGATCGGTCTCATGGGAGGAACTGTGGGCAGACTAAAAAAAA";
    const char * p3 = "AAACCCTGTATATCCTATGGACAACGGTATGACATGATCAGCGTCGCCCCAAGGTATTCGAGCGTAAACCATATACGGATCACACCTGTATTGCGTCTCTGGGAGGAACGAGGGCAACTA";

    const char * t3 = "ACAACCCTTAGACCGTGATGTATATCCTCATGGACAACGGGTTAATCGGAAACATGATCAGCGTCGGCCCCATAGGTATTTTCGAGCGTAAACCATATAACCGGATACGCACGACTACACCTGTATTGCGTCTCTGGGGAATCAGAGGAACGAGGGCAACTAAA";
    
    const char * texts[] = { t2, t1, t3, NULL };
    const char * pats[]  = { p1, p2, p3, NULL };

    in . max_gap        = strlen ( t2 ) - 1;
    in . scoring_matrix = 0;
    in . gap_open_pen   = 10;
    in . gap_extend_pen = 2;

    /* ONE_TO_ONE test */
    gapmis_one_to_one ( p1, t2, &in, &out );
    print ( "ONE_TO_ONE\n", &out );
    printf ( "\n\n" );
//    return 0;
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

    return ( 0 );
}


