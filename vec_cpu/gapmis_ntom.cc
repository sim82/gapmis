#include <vector>
#include <string>
#include <iostream>

#include <boost/thread.hpp>
#include "fasta.h"
#include "gapmis_vec.h"

typedef std::vector<uint8_t> seq;
    

struct work {
    
    size_t b_start_;
    size_t b_end_;
    
    std::vector<gapmis_align> out_;
    
    
};

class worker {
public:
    
    worker( const std::vector<const char*> &gmi_a, const std::vector<seq> &seqb, work *w ) 
    : gmi_a_(gmi_a), lenb_(-1), work_(w)
    {
        assert( work_->b_end_ > work_->b_start_ );
        
        const size_t num_b = work_->b_end_ - work_->b_start_;
        std::cout << "worker: " << work_->b_start_ << " - " << work_->b_end_ << "\n";
        work_->out_.resize( (gmi_a_.size() - 1) * num_b );
        
        gmi_b_.reserve( num_b + 1 );
        
        
        for( size_t i = work_->b_start_; i < work_->b_end_; ++i ) {
        
            gmi_b_.push_back( (char *) seqb.at(i).data() );
        
            if( lenb_ == size_t(-1) ) {
                lenb_ = seqb[i].size() - 1;
            } else {
                if( lenb_ != seqb[i].size() - 1 ) {
                        
                    throw std::runtime_error( "lenb != it->size() - 1\n" );
                }
            }
        
        }
        gmi_b_.push_back(0);
        
    }
    
    
    void operator()() {
        gapmis_params parm;
        
        parm.scoring_matrix = 0;
        parm.gap_open_pen = 10;
        parm.gap_extend_pen = 1;
        parm.max_gap = lenb_ - 1;
        
        
        
        //std::vector<gapmis_align> out( (gmi_a_.size() - 1) * numb );
        gapmis_many_to_many( const_cast<const char **>(gmi_a_.data()), gmi_b_.data(), &parm, work_->out_.data() );
        
//         for( size_t i = 0; i < out.size(); ++i ) {
//             std::cout << "out: " << i << " " << out[i].max_score << "\n";
//         }
    }
    
//     const std::vector<gapmis_align> &get_out() const {
//            
//         return out_;
//     }
private:
    const std::vector<const char*> &gmi_a_;//( seqa.size() );
//     const std::vector<seq> &seqb_;//( seqb.size() );
    
    std::vector<const char*> gmi_b_;//( seqb.size() );
    
    
    size_t lenb_;
    
    work *work_;
    
};

int main( int argc, char *argv[] ) {
    assert( argc == 3 );
    
    std::vector<seq> seqa;
    std::vector<seq> seqb;
    
    std::vector<std::string> namea;
    std::vector<std::string> nameb;
    
    {
        std::ifstream is( argv[1] );
        read_fasta( is, namea, seqa );
    }
    {
        std::ifstream is( argv[2] );
        read_fasta( is, nameb, seqb );
    }
    
    std::cout << "align " << seqa.size() << " x " << seqb.size() << "\n";
    
    
    // build char** input for gapmis
    std::vector<const char*> gmi_a;//( seqa.size() );
    std::vector<const char*> gmi_b;//( seqb.size() );
    
    
    for( std::vector< seq >::iterator it = seqa.begin(); it != seqa.end(); ++it ) {
        it->push_back(0);
        gmi_a.push_back( (char *)it->data() );
    }
    gmi_a.push_back(0);
    
    for( std::vector< seq >::iterator it = seqb.begin(); it != seqb.end(); ++it ) {
         it->push_back(0);
    }
    
    const size_t num_threads = 2;
    
    const size_t b_per_thread = seqb.size() / num_threads;
    
    boost::thread_group tg;
    size_t b_end = 0;
    
    std::vector<work> works(num_threads);
    
    for( size_t i = 0; i < num_threads - 1; ++i ) {
        size_t b_start = b_end;
        b_end += b_per_thread;
        works[i].b_start_ = b_start;
        works[i].b_end_ = b_end;
        
        
       
        tg.create_thread( worker(gmi_a, seqb, &works[i] ) );
    }
    
    works.back().b_start_ = b_end;
    works.back().b_end_ = seqb.size();
    
    worker w(gmi_a, seqb, &works.back());
    w();
    
    tg.join_all();
    
    
    for( size_t i = 0; i < works.size(); ++i ) {
     
        std::cout << "worker: " << i << "\n";
        
        const std::vector<gapmis_align> &out = works.at(i).out_;
        
//         for( size_t j = 0; j < out.size(); ++j ) {
//             std::cout << "out: " << i << " " << j << " " << out[j].max_score << "\n";
//         }
        
    }
    
    
//     size_t lenb = -1;
//     
//     for( std::vector< seq >::iterator it = seqb.begin(); it != seqb.end(); ++it ) {
//         it->push_back(0);
//         gmi_b.push_back( (char *)it->data() );
//         
//         if( lenb == size_t(-1) ) {
//             lenb = it->size() - 1;
//         } else {
//             if( lenb != it->size() - 1 ) {
//                     
//                 throw std::runtime_error( "lenb != it->size() - 1\n" );
//             }
//         }
//     }
//     gmi_b.push_back(0);
//     gapmis_params parm;
//     
//     parm.scoring_matrix = 0;
//     parm.gap_open_pen = 10;
//     parm.gap_extend_pen = 1;
//     parm.max_gap = lenb - 1;
//     
//     std::vector<gapmis_align> out( seqa.size() * seqb.size() );
//     gapmis_many_to_many( gmi_a.data(), gmi_b.data(), &parm, out.data() );
//     
//     for( size_t i = 0; i < out.size(); ++i ) {
//         
//        
//         std::cout << "out: " << i << " " << out[i].max_score << "\n";
//     }
}