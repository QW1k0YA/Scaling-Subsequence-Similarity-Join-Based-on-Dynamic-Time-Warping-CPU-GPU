
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
using namespace std;
#define POSITIVE_INFINITY 1E20

inline std::size_t cap_start_index_to_window(std::size_t index, std::size_t window){
    if(index>window){ return index-window; } else { return 0; }
}

inline std::size_t cap_stop_index_to_window_or_end(std::size_t index, std::size_t window, std::size_t end){
    
    if(window<end && index+1<end-window){ return index + window + 1; } else { return end; }
}

double MON_dtw(
        const vector<double>& lines,
        const vector<double>& cols,
        const vector<double>& cb,
        int l,
        int w,
        double bsf
)
{
    
    double UB= bsf - cb[w+1];

    const size_t nbcols = l;
    const size_t nblines = l;

    if (w > nblines) { w = nblines; }

    std::vector<double> buffers_v((1+nbcols) * 2, POSITIVE_INFINITY);
    double *buffers = buffers_v.data();
    size_t c{1}, p{nbcols+2};                 

    buffers[c-1] = 0;
    size_t next_start{0};
    size_t pruning_point{0};

    for(size_t i=0; i<nblines; ++i) {
        
        std::swap(c, p);
        if(i+w < nbcols){ UB = bsf - cb[i+w+1]; }
        const double li = lines[i];
        const std::size_t jStop = cap_stop_index_to_window_or_end(i, w, nbcols);
        
        const std::size_t jStart = std::max(cap_start_index_to_window(i, w), next_start);
        
        std::size_t next_pruning_point = jStart; 
        std::size_t j = jStart;
        next_start = jStart;
        
        buffers[c+j-1] = POSITIVE_INFINITY;
        double cost = POSITIVE_INFINITY;
        
        for(; j==next_start && j < pruning_point; ++j) {
            const auto d = DIST(li, cols[j]);
            cost = std::min(buffers[p + j - 1], buffers[p + j]) + d;
            buffers[c + j] = cost;
            if(cost<=UB){ next_pruning_point = j + 1;} else { ++next_start; }
        }
        
        for(; j < pruning_point; ++j) {
            const auto d =  DIST(li, cols[j]);
            cost = std::min(cost, std::min(buffers[p + j - 1], buffers[p + j])) + d;
            buffers[c + j] = cost;
            if(cost<=UB){ next_pruning_point = j + 1;}
        }

        if(j<jStop){
            const auto d = DIST(li, cols[j]);
            if(j==next_start){ 
                cost = buffers[p + j - 1] + d;
                buffers[c + j] = cost;
                if(cost<=UB){ next_pruning_point = j + 1;}
                else
                {return POSITIVE_INFINITY; }
            } else { 
                cost = std::min(cost, buffers[p + j - 1]) + d;
                buffers[c + j] = cost;
                if(cost<=UB){ next_pruning_point = j + 1;}
            }
            ++j;
        } else if(j==next_start)
        { return POSITIVE_INFINITY; }

        for(;j==next_pruning_point && j<jStop;++j){
            const auto d =  DIST(li, cols[j]);
            cost = cost + d;
            buffers[c + j] = cost;
            if(cost<=UB){ ++next_pruning_point; }
        }
        
        pruning_point=next_pruning_point;

    }

    if(pruning_point != nbcols)
    { return POSITIVE_INFINITY; }
    else {
        return buffers[c+nbcols-1];
    }
}
