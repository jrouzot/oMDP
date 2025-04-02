#ifndef DOWNLINK_COUNT
#define DOWNLINK_COUNT

#include <vector>
#include <numeric>
#include <random>
#include <algorithm>
#include "InstanceSingleWindow.hpp"
#include "RoundAndCast.hpp"

using namespace std;
using namespace dataflow;

class DownlinkCount
{
    public:

    DownlinkCount(vector<InstanceSingleWindow *> ins) : ins_(ins), N_(ins[0]->getNumBuffers()), W_(ins.size()) {
        for(size_t b(0); b < N_; ++b)
        {
            capacities_.push_back(round_and_cast_up<double, int64_t>(ins[0]->getBufferCapacity(b)));
        }
    }


    /*
     *  Run downlinkcount heuristic to rank the priority variables according to the fill rates and the memory state.
     *  @mem : memory state
     *  @window_index : the window_index to run the heuristic from
     *  @rank : the rank vector
     *  @ratio : The overflow is reach at capacity * ratio
     * 
     *  returns a priority variable index
     */
    void GetRank(vector<int64_t> &mem, int window_index, double ratio, vector<vector<int>> &rank)
    {
        rank = {};  // Reset rank

        vector<int> indexes = {};
            
        for(int i(0); i < N_; ++i)  // Init indexes
        {
            indexes.push_back(i);
        }

        while(window_index < W_)
        {
            vector<int> r = {};  // Current rank spot

            for(auto i : indexes)
            {
                mem[i] += round_and_cast_up<double, int64_t>(ins_[window_index]->getTotalFill(i));
                if(mem[i] > capacities_[i]*ratio)
                {
                    r.push_back(i);                                                                 // Add index to the current rank
                    indexes.erase(std::remove(indexes.begin(), indexes.end(), i), indexes.end());   // Remove index
                }
            }
            if(r.size() != 0)
                rank.push_back(r);
            ++window_index;
        }
    }


    /*
     *  Run downlinkcount heuristic to find a the most prioritary variable according to the fill rates and a subset of variables.
     *  @mem : The memory state
     *  @ratio : The overflow is reach at capacity * ratio
     *  @window_index : The window_index to run the heuristic from
     *  @buffer_index : The buffer for which we want the priority assignement
     * 
     *  returns a priority variable index
     */
    int GetNextPriority(vector<int64_t> &mem, double ratio, int window_index, int buffer_index)
    {
        vector<int> indexes = {};
            
        for(int i(0); i < N_; ++i)  // Init indexes
        {
            indexes.push_back(i);
        }

        int rank = 0;

        while(window_index < W_)
        {
            bool assign = false; 
            for(auto i : indexes)
            {
                mem[i] += round_and_cast_up<double, int64_t>(ins_[window_index]->getTotalFill(i));
                if(mem[i] > capacities_[i]*ratio)
                {
                    if(i == buffer_index)
                        return rank;
                    assign = true;  
                    indexes.erase(std::remove(indexes.begin(), indexes.end(), i), indexes.end());   // Remove index
                }
            }
            if(assign)
                ++rank;
            ++window_index;
        }

        // Return first index if no overflow
        return indexes[0];
    }


    /*
     *  Run downlinkcount heuristic to find a solution.
     *  @mem : memory state
     *  @ratio : The overflow is reach at capacity * ratio
     *  @window_index : the window_index to run the heuristic from
     *  @solution : The solution matrix
     */
    void GetSolution(vector<int64_t> &mem, double ratio, int window_index, vector<vector<int>> &solution)
    {
        solution = {};              // Init solution

        vector<int> indexes = {};
        vector<int> indexes_tmp = {};

        for(int i(0); i < N_; ++i)  // Init indexes
        {
            indexes.push_back(i);
        }

        int rank = 0;
        bool assign;

        while(window_index < W_)
        {
            bool assign = false;
            indexes_tmp = indexes;
            for(auto i : indexes_tmp)
            {
                mem[i] += ins_[window_index]->getTotalFillInt(i);
                if(mem[i] > capacities_[i]*ratio)
                {
                    if(!assign)
                        solution.push_back({});
                    solution[rank].push_back(i);
                    assign = true;
                    indexes.erase(std::remove(indexes.begin(), indexes.end(), i), indexes.end());   // Remove index
                }
            }
            if(assign)
                ++rank;
            ++window_index;
        }
        if(indexes.size())
            solution.push_back({});

        // Push last indexes (buffers that doesnt overflow)
        for(auto i : indexes)
        {
            solution[rank].push_back(i);
        }
    }

    private:
        vector<InstanceSingleWindow *> ins_;
        vector<int64_t> capacities_;
        int N_;
        int W_;
};


    /*
     *  Extract the priority of buffer b for a given solution
     *  @b : the buffer index of which we want the priority
     *  @solution : The solution matrix
     */
    int GetPriorityFromSolution(int b, vector<vector<int>> &solution)
    {
        for(auto rank(0); rank < solution.size(); ++rank)
        {
            auto it = std::find(solution[rank].begin(), solution[rank].end(), b);

            if (it != solution[rank].end())
                return rank;
        }
        std::cout << "The buffer is not in the solution" << std::endl;
        return -1;
    }


#endif