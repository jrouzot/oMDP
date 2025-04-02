#ifndef DOWNLINK_COUNT_DECISION_BUILDER
#define DOWNLINK_COUNT_DECISION_BUILDER

#include <vector>
#include <random>

#include "ortools/constraint_solver/constraint_solver.h"
#include "ortools/constraint_solver/constraint_solveri.h"
#include "../utils/header/DownlinkCount.hpp"
#include "../utils/header/InstanceSingleWindow.hpp"
#include "../utils/header/PrintVector.hpp"

using namespace std;

namespace operations_research
{
    class DownlinkCountDecisionBuilder : public DecisionBuilder
    {
    public:
        explicit DownlinkCountDecisionBuilder(
                                    vector<InstanceSingleWindow *> ins, 
                                    vector<vector<IntVar *>> &prio, 
                                    vector<vector<IntVar *>> &mem, 
                                    IntVar *maxpeak) : 
                                    ins_(ins), 
                                    prio_(prio), 
                                    mem_(mem), 
                                    maxpeak_(maxpeak), 
                                    N_(prio[0].size()), 
                                    W_(prio.size()), 
                                    dc_heuristic_(DownlinkCount(ins)), 
                                    generator_(std::random_device{}())
        {
            memory_state_.resize(N_);
            std::random_device rd;                                          // Obtain a random seed from the hardware
            std::mt19937_64 gen(rd());                                      // Initialize the random number generator engine
            std::uniform_real_distribution<double> dis(min_, max_);
            treshold_ = dis(gen)*maxpeak_->Min()*PRECISION_SOLVE;           // Get random number in range min_ max_

            for(size_t b(0); b < N_; ++b)
            {
                memory_state_[b] = mem_[0][b]->Min();                       // Pick mem min as memory state
            }

            dc_heuristic_.GetSolution(memory_state_, treshold_, 0, current_solution_);

            ranking_ = {};
            // Get ranking from the solution: we want to assign the more prioritary buffers first.
            for(auto line : current_solution_)
            {
                for(auto item : line)
                {
                    ranking_.push_back(item);
                }
            }
        }


        /*
         *  Find next priority variables assigment using DownlinkCount heuristic or first unbound with min value on backtrack
         */
        Decision* Next(Solver* const solver) override
        {
            int window_index(0);
            int buffer_index(0);

            // Find next unbound variable            
            while(window_index < W_ && prio_[window_index][ranking_[buffer_index]]->Bound())
            {
                ++buffer_index;
                if(buffer_index == N_)
                {
                    buffer_index = 0;
                    ++window_index;
                }
            }

            // All variables are bounded, bottom of the search tree, no more decision to take
            if(window_index == W_)
            {
                return nullptr;
            }

            /* 
             *  Compute next solution for the next window
             */
            if(buffer_index == 0)
            {
                // std::cout << "Computing new solution and ranking" << std::endl;
                std::random_device rd;                                          // Obtain a random seed from the hardware
                std::mt19937_64 gen(rd());                                      // Initialize the random number generator engine
                std::uniform_real_distribution<double> dis(min_, max_);
               
                // std::cout << "current rmax: " << maxpeak_-> Max() << std::endl;
                treshold_ = dis(gen)*maxpeak_->Max()*PRECISION_SOLVE;           // Get random number in range min_ max_
                // std::cout << "Treshold: " << treshold_ << std::endl;

                for(size_t b(0); b < N_; ++b)
                {
                    memory_state_[b] = mem_[window_index][b]->Min();            // Pick mem min as memory state
                }
                // std::cout << "mem_state: ";
                // print_vec(memory_state_);
                dc_heuristic_.GetSolution(memory_state_, treshold_, window_index, current_solution_);

                ranking_ = {};
                // Get ranking from the solution: we want to assign the more prioritary buffers first.
                for(auto line : current_solution_)
                {
                    for(auto item : line)
                    {
                        ranking_.push_back(item);
                    }
                }
                // print_vec(ranking_);
                // print_2Dvec(current_solution_);
            }

            // Make actual decision
            int priority = GetPriorityFromSolution(ranking_[buffer_index], current_solution_);   // Get the priority to assign

            if(priority > prio_[window_index][ranking_[buffer_index]]->Max())
            {
                priority = prio_[window_index][ranking_[buffer_index]]->Max();
            } 
            else if(priority < prio_[window_index][ranking_[buffer_index]]->Min())
            {
                priority = prio_[window_index][ranking_[buffer_index]]->Min();
            }
            else
            {
                while(!prio_[window_index][ranking_[buffer_index]]->Contains(priority))              // Find the closest priority value in the current variable domain
                    priority = (priority + 1) % N_;
            }

            // std::cout << "Decision: p[" << window_index << "][" << ranking_[buffer_index] << "]=" << priority << std::endl;
            return solver->MakeAssignVariableValue(prio_[window_index][ranking_[buffer_index]], priority);
        }

    private:
        DownlinkCount dc_heuristic_;           // DownlinkCount heuristic
        vector<InstanceSingleWindow *> ins_;   // Instance
        vector<vector<IntVar *>> mem_;         // Memory variables
        vector<vector<IntVar *>> prio_;        // Priority variables
        IntVar *maxpeak_;                      // Heighest peak variable
        vector<vector<int>> current_solution_; // Downlink count solution for the current window
        vector<int64_t> memory_state_;         // Memory state for DC heuristic
        vector<int> ranking_;                  // Ranking returned from the downlink count heuristic 
        int N_;                                // Number of buffers
        int W_;                                // Number of windows
        const double min_ = 0.5;               // Minimum value of the range
        const double max_ = 2;                 // Maximum value of the range
        double treshold_ = 0;                  // Value of the treshold (what we considerer beeing a buffer overflow) 
        const double PRECISION_SOLVE = 0.001;  // How much do we approximate the peaks
        const double PRECISION_INSTANCE = 1.0; // How do we approximate the instance
        std::mt19937 generator_;               // Mersenne Twister engine
    }; 
}

#endif