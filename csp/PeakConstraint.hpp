#ifndef PEAK_CONSTRAINT
#define PEAK_CONSTRAINT

#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>

#include "ortools/base/commandlineflags.h"
#include "ortools/base/init_google.h"
#include "ortools/base/map_util.h"
#include "ortools/constraint_solver/constraint_solver.h"
#include "ortools/constraint_solver/constraint_solveri.h"
#include "ortools/util/bitset.h"

#include "../utils/header/InstanceSingleWindow.hpp"
#include "../utils/header/AlgorithmSingleWindow.hpp"
#include "../utils/header/PrintVector.hpp"

using namespace std;
using namespace dataflow;


namespace operations_research
{

// Utility to print the debug string for each priority var
void printDomainsPeak(vector<IntVar*> priorities) {
    for(auto p : priorities) {
        string debug = p->DebugString();
        cout << debug << endl; 
    }
    cout << "---------------------------" << endl;
}

    /*
     *  This constraint links the peak and memory variables with the priority assignment for a single window.
     *  We run the polynomial simulation algorithm to compute the peak memory usage from a given solution.
     *  We run the polynomial singlewindow algorithm to compute a lower bound on the peak for the instance and some priority domain filtering.
     */
    class PeakConstraint : public Constraint
    {

    public:
        /*
         *  Instanciate the constraint.
         */
        PeakConstraint( Solver *const solver, vector<IntVar *> &memory_start, vector<IntVar *> &memory_end, vector<IntVar *> &priorities,
                        vector<IntVar *> &peaks, IntVar* &maxpeak, InstanceSingleWindow *&mtp, bool global_constraint, bool singlewindow_constraint, double precision) : 
                        Constraint(solver), memory_start_(memory_start), memory_end_(memory_end), priorities_(priorities), peaks_(peaks), mtp_(mtp), cnt_(0),
                        maxpeak_(maxpeak), N_(mtp_->getNumBuffers()), prioritytransfer_ct_(global_constraint), singlewindow_ct_(singlewindow_constraint), precision_(precision)
        {
            algo_ = new AlgorithmSingleWindow<double, double>(*mtp_); // Init algorithm

            // Init priorities_algo_
            for (short int i(0); i < N_; ++i)
            {
                priorities_algo_.push_back({});
            }
            peaks_state_algo_.resize(N_);  
            memory_state_algo_.resize(N_);  
            memory_limit_algo_.resize(N_);

            old_prio_min_.resize(N_);
            old_prio_max_.resize(N_);
            old_mem_min_.resize(N_);
            old_mem_max_.resize(N_);

            for (short int i(0); i < N_; ++i)
            {
                old_prio_min_[i] = priorities_[i]->Min();
                old_prio_max_[i] = priorities_[i]->Max();
                old_mem_min_[i] = memory_start_[i]->Min();
                old_mem_max_[i] = memory_end_[i]->Max();
            }
        }

        ~PeakConstraint() override
        {
            delete (algo_);
        }

        /*
         *  Setup the demon for the propagation algorithm of this constraint.
         */
        void Post() override
        {
            for (int i(0); i < N_; ++i)
            {
                Demon* const priority_peak_demon = MakeDelayedConstraintDemon1(
                    solver(), this, &operations_research::PeakConstraint::PriorityPeakPropagator, "priority->peak propagator", i);
                Demon* const memory_peak_demon = MakeDelayedConstraintDemon1(
                    solver(), this, &operations_research::PeakConstraint::MemoryPeakPropagator, "memory->peak propagator", i);
                // TODO: Uncomment for method 1 : always propagate when domain change
                // priorities_[i]->WhenRange(priority_peak_demon);         // Attach to the priority variables
                memory_start_[i]->WhenRange(memory_peak_demon);         // Attach to the memory variables
                // TODO: Uncomment for method 2 : only propagate when a variable is bounded
                priorities_[i]->WhenBound(priority_peak_demon);         // Attach to the priority variables
                // memory_start_[i]->WhenBound(memory_peak_demon);         // Attach to the memory variables
            }
        }

        /*
         *  Run simulation algorithm with priorities and memory corresponding to the best case scenario for buffer @index
         */
        void SimulationOptimistic(int index)
        {
            // std::cout << "Running SimulationOptimistic.." << std::endl;
            // Clear last priorities
            for (short int i(0); i < N_; ++i)
            {
                priorities_algo_[i].clear();
            }

            // Build optimistic priorities
            for (short int i(0); i < N_; ++i)
            {
                if (i != index)
                {
                    priorities_algo_[N_-1 - priorities_[i]->Max()].push_back(i);
                }
                else
                {
                    priorities_algo_[N_-1 - priorities_[i]->Min()].push_back(i);
                }
            }

            for(int i(0); i < N_; ++i)
            {
                memory_state_algo_[i] = memory_start_[i]->Min();
            }

            // std::cout << "Running simulation with:" << std::endl;
            // std::cout << "mem=";
            // print_vec(memory_state_algo_);
            // std::cout << "prio=" << std::endl; 
            // print_2Dvec(priorities_algo_);

            // Run simulation and update memory and peak values and update memory start
            algo_->initialise(); 

            // update peaks and memory_end
            algo_->playDownlink(priorities_algo_, memory_state_algo_, peaks_state_algo_, false);

            // update memory_end lb and peak lb for the buffer @index
            memory_end_[index]->SetMin(max(static_cast<int64_t>(0), round_and_cast_up<double, int64_t>(memory_state_algo_[index])));
            peaks_[index]->SetMin(round_and_cast_up<double, int>(peaks_state_algo_[index] / mtp_->getBufferCapacity(index) * (1 / precision_)));

            // std::cout << memory_end_[index]->DebugString() << std::endl;
            // std::cout << peaks_[index]->DebugString() << std::endl; 
        }

        /*
         *  Run single window algorithm with priorities domains to check feasibility. 
         *  Makes the solver fails if no solution.
         *  SingleWindow algorithm also filter some values of the priority variables to eliminate impossible values according to the current domains. 
         */
        void SingleWindow()
        {
            // std::cout << "Running SingleWindow.." << std::endl;

            maxpeak_algo_ = maxpeak_->Max() * precision_;

            // std::cout << "rmax=" << maxpeak_algo_ << std::endl; 

            for(int i(0); i < N_; ++i)
            {
                memory_state_algo_[i] = memory_start_[i]->Min();
                memory_limit_algo_[i] = memory_end_[i]->Max();
            }



            // std::cout << "mem=";
            // print_vec(memory_state_algo_);

            // for(auto p : priorities_)
            // {
            //     std::cout << p->DebugString() << std::endl;
            // }

            algo_->initialise();

            bool fail = algo_->singleWindowWithVar(memory_state_algo_, memory_limit_algo_, priorities_, maxpeak_algo_, false);

            // std::cout << "fail=" << fail << std::endl;

            if(fail)
            {
                solver()->Fail();
            }

            // For each prio var we keep track of variables that become bounded because of single window
            // If one is bounded it will trigger a new call to single window constraint, that will be useless. 
            // for(auto p : priorities_)
            // {
            //     if(p->Bound() && p->OldMax() != p->Max())
            //     {
            //         ++cnt_;
            //         // std::cout << "Updated " << p->DebugString() << std::endl; 
            //     }
            // }
        }


        /*
         *  Simulation algorithm can be run if priorities and memory start are fixed
         */
        void Simulation()
        {
            // std::cout << "Running Simulation.." << std::endl;
            // Clear last priorities
            for (short int i(0); i < N_; ++i)
            {
                priorities_algo_[i].clear();
            }

            // Build priorities
            for (short int i(0); i < N_; ++i)
            {
                priorities_algo_[N_-1 - priorities_[i]->Value()].push_back(i);
            }

            for(int i(0); i < N_; ++i)
            {
                memory_state_algo_[i] = memory_start_[i]->Value();
            }

            // Initialise simulation
            algo_->initialise();

            // Run simulation
            algo_->playDownlink(priorities_algo_, memory_state_algo_, peaks_state_algo_, false);

            // Update peaks and memory_end
            for (short int i(0); i < N_; ++i)
            {
                memory_end_[i]->SetValue(round_and_cast_up<double, int64_t>(memory_state_algo_[i]));
                peaks_[i]->SetValue(round_and_cast_up<double, int>(peaks_state_algo_[i] / mtp_->getBufferCapacity(i) * (1 / precision_)));
                // std::cout << memory_end_[i]->DebugString() << std::endl;
                // std::cout << peaks_[i]->DebugString() << std::endl;
            }

        }


        /*
         *  Check if all the priority and memory variables are fixed, which mean that we can run simulation algorithm
         */
        bool VariablesAreFixed() {
            for (auto p: priorities_) {
                if (!(p->Bound())) {
                    return false; // Current priority variable is not fixed
                }
            }
            for (auto m: memory_start_) {
                if (!(m->Bound())) {
                    return false; // Current memory variable is not fixed
                }
            }
            return true; // If we reach this state, this mean that all variables are fixed
        }


        /*
         *  Propagator that link the memory variables and the max peak usage.
         *  This propagator call simulation algorithm if all priority and meory variables are fixed.
         */
        void MemoryPeakPropagator(int i)
        {
            // We want to run the simulation only when all variables are fixed
            if(VariablesAreFixed())
            {
                Simulation();
            }   
            else if(prioritytransfer_ct_)
            {

                // TODO : Uncomment for method 1 : propagate all variables
                // for (int k(0); k < N_; ++k)
                // {
                //     SimulationOptimistic(k);
                // }

                // // TODO : Uncomment for method 2 : only propagate current variable (i)
                // if (memory_start_[i]->Min() > memory_start_[i]->OldMin())
                // {
                //     for(int j(0); j < N_; ++j)
                //     {
                //         SimulationOptimistic(j);
                //     }
                // }

                // // TODO : Uncomment for method 3 : do nothing 
                // // ;

                // TODO : Uncomment for method 4 : 
                if (priorities_[i]->Min() > old_prio_min_[i])
                {
                    SimulationOptimistic(i);
                }
                for(int j(0); j < N_; ++j)
                {
                    if (priorities_[j]->Min() > old_prio_min_[j] && priorities_[i] <= priorities_[j])
                    {
                        SimulationOptimistic(j);
                    }
                }
            }
            old_prio_min_[i] = priorities_[i]->Min();
            old_prio_max_[i] = priorities_[i]->Max();
            old_mem_min_[i] = memory_start_[i]->Min();
            old_mem_max_[i] = memory_end_[i]->Max();
        }


        /*
         *  Propagator that link the priority variables and the max peak usage.
         *  This propagator call simulation algorithm if all priority variables are fixed.
         */
        void PriorityPeakPropagator(int i)
        {
            if(VariablesAreFixed())
            {
                Simulation();
            }
            else if(prioritytransfer_ct_)
            {
                // 2. SimulationOptimistic propagator
                
                // TODO : Uncomment for method 2.1 : propagate all variables
                // for (int k(0); k < N_; ++k)
                // {
                //     SimulationOptimistic(k);
                // }

                // TODO : Uncomment for method 2.2 : only propagate variable
                // SimulationOptimistic(i);

                // TODO : Uncomment for method 2.3 : propagate for all variable if best case scenario is worst because of new variable i 
                // SimulationOptimistic(i);
                // for(int j(0); j < N_; ++j)
                // {
                //     // if(i != j && priorities_[i]->Max() <= priorities_[j]->Min() && priorities_[i]->OldMax() > priorities_[j]->Min())
                //     // {
                //     if(i != j && priorities_[i]->Max() <= priorities_[j]->Min())
                //     {
                //         SimulationOptimistic(j);
                //     }  
                // }

                if(priorities_[i]->Min() > old_prio_min_[i])
                {
                    SimulationOptimistic(i);
                }

                for(int j(0); j < N_; ++j)
                {
                    if  (   priorities_[i]->Max() <= priorities_[j]->Max() && old_prio_max_[i] > old_prio_max_[j] ||
                            priorities_[i]->Max() < priorities_[j]->Max() && old_prio_max_[i] == old_prio_max_[j]
                        )
                    {
                        SimulationOptimistic(j);
                    }
                }

                // 1. SingleWindow propagator

                // TODO : Uncomment for method 1.0 : Run single window at every prio var change 
                // SingleWindow(); 

                // TODO : Uncomment for method 1.1 : Run single window at every prio var change and check for redundant calls
                // if(cnt_ > 0)
                //     --cnt_;
                // else
                //     SingleWindow();                  // Run SingleWindow algorithm to see if a priority assignement is possible
                
                // TODO : Uncomment for method 1.2 : Run single window if memory state exceed a certain treshold and check redundant calls
                // std::cout << "Checking single window: memory + fill=" << (memory_start_[i]->Min() + mtp_->getTotalFillInt(i)) / mtp_->getBufferCapacity(i) * (1 / precision_) << ", rmax=" <<  maxpeak_->Max() << std::endl;
            }
            if(singlewindow_ct_)
            {
                // if(cnt_ > 0)
                //     --cnt_;
                // else
                if((memory_start_[i]->Min() + mtp_->getTotalFillInt(i)) / mtp_->getBufferCapacity(i) * (1 / precision_) >= maxpeak_->Max())
                {
                    SingleWindow();                 // Run SingleWindow algorithm to see if a priority assignement is possible
                }
            } 
            old_prio_min_[i] = priorities_[i]->Min();
            old_prio_max_[i] = priorities_[i]->Max();
            old_mem_min_[i] = memory_start_[i]->Min();
            old_mem_max_[i] = memory_end_[i]->Max();
        }

        /*
         *  Propagation method for the peak variable.
         *
         *  We compute the maximum peak and memory usage among the buffers according to a set of fixed priorities.
         */
        void InitialPropagate() override
        {
            if(VariablesAreFixed())
            {
                Simulation();
            }
        }

    private:
        vector<IntVar *> memory_start_;                 // The memory usage at the start of the window
        vector<IntVar *> memory_end_;                   // The memory usage at the end of the window
        vector<IntVar *> priorities_;                   // The priority variables
        vector<int> old_prio_min_;                      // Keep track of the last prio lb 
        vector<int> old_prio_max_;                      // Keep track of the last prio ub
        vector<int> old_mem_min_;                       // Keep track of the last mem lb 
        vector<int> old_mem_max_;                       // Keep track of the last mem ub
        vector<IntVar *> peaks_;                        // The peak variables for the current window
        IntVar * maxpeak_;                              // Max peak variable
        InstanceSingleWindow *mtp_;                     // The instance to solve
        int N_;                                         // Number of buffers
        bool prioritytransfer_ct_;
        bool singlewindow_ct_;
        double precision_;                              // Define how much the results are rounded
        AlgorithmSingleWindow<double, double> *algo_;   // AlgorithmSingleWindows class
        vector<vector<int>> priorities_algo_;           // Priorities in a format that is readable for algorithm
        vector<double> memory_state_algo_;              // Memory state to be updated by the algo
        vector<double> memory_limit_algo_;              // Memory limit for single window
        vector<double> peaks_state_algo_;               // The peaks to be updated by the algo
        double maxpeak_algo_;                           // Max peak in algo format
        int cnt_;                                       // Counter to avoid useless triggers of singlewindow global constraint
    };
}
#endif