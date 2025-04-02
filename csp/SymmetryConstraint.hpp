#ifndef SYMMETRY_CONSTRAINT
#define SYMMETRY_CONSTRAINT

#include <vector>

#include "ortools/constraint_solver/constraint_solver.h"
#include "ortools/constraint_solver/constraint_solveri.h"

#include "../utils/header/InstanceSingleWindow.hpp"
#include "../utils/header/PrintVector.hpp"

namespace operations_research {

// Utility to print the debug string for each priority var
void printDomainsSymmetry(vector<IntVar*> priorities) {
    for(auto p : priorities) {
        string debug = p->DebugString();
        cout << debug << endl; 
    }
    cout << "---------------------------" << endl;
}


/*
 *  This constraint enable symmetry breaking for the priority domains.
 *  When the higher priority buffers take all the bandwidth, there is two cases for another buffer : 
 *  - The buffer is also in the highest priority pool and share the bandwidth
 *  - The buffer is in a lowest priority pool and will get no bandwidth
 * 
 *  The priority domains can be reduced as any pool with lowest priority will be equivalent. 
 *
 */
class SymmetryConstraint : public Constraint {

    public:

        /*
         *  Instanciate the constraint and performs some sanity checks.
         */
        SymmetryConstraint(Solver * const solver, const std::vector<IntVar *>& priorities, const std::vector<IntVar *>& memory, InstanceSingleWindow * &mtp) :
            Constraint(solver), priorities_(priorities), memory_(memory), mtp_(mtp), N_(mtp_->getNumBuffers()) 
            {

                algo_ = new AlgorithmSingleWindow<double, double>(*mtp_); // Init algorithm
                // Init priorities_algo_
                for (short int i(0); i < N_; ++i)
                {
                    priorities_algo_.push_back({});
                }
                memory_state_algo_.resize(N_);  
                hasBandwidth_.resize(N_);
            }

        /*
         *  Setup the demon for the propagation algorithm of this constraint.
         */
        void Post() override 
        {
            for (int i(0); i < N_; ++i)
            {
                Demon* const symmetry_demon = MakeConstraintDemon1(
                    solver(), this, &operations_research::SymmetryConstraint::SymmetryPropagator, "symmetry propagator", i);
                priorities_[i]->WhenBound(symmetry_demon);         // Attach to the priority variables
            }
        }

        /*
         *  Propagator for symmetry
         */
        void SymmetryPropagator(int i) {

            // std::cout << "Propagate " << i << "..." << std::endl;

            // for(auto p : priorities_) {
            //     std::cout << p->DebugString() << std::endl;
            // }

            for(int i{0}; i < N_; ++i) {
                hasBandwidth_[i] = false;
            }

            // Clear last priorities
            for (short int i(0); i < N_; ++i)
            {
                priorities_algo_[i].clear();
            }

            for (short int i(0); i < N_; ++i)
            {
                priorities_algo_[N_-1 - priorities_[i]->Max()].push_back(i);
                memory_state_algo_[i] = memory_[i]->Min();
            }

            // Run simulation and update memory and peak values and update memory start
            algo_->initialise(); 

            // update peaks and memory_end
            algo_->playDownlinkAndTrack(priorities_algo_, memory_state_algo_, hasBandwidth_);

            // print_vec(hasBandwidth_);

            // Check the highest priority that has bandwidth
            int max = 0;
            for(int i{0}; i < N_; ++i) {
                
                if(hasBandwidth_[i] and priorities_[i]->Max() > max) {
                    max = priorities_[i]->Max();
                }
            } 

            if(max < N_) {
                // Buffers with priority at least @max take all the bandwidth   
                for(auto p : priorities_) {
                    if(p->Max() > max+1) {
                        p->SetMax(max+1);
                    }
                }
            }    

            // for(auto p : priorities_) {
            //     std::cout << p->DebugString() << std::endl;
            // }

        }

        void InitialPropagate() override
        {
            return;
        }


    private: 
        std::vector<IntVar *> priorities_;
        std::vector<IntVar *> memory_;
        InstanceSingleWindow * mtp_;
        int N_; 
        AlgorithmSingleWindow<double, double> *algo_;   // AlgorithmSingleWindows class
        vector<vector<int>> priorities_algo_;           // Priorities in a format that is readable for algorithm
        vector<double> memory_state_algo_;              // Memory state to be updated by the algo
        vector<bool> hasBandwidth_;                     // Keep track of buffer that are allocated bandwidth during simulation
};

}

#endif