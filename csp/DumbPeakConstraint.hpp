#ifndef DUMB_PEAK_CONSTRAINT
#define DUMB_PEAK_CONSTRAINT

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
    /*
     *  This constraint links the peak and memory variables with the priority assignment for a single window.
     *  We run the polynomial simulation algorithm to compute the peak memory usage from a given solution.
     *  We don't do any propagation, we just compute the memory and the peak for each buffer for a given start memory and priorities
     */
    class DumbPeakConstraint : public Constraint
    {

    public:
        /*
         *  Instanciate the constraint.
         */
        DumbPeakConstraint(Solver *const solver, vector<IntVar *> &memory_start, vector<IntVar *> &memory_end, vector<IntVar *> &priorities,
                       vector<IntVar *> &peaks, InstanceSingleWindow *&mtp, double precision) :
                       Constraint(solver), memory_start_(memory_start), memory_end_(memory_end), priorities_(priorities), peaks_(peaks),
                       mtp_(mtp), N_(mtp_->getNumBuffers()), precision_(precision)
        {

            algo_ = new AlgorithmSingleWindow<double, double>(*mtp_); // Init algorithm

            // Init priorities_algo_
            for (short int i(0); i < N_; ++i)
            {
                priorities_algo_.push_back({});
            }

            memory_state_algo_.resize(N_); // Init memory_state_algo_
            peaks_state_algo_.resize(N_);  // init peaks_state_algo_
        }

        ~DumbPeakConstraint() override
        {
            delete (algo_);
        }

        /*
         *  Setup the demon for the propagation algorithm of this constraint.
         */
        void Post() override
        {
            Demon *const priority_peak_demon = MakeDelayedConstraintDemon0(
                    solver(), this, &operations_research::DumbPeakConstraint::PriorityPeakPropagator, "priority->peak propagator");

            for (short i(0); i < N_; ++i)
            {
                priorities_[i]->WhenBound(priority_peak_demon);     // Attach to the priority variables
                memory_start_[i]->WhenBound(priority_peak_demon);   // Attach to the memory variables
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
         *  Propagator that link the memory start/priority variables and the peaks usage.
         *  This propagator call simulation algorithm if all priority variables are fixed.
         */
        void PriorityPeakPropagator() {

            // We only compute memory end and peaks when the memory start and priorities are fixed
            if (VariablesAreFixed()) {

                // Clear last priorities
                for (short int i(0); i < N_; ++i) {
                    priorities_algo_[i].clear();
                }

                // Build priorities
                for (short int i(0); i < N_; ++i) {
                    priorities_algo_[priorities_[i]->Value()].push_back(i);
                }

                reverse(priorities_algo_.begin(), priorities_algo_.end()); // Priorities must be reversed for simulation algo

                // Build memory state
                for (short int i(0); i < N_; ++i) {
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
                }
            }
        }

        void InitialPropagate() override {
            PriorityPeakPropagator();
        }


    private:
        vector<IntVar *> memory_start_;             // The memory usage at the start of the window
        vector<IntVar *> memory_end_;               // The memory usage at the end of the window
        vector<IntVar *> priorities_;               // The priority variables
        vector<IntVar *> peaks_;                    // The peak variables for the current window
        InstanceSingleWindow *mtp_;                 // The instance to solve
        int N_;                                     // Number of buffers
        double precision_;                          // Define how much the results are rounded
        AlgorithmSingleWindow<double, double> *algo_;           // Algorithms class
        vector<vector<int>> priorities_algo_;       // Priorities in a format that is readable for algorithm
        vector<double> memory_state_algo_;          // Memory state to be updated by the algo
        vector<double> peaks_state_algo_;           // The peaks to be updated by the algo
    };

}
#endif