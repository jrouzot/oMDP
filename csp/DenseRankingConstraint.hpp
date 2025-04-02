#ifndef DENSE_RANKING_CONSTRAINT
#define DENSE_RANKING_CONSTRAINT

#include <vector>

#include "ortools/base/commandlineflags.h"
#include "ortools/base/init_google.h"
#include "ortools/base/map_util.h"
#include "ortools/constraint_solver/constraint_solver.h"
#include "ortools/constraint_solver/constraint_solveri.h"
#include "ortools/util/bitset.h"

#include "../utils/header/InstanceSingleWindow.hpp"
#include "../utils/header/PrintVector.hpp"
#include "../utils/header/BinaryHeap.hpp"

using namespace std;
using namespace dataflow;

namespace operations_research {


// Utility to print the debug string for each priority var
void printDomainDense(vector<IntVar*> &priorities) {
    for(int i(0); i < priorities.size(); ++i) {
        string debug = priorities[i]->DebugString();
        cout << debug << endl; 
    }
    cout << "---------------------------" << endl;
}


// We want to order the priority var according to their upper bound in the binary heap. 
class IntVarComparator
{
public:
    bool operator()(IntVar* a, IntVar* b) {
        return a->Max() < b->Max();
    }
};


// We want to order the priority var according to their lower bound before to insert them in the binary heap. 
bool priorityByLowerBound(IntVar* a, IntVar* b) {
    return a->Min() < b->Min();
}


/*
 *  Check if the root upper bound is equal to the current priority.
 *  In that case, the variable can be assigned to its upper bound and is removed from the heap.
 *  Returns true if the root upper bound is equal to the current priority.  
 */
bool checkUpperBound(Heap<IntVar*, IntVarComparator> &heap, int current_prioprity) {
    if(heap.empty()) {
        return false;
    }
    if(heap.pick()->Max() == current_prioprity) {
        heap.deleteMin();
        return true;
    }
    return false;
}


/*
 *  Affect the current priority to the root variable.
 *  If this affectation is possible remove the root element and returns true.
 *  Returns false otherwise.
 */
bool affectRoot(Heap<IntVar*, IntVarComparator> &heap, int current_prioprity) {
    if(heap.empty()) {
        return false;
    }
    if(heap.pick()->Contains(current_prioprity)) {
        heap.deleteMin();
        return true;
    }
    return false;
}


/*
 *  Try to build a valid solution for all priority variable values assuming a value for a variable.
 *  If we can find a dense ranking with the prio vars and prio[@index] = @value, we return true, otherwise false.
 */
bool prioSat(vector<IntVar*> &priorities, vector<int> &indexes, int index, int value) {

    IntVarComparator cmp;                                   // Comparator for the heap (order by upper bound)
    Heap<IntVar*, IntVarComparator> heap;                   // Heap structure to process the priorities efficiently
    heap = Heap<IntVar*, IntVarComparator>(cmp);            // Init the heap to handle IntVar*, order them by increasing upper bound 

    auto i = indexes.begin();                               // Iterator on the indexes
    int current_priority = -1;                              // Init current priority to -1 as we increment at the beginning of the loop
    bool affected = false;                                  // Have we found a variable to satisfy the current priority level ? 
    bool target_affected = false;                           // Have we affected @priorities[@index] to @value in this dense ranking ?                    

    while(!(heap.empty() && i == indexes.end())) {          // We loop until the heap is empty AND we check all the prio variables
        ++current_priority;                                 // Increment the current priority to proccess   
                                                            // Add prio var to the heap if they contains current prio (the prio var are ordered by lower bound)
        while(i != indexes.end() && (priorities[*i]->Min() == current_priority)) {
            if(*i != index)                                 // We don't want to add the target variable to the heap as we want to check if a dense ranking is possible with target == value
                heap.add(priorities[*i]);                   // Add var to the heap
            ++i;                                            // Increment the iterator on the prio vars
        }

        if (!target_affected) {
            target_affected = (value == current_priority); // If the current priority is the target @value for variable priorities[@index], we use this variable
            affected = target_affected;                    // The current priority is satisfied is target_affected
        }

        while(checkUpperBound(heap, current_priority))      // Remove all variables that can be set to their upper bound.
            affected = true;
        
        if(!affected)                                       // If we couldn't find a variable for the current priority, we use the heap root
            affected = affectRoot(heap, current_priority);  // Try to affect the current priority to the root.

        if(!affected)                                       // At this point, it means that we can't find a variable for current priority
            return false;
    }

    if (!target_affected)
        target_affected = (value == current_priority + 1);  // The target variable can be the last to be affected

    return target_affected;                                 // At this point we found a value (with regard to the dense ranking constraint) for all variables if target_affected.
}


/*
 *  Reduce the domains of a set of priorities variables.
 *  For each priority variable, we fix the variable to a value in the domain. 
 *  and try to build a solution with the new domains.
 *  If we can't find a solution, it mean that this value must be removed from the domain.
 */
void reduceDomains(vector<IntVar*>& priorities, vector<int> &indexes) {

    std::iota(indexes.begin(), indexes.end(), 0);

    std::sort(indexes.begin(), indexes.end(), [&priorities](int a, int b) {
            return priorities[a]->Min() < priorities[b]->Min();
        });

    for(int i(0); i < priorities.size(); ++i) {
        for(int j(priorities[i]->Min()); j < priorities[i]->Max()+1; ++j) {
            if(!prioSat(priorities, indexes, i, j))
                priorities[i]->RemoveValue(j);
        }
    }
}


/*
 *  This constraint maintains the "left first" priority assignment.
 *
 *  Example : 
 *  3 buffers -> 3 priority level.
 *  {0, 0, 1} is a valid priority assignment.
 *  {0, 2, 2} is not a valid priority assignment (should be {0, 1, 1}).
 */
class DenseRankingConstraint : public Constraint {

    public:

        /*
         *  Instanciate the constraint and performs some sanity checks.
         */
        DenseRankingConstraint(Solver * const solver, const vector<IntVar *>& priorities, const vector<IntVar *>& memory, InstanceSingleWindow * &mtp) :
        Constraint(solver), priorities_(priorities), memory_(memory), mtp_(mtp), indexes_(priorities_.size()), N_(priorities_.size()) {
            for(auto p : priorities_) { // The priorities are bounded by the number of buffers
                CHECK_LE(p->Max(), priorities_.size() - 1); 
                CHECK_GE(p->Min(), 0);
            }
        }

        ~DenseRankingConstraint() override {}

        /*
         *  Setup the demon for the propagation algorithm of this constraint.
         */
        void Post() override {
            Demon* const global_demon = // Create the global demon that bind events on variables
                // solver()->MakeDelayedConstraintInitialPropagateCallback(this);
                solver()->MakeConstraintInitialPropagateCallback(this);
            for(size_t i = 0; i < priorities_.size(); ++i) { 
                priorities_[i]->WhenRange(global_demon); // Attach to all variables
            }
        }


        /*
         *  Propagation method for the priority variables.
         */
        void InitialPropagate() override {
            // std::cout << "DenseRanking.." << std::endl;
            reduceDomains(priorities_, indexes_);
            // for (auto p: priorities_)
            // {
            //     std::cout << p->DebugString() << std::endl;
            // }
        }

    private:
        vector<IntVar *> priorities_;
        vector<IntVar *> memory_;
        InstanceSingleWindow * mtp_;
        int N_; 
        vector<int> indexes_;
};

}

#endif