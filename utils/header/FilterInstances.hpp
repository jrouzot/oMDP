#ifndef FILTER_INSTANCES_H
#define FILTER_INSTANCES_H

#include "Instance.hpp"
#include "InstanceSingleWindow.hpp"
#include "Options.hpp"
#include "Solver.hpp"

class FilterInstances
{
private:
    Options opt_;
    Instance<double, double> mtp_;
    vector<InstanceSingleWindow*> instances_;

public:
    FilterInstances(Options opt, Instance<double, double> mtp);
    bool Filter();
    bool FilterLb();
};

FilterInstances::FilterInstances(Options opt, Instance<double, double> mtp) : opt_(opt), mtp_(mtp) {
    // Create the single window instances
    for (short w(0); w <= mtp.numDownlinks(); ++w)
    {
        instances_.push_back(new InstanceSingleWindow(mtp_, opt_.instance_file, w, opt.factor));
    }
}

bool FilterInstances::Filter()
{
    // Instanciate the solver
    mysolver::Solver solver = mysolver::Solver(mtp_, instances_, opt_);
    double obj_dummy = solver.playSolution(solver.dummySolution2D()); 
    double obj_lb = solver.getRealLb()*100.0;

    if(obj_dummy == obj_lb)
    {
        std::cout << "Solved to instance " << opt_.instance_file << " to " << obj_dummy << std::endl;
        std::cout << "Lower bound: " << obj_lb << std::endl;
        return true;
    }
    return false;
}

bool FilterInstances::FilterLb()
{
    // Instanciate the solver
    mysolver::Solver solver = mysolver::Solver(mtp_, instances_, opt_);
    auto res = solver.solveRepairDescent(); 
    double obj = get<1>(res);
    double obj_lb = solver.getRealLb()*100.0;

    if(obj == obj_lb)
    {
        std::cout << "Solved instance " << opt_.instance_file << " to lb" << std::endl;
        return true;
    }
    return false;
}

#endif