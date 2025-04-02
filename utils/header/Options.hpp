

#ifndef __CMDLINE_HPP
#define __CMDLINE_HPP

#include <tclap/CmdLine.h>

using namespace std;

namespace dataflow {

class Options {

public:
  // the actual options
  string cmdline; // for reference
  string instance_file;
  bool debug;
  bool dc_heuristic;
  bool firstfail_heuristic;
  bool random_heuristic;
  bool dense_ranking;
  bool global;
  bool singlewindow;
  bool delta;
  bool symmetry;
  int htime_limit;
  int time_limit;
  int luby_factor;
  bool warm_start;
  bool lns;
  bool randomwalk;
  bool iterativeleveling;
  bool repairdescent;

  // Options for heuristic
  enum verbosity { SILENT = 0, QUIET, NORMAL, YACKING, SOLVERINFO };
  int verbosity;
  string algorithm;
  int seed;
  double factor;
  int runs;
  double epsilon;
  int heuristic;
  int lds_max;

  Options(){};
  Options(const Options &opt)
      : cmdline(opt.cmdline), instance_file(opt.instance_file),
        debug(opt.debug), dc_heuristic(opt.dc_heuristic), firstfail_heuristic(opt.firstfail_heuristic), random_heuristic(opt.random_heuristic),
        dense_ranking(opt.dense_ranking), global(opt.global), singlewindow(opt.singlewindow), symmetry(opt.symmetry),
        time_limit(opt.time_limit), luby_factor(opt.luby_factor), verbosity(opt.verbosity), algorithm(opt.algorithm), seed(opt.seed),
        factor(opt.factor), runs(opt.runs), epsilon(opt.epsilon), randomwalk(opt.randomwalk), iterativeleveling(opt.iterativeleveling), repairdescent(opt.repairdescent),
        heuristic(opt.heuristic), lds_max(opt.lds_max), delta(opt.delta), warm_start(opt.warm_start), htime_limit(opt.htime_limit),
        lns(opt.lns) {}
};

Options parse(int argc, char *argv[]);

static Options no_option;
}

// minischeduler::Options no_option;

#endif // __CMDLINE_HPP