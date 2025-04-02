
#include <numeric>

#include "../header/Options.hpp"

using namespace std;
// using namespace dataflow;

struct argbase {
  virtual ~argbase() {}
  virtual void assign() = 0;
};

template <typename Opt, typename ClapArg, typename E = void>
struct arg : public argbase {
  ClapArg carg;
  Opt &opt;

  template <typename... T>
  arg(TCLAP::CmdLine &cmd, Opt &opt, T &&... args)
      : carg(std::forward<T>(args)...), opt(opt) {
    cmd.add(carg);
  }

  virtual void assign() override { opt = carg.getValue(); }
};

template <typename Opt, typename ClapArg>
struct arg<Opt, ClapArg, typename std::enable_if<std::is_enum<Opt>{}>::type>
    : public argbase {
  ClapArg carg;
  Opt &opt;

  template <typename... T>
  arg(TCLAP::CmdLine &cmd, Opt &opt, T &&... args)
      : carg(std::forward<T>(args)...), opt(opt) {
    cmd.add(carg);
  }

  virtual void assign() override {
    opt =
        static_cast<typename std::remove_reference<Opt>::type>(carg.getValue());
  }
};

struct cmdline {
  TCLAP::CmdLine cmd;
  std::vector<std::unique_ptr<argbase>> args;

  cmdline(const std::string &message, const char delimiter = ' ',
          const std::string &version = "none", bool helpAndVersion = true)
      : cmd(message, delimiter, version, helpAndVersion) {}

  template <typename ClapArg, typename Opt, typename... T>
  void add(Opt &opt, T &&... clapargs) {
    args.emplace_back(std::move(std::make_unique<arg<Opt, ClapArg>>(
        cmd, opt, std::forward<T>(clapargs)...)));
  }

  void parse(int argc, char *argv[]) {
    cmd.parse(argc, argv);
    for (auto &arg : args)
      arg->assign();
  }
};

dataflow::Options dataflow::parse(int argc, char *argv[]) {
  using namespace TCLAP;
  using namespace string_literals;
  cmdline cmd("transfer-scheduling", ' ');

  Options opt;

  cmd.add<UnlabeledValueArg<std::string>>(opt.instance_file, "file",
                                          "instance file", true, "", "string");

  cmd.add<ValueArg<bool>>(opt.debug, "d", "debug", "Displays all logs", false, false, "bool");

  cmd.add<ValueArg<bool>>(opt.dc_heuristic, "", "downlink_count", "Enable downlink count heuristic for branching", false, true,
                         "bool");

  cmd.add<ValueArg<bool>>(opt.firstfail_heuristic, "", "firstfail", "Enable firstfail heuristic for branching", false, false,
                         "bool");
  
  cmd.add<ValueArg<bool>>(opt.random_heuristic, "", "random", "Enable random heuristic for branching", false, false,
                         "bool");
                        
  cmd.add<ValueArg<bool>>(opt.dense_ranking, "", "dense_ranking", "Enable Dense Ranking constraint", false, true,
                         "bool");

  cmd.add<ValueArg<bool>>(opt.global, "", "global", "Enable global constraint for priority transfer", false, true,
                         "bool");

  cmd.add<ValueArg<bool>>(opt.singlewindow, "", "single_window", "Enable singleWindow constraint", false, true,
                         "bool");

  cmd.add<ValueArg<bool>>(opt.symmetry, "", "symmetry", "Enable symettry breaking for priority transfer", false, true,
                          "bool");

  cmd.add<ValueArg<bool>>(opt.lns, "", "lns", "Enables LNS", false, false,
                          "bool");

  cmd.add<ValueArg<int>>(opt.time_limit, "t", "time_limit", "Set the time limit for the solver", false, 3600000,
                          "int");

  cmd.add<ValueArg<int>>(opt.luby_factor, "", "luby", "Scale factor for the luby sequence to perform the restarts", false, 100,
                          "int");

  cmd.add<ValueArg<int>>(opt.randomwalk, "", "randomwalk", "Enable randomwalk", false, false,
                          "bool");

  cmd.add<ValueArg<int>>(opt.iterativeleveling, "", "iterativeleveling", "Enable iterativeleveling", false, false,
                          "bool");

  cmd.add<ValueArg<int>>(opt.repairdescent, "", "repairdescent", "Enable repairdescent", false, false,
                          "bool");

  cmd.add<ValueArg<int>>(
      opt.verbosity, "v", "verbosity",
      "verbosity level (0:silent,1:quiet,2:improvements only,3:verbose", false,
      0, "int");

  cmd.add<ValueArg<int>>(opt.seed, "", "seed", "random seed", false, 12345,
                         "int");

  cmd.add<ValueArg<int>>(opt.runs, "", "runs",
                         "number of randomized runs on each dicho step", false,
                         100, "int");

  cmd.add<ValueArg<double>>(opt.factor, "", "factor",
                            "capacity relaxation factor", false, 1, "double");

  cmd.add<ValueArg<string>>(opt.algorithm, "", "algorithm",
                            "choice of algorithm (roc, dcount, ileveling)",
                            false, "roc", "string");

  cmd.add<ValueArg<double>>(opt.epsilon, "", "epsilon", "epsilon increment",
                            false, 0.00001, "double");

  cmd.add<ValueArg<int>>(opt.heuristic, "", "heuristic", "init heuristic",
                         false, 1, "int");

  cmd.add<ValueArg<int>>(opt.lds_max, "", "lds_max", "Set lds max differences",
                         false, 5, "int");

  cmd.add<ValueArg<bool>>(opt.delta, "", "delta", "Enable delta constraint",
                         false, true, "bool");

  cmd.add<ValueArg<bool>>(opt.warm_start, "", "warm_start", "Enable warm start",
                         false, false, "bool");

  cmd.add<ValueArg<int>>(opt.htime_limit, "", "htime_limit", "Set the time limit for the start heuristic", false, 10000,
                          "int");

  cmd.parse(argc, argv);
  return opt;
}

