<h1>Reproduce experimental results presented in "Scheduling Data Transfers in Space Missions with Constraint Programming"</h1>

This branch is dedicated to the reproduction of our experimental results. All instances are available in "instances" folder and the raw results in the "results" folder. However, you can regenerate those results by running the experiments again in this repository.

<h2>Installation</h2>

* Install OR-Tools for C++ : https://developers.google.com/optimization/install/cpp
* Install Cmake 3.18 or above : https://cmake.org/download/
* Install boost library : 

```bash
sudo apt-get install libboost-all-dev
```

* Edit CMakeLists.txt and change or-tools path to your own path :

```c
list(APPEND CMAKE_PREFIX_PATH "/home/jrouzot/or-tools/cmake")
// ===>
list(APPEND CMAKE_PREFIX_PATH "<your/path/to/or-tools>")
```

<h2>Usage</h2>

In the root folder : 

The first time :

```bash
mkdir build
cmake . -B build
```

then :

```bash
cmake --build build
./build/bin/global-constraint-transfert-scheduling -h
```

<h2>Options</h2>

All options expect a number (0/1 for bool options, any integer otherwise)

| **Option** | **Description** | **Default** |
|------------|-----------------|--------------|
| -d (bool)  | Enables debug mode | 0 |
| -t         | Time limit (ms) | 36000000 (1h) |
| --iterativeleveling (bool) | Solve the instance with iterated leveling heuristic              | 0 |
| --repairdescent (bool) | Solve the instance with repair descent heuristic                | 0 |
| --downlink_count (bool) | Enables downlink count branching strategy | 1 |
| --firstfail (bool) | Enables min-dom branching strategy | 0 |
| --random (bool) | Enables random branching strategy | 0 |
| --global (bool) | Enables simulation lower bound propagator | 1 |
| --single_window (bool) | Enables single window propagator | 1 |
| --dense_ranking (bool) | Enables dense ranking constraint | 1 |
| --symmetry (bool) | Enables priority symmetry propagator | 1 |

<h2>Run experiments</h2>

We used [slurm](https://slurm.schedmd.com/documentation.html) for running our experiments. All the scripts for reproducing our experiments are stored in the slurm-scripts folder. We describe the use of each slurm script below:

<ul>
    <li>run-baseline.sh: Run the CP model with all propagators and constraints (but Simulation Check) disabled for 20 small instances</li>
    <li>run-all.sh: Run the CP model for 20 small instances</li>
    <li>run-no-slb.sh: Run the CP model with Simulation Lower Bound disabled for 20 small instances</li>
    <li>run-no-dr.sh: Run the CP model with Dense Ranking disabled for 20 small instances</li>
    <li>run-no-ps.sh: Run the CP model with Priority Symmetry disabled for 20 small instances</li>
    <li>run-no-sw.sh: Run the CP model with Single Window for 20 small instances</li>
    <li>run-small.sh: Run the CP on 240 small instances</li>
    <li>run-large.sh: Run the CP model on 240 large instances</li>
    <li>run-mtp.sh: Run the CP model on the 4 MTP instances</li>
    <li>run-iterative-small.sh: Run IterativeLeveling heuristic on 240 small instances</li>
    <li>run-iterative-large.sh: Run IterativeLeveling heuristic on 240 large instances</li>
    <li>run-repair-small.sh: Run RepairDescent heuristic on 240 small instances</li>
    <li>run-repair-small.sh: Run RepairDescent heuristic on 240 large instances</li>
</ul>