
This repository contains Matlab code used for the numerical experiments described in https://arxiv.org/abs/2312.04217. We provide code to perform Hybrid runs as well as monolithic SN-runs.

When starting from the file interactive_run.m, the user will be prompted to enter the necessary parameters for the chosen run into the console. The script will then automatically start such a run and return plots of the results along with a reference solution.

After the run, the solution will be stored in the Matlab variable Phi_hybrid_mc for a hybrid run and in Phi_dg for an SN run. The reference solution is stored in the variable Phi_analytical.

Alternatively, user-specified parameters can be entered in the files input_dg.m for an SN-run or input_hybrid.m for a hybrid run. The runs will start automatically after running the corresponding input file. Running the code through the input file allows the user to enter lists of parameters to compare multiple runs at once.

The files run_hybrid.m and run_dg.m contain the actual body of code that will produce solutions, while the files hohlraum.txt and lattice.txt contain the reference solutions mentioned in the paper.

The directory 'paper_runs' contains branching directories leading to various input files reproducing the runs in the paper.

Note that these implementations are far from runtime optimized and can take a substantial amount of time to execute. We recommend experimenting with relatively small values of the parameters at first and continuously increasing their value to get a feel for the expected time spent before moving on to more serious runs. This is particularly true for the input files in 'paper_runs' as these start multiple long runs in sequence, which can take hours.
