This repo contains Matlab-code used for the numerical experiments described in https://arxiv.org/abs/2312.04217
We provide code to do Hybrid runs as well as monolithic SN-runs.

When starting from the file interactive_run.m the user will be prompted to enter the necessary parameters
 for the run of their chosing into the console. The script will start such a run automatically and return
 plots of the results along with a reference solution.

Alternatively user specified paramters can be entered in the files input_dg.m for a SN-run or input_hybrid.m for a hybrid run.
The runs will start automatically after running the corresponding hybrid file. Running the code through the input file allows
the user to enter lists of parameters to compare multiple runs at once.


The files run_hybrid.m and run_dg.m contain the actual body of code that will produce solutions, while the files hohlraum.txt and lattice.txt
 contain the  reference solutions mentioned in the paper.

Note that these implementations are far from runtime optimized and can take a substantial amount of time to execute. We recommend experimenting
 with relatively small values of the parameters at first and continuously increasing their value to get a feeling for the expected time spend, 
before moving on to more serious runs.
