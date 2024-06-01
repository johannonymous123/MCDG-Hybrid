%%Input file for DG-runs
%This file allows the user to change the user-defined parameters used in a
%SN run and starts the run afterwards.



mode = 'line';   %Chose test problem. Options: line, hohlraum,lattice
% Following inputs can be list and all combinations  of them will be used.
n_cell = [51];  %# of cells along one dimension (make odd for line, multiples of 26 < 208 for hohlraum, multiples of 28 <  224 for lattice)
n_mus = [ 4 8 16];% 32];   %# of ordinates, should be even and >=4


%Starts SN run
cwd = pwd;
cd ../../../..
run_dg;
cd(cwd)