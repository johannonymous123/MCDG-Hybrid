%%Input file for Hybrid-runs
%This file allows the user to change the user-defined parameters used in a
%Hybrid run and starts the run afterwards.


mode = 'line'; %Chose test problem. Options: line, hohlraum,lattice

% Following inputs can be list and all combinations  of them will be used.
n_cell = [201];   %number of spatial cells along each dimension 
N_Qs = [1e4 1e5 1e6]; % number of new MC particles per MC step (including relabling)
w_mins = [1e-15 ];  %killing weight
n_mus = [4];        %number of ordinates for collided equation


%Starts hybrid run
cwd = pwd;
cd ../../../../../
run_hybrid;
cd(cwd)