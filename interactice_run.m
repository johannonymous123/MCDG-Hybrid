%This file asks the user to chose parameters for either Hybrid or
%monolithic SN runs and enter them in the console. 
% It checks, if said parameters make sense and starts a corresponding run.


validParams = false;

while ~validParams    
    param1 = input('Enter 1 for a Hybrid run and 2 for a monolithic SN run:');
    if param1~=1 & param1~=2
        disp('Invalid input! Try again.');
    else
        validParams = true;
    end
end
 validParams = false;
while ~validParams
    validParams =true;
    param2 = input('What test problem would you like to run: Enter 1 for Line source, 2 for Lattice, 3 for Hohlraum:');
    switch param2
        case 1
            mode = 'line';
        case 2
            mode = 'lattice';
        case 3
            mode = 'hohlraum';
        otherwise
        validParams=false;
        disp('Invalid input! Try again.');
    end  
end

validParams = false;
while ~validParams
    n_cells = [input('Pick the number of carthesian cells along each dimension. Enter n_cells>0:')];
    if n_cells<1
        disp('Invalid input! Try again.');
    else 
        validParams = true;
    end  
end

validParams = false;
while ~validParams
    if param1 ==1
    n_mu = [input('Pick the number of ordinates for the collided equation. Enter n_mu>3 and even:')];
    else
    n_mu = [input('Pick the number of ordinates. Enter n_mu>4:')];
    end
    if n_mu<4 || mod(n_mu,2)~=0
        disp('Invalid input! Try again.');
    else 
        validParams = true;
    end  
end


if param1 == 1 
    validParams = false;
    while ~validParams
        n_Q = [input('Pick the number particles added in each MC step. Enter n_Q>0:')];
        w_mins = [1e-10];
        if n_Q<1
        disp('Invalid input! Try again.');
        else 
            validParams = true;
        end  
    end
    run_hybrid;
else
    run_dg;
end