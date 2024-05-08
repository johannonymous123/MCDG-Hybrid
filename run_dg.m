                                               %Read in input

%%Everything till line 76 sets up different cases (line, lattice) for
%%parallel execution
pars = combvec_no_add_ons_needed(n_cell,n_mus); 

switch mode                     
    case 'line'
        analytical_solution = @(n_x) line_source_solution(201,false); 
        CFL = 0.5;
        length = 3; 
    case 'lattice'
        Phi_analytical_out = readmatrix('lattice.txt');
        L2_ana_out = sqrt(sum((Phi_analytical_out).^2,'all')*49/448^2);
        CFL = 25.6;
        length = 7;
    case 'hohlraum'
        Phi_analytical_out = readmatrix('hohlraum.txt');
        L2_ana_out = sqrt(sum((Phi_analytical_out).^2,'all')*1.3^2/416^2);
        CFL = 52;
        length = 1.3;
    otherwise
        return
end
plot_bool = true; 
dg_solution = @(dt,m,n) dg(dt,m,m,n,mode,plot_bool);  %dg() is actual solver




N_jobs = size(pars,2);
for i = 1:size(pars,2)                             %Comment for and use parfor for multiple runs
%parfor i = 1:size(pars,2)
disp(strcat(['DG:',mode,'   ',num2str(i),' of ',num2str(N_jobs),'.    started at:',datestr(datetime('now'))]))

    n_x = pars(1,i);
    n_mu = pars(2,i);
    dt = CFL*length/n_x;


   [Phi_dg,run_time,moves,moves_tot,iterations] = dg_solution(dt,n_x,n_mu);
        switch mode
        case 'line'
            Phi_analytical = analytical_solution(n_x);
          

        case 'lattice'
            Phi_analytical = Phi_analytical_out;
           
        case 'hohlraum'
            Phi_analytical = Phi_analytical_out;
                
    end
  
     disp(strcat(['DG:',mode,'   ',num2str(i),' of ',num2str(N_jobs),'    done!']))   

switch mode
        case 'line'
            Phi_analytical = analytical_solution(n_x);
            if plot_bool
            figure
            surf(Phi_analytical)
            colorbar
            caxis([0,0.5])
            xlabel('x')
            ylabel('y')
            title('Reference')
            shading interp
            end
          
        case 'lattice'
            Phi_analytical = Phi_analytical_out;
            if plot_bool
            figure
            surf(log10(max(Phi_analytical,1e-15)))
            colorbar
            caxis([-6,-1])
            xlabel('x')
            ylabel('y')
            title('log(Reference)')
            shading interp
            end
        case 'hohlraum'
            Phi_analytical = Phi_analytical_out;
            if plot_bool
            figure
            surf(log10(max(Phi_analytical,1e-15)))
            colorbar
            caxis([-6,1])
            xlabel('x')
            ylabel('y')
            title('log(Reference)')
            shading interp
            end
end
if N_jobs>1
    disp('Paused before next example: press any button in console to continue')
pause
end
end


delete(gcp('nocreate'));



%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pars = combvec_no_add_ons_needed(varargin)
    %varargin is a comma separated list of input vectors
    %This function generates all combinations of these input vectors as
    %output. Input ([a,b],[c,d])-> [ac,ad;bc,bd]    
    n = nargin;
    pars=[];

    pars = varargin{1};
    for i=2:n
        j=1:size(pars,2);
        [a,b] = meshgrid(varargin{i},j);
        c=cat(2,a',b');
        d=reshape(c,[],2);
        pars = [pars(:,d(:,2)');d(:,1)'];
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%


function [phi,Px,Py] = create_plot(degree_x,degree_y,C,w_mu,delta_x,delta_y,plot_bool)
%Just plots cell averages
plot_points_x = [];
plot_points_y = [];
sub = 1;

In = 0;;
[deg_tot,~,n_x,n_y]=size(C);
 phi = zeros(sub*n_x,sub*n_y);     
for j = (1:n_x)
    for i = (1:n_y)
        C_bar(:,j,i) = (C(:,:,j,i)*w_mu);
    end
end

Legendre_at_zero = @(deg) (mod(deg,2)==0)*(-1)^(deg/2)*prod(1:deg)/2^deg/(prod(1:deg/2))^2;

for dx= 0:degree_x
    for dy = 0:degree_y 
        P(:,:,dy+dx*(degree_y+1)+1) = Legendre_at_zero(dx)*Legendre_at_zero(dy);
     %P(:,:,dy+dx*(degree_y+1)+1) = legendreP(dx,In)'*legendreP(dy,In);
    end
 end
for i=0:n_x-1
    plot_points_x = [plot_points_x, In*delta_x*0.5+delta_x*0.5+delta_x*i];    
end
for i=0:n_y-1
    plot_points_y = [plot_points_y, In*delta_y*0.5+delta_y*0.5+delta_y*i];    
end
[Px,Py]=meshgrid(plot_points_x,plot_points_y);
for i = 0:n_x-1
    for j = 0:n_y-1
        for d = 1:deg_tot
        phi(((1+i*sub):(sub+i*sub)),((1+j*sub):(sub+j*sub))) =  phi(((1+i*sub):(sub+i*sub)),((1+j*sub):(sub+j*sub))) + C_bar(d,i+1,j+1)*P(:,:,d);
        end
    end    
end
if plot_bool
figure()

surf(Px,Py,phi)
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






function [Phi,run_time,moves_cell,moves_tot,iterations] = dg(delta_t,n_x,n_y,n_mu_alt,mode,plot_bool)
%This function is the actual solver: Inputs are time steps, cells in x and
%y direaction, number of ordinates, mode(line,lattice, hohlraum), and TF(plot or not)


n_mu_az = 2*n_mu_alt;   %# of azumuth angles
epsilon = 1e-4;         %source iteration error-threshold    
degree_x = 1;           %degrees of product basis-functions
degree_y = 1;           %

switch mode      %Initialize material parameters depending on modes
    case 'line'
        [sigma_s,sigma_a,q,~,x_left,x_right,y_left,y_right,t_final,~] = Line_Init(delta_t);
        res = n_x;
    case 'lattice'
        [sigma_s,sigma_a,q,~,x_left,x_right,y_left,y_right,t_final,~] = Lattice_Init(delta_t);
        res = 448;    
    case 'hohlraum'
        res = 416;
        [sigma_s,sigma_a,q,~,x_left,x_right,y_left,y_right,t_final,~] = Hohlraum_Init(delta_t);
    otherwise 
        disp('no valid mode')
        return
end


%Get discretization

n_mu = n_mu_az*n_mu_alt/2;    %total number of angles
[mu_alt,w] = lgwt(n_mu_alt,-1,1);  %get ordinates and weighs
mu_alt = mu_alt(1:n_mu_alt/2);     %symmetry in z -directtion
azimuth = (0:n_mu_az-1)*pi/n_mu_alt+pi/2/n_mu_alt;
w = repelem(w(1:(n_mu_alt/2)),n_mu_az);
w = w/sum(w);
[angle1,angle2]= meshgrid(mu_alt,azimuth);
mu = [reshape(angle1,[n_mu,1]), reshape(angle2,[n_mu,1])];
Omega = [cos(mu(:,2)).*sqrt(1-mu(:,1).^2), sin(mu(:,2)).*sqrt(1-mu(:,1).^2)];
delta_x = (x_right-x_left)/n_x;
delta_y = (y_right-y_left)/n_y;
time_steps = ceil(t_final/delta_t);
n_deg = (degree_x+1)*(degree_y+1);
[degree1,degree2] = meshgrid(0:degree_x,0:degree_y);
degrees = [reshape(degree1,[n_deg,1]),reshape(degree2,[n_deg,1])];

%Initialize Matrices used in linear system resulting from sweep
matrix = matrix_init(degree_x,degree_y,delta_x,delta_y);


x_center = (1:n_x)*delta_x-0.5*delta_x;
y_center = (1:n_y)*delta_y-0.5*delta_y;
Sigma_a = sigma_a(x_center,y_center);
Sigma_s = sigma_s(x_center,y_center);

Q = zeros(n_deg,n_mu,n_x,n_y);

for i = 1:n_mu
   Q(1,i,:,:) = reshape(q(Omega(i,:),x_center,y_center),[ 1 1 n_x n_y]);
end

Sigma_t = Sigma_a + Sigma_s + 1/delta_t;

%initialize coefficients
switch mode 
    case 'line'        
        Coef_new = gauss_line_init(x_center,y_center,n_mu,n_deg); %Gaussian init
        sweeping = @(Omega,Q_new, Sigma_t, matrix, n_deg) sweep(Omega,Q_new, Sigma_t, matrix, n_deg);
    case 'hohlraum'
        Coef_new = zeros(n_deg,n_mu,n_x,n_y);
        sweeping = @(Omega,Q_new, Sigma_t, matrix, n_deg) hohl_sweep(Omega,Q_new, Sigma_t, matrix, n_deg);
    otherwise
        sweeping = @(Omega,Q_new, Sigma_t, matrix, n_deg) sweep(Omega,Q_new, Sigma_t, matrix, n_deg);
        Coef_new = zeros(n_deg,n_mu,n_x,n_y);
end




%Time evolution
iterations=0;
tic
for t = 0: (time_steps-1) 
    [Q_new, Q_old] = source(Q,Coef_new,Sigma_s,w,delta_t,delta_x,delta_y,matrix.N);  %Initial source
    Coef_new = sweeping(Omega,Q_new, Sigma_t, matrix, n_deg);                        %Initial coeficients
   
    
    iterate = true;
    iterations=iterations+1;
    while iterate
        Coef_old = Coef_new;
        
        Q_new = re_source(Coef_new,Q_old,Sigma_s,w,matrix.N);  %Regenerate RHS with new Phi
        Coef_new = sweeping(Omega,Q_new, Sigma_t, matrix, n_deg);
        delta_C = max(abs(Coef_new-Coef_old),[],'all');
        iterate = delta_C>epsilon;  
        iterations = iterations+1;
    end
    
    Coef_new = Coef_new;
end
run_time=toc;
iterations = iterations/time_steps;
moves_cell = 4*iterations*n_mu_az^2;
moves_tot = moves_cell*n_x*n_y;
%[Phi,Px,Py] = create_plot(degree_x,degree_y,Coef,w,delta_x,delta_y,plot_bool);
Phi = make_phi(Coef_new,w,n_x,delta_x,delta_y,false);

if plot_bool
    figure
    %surf(Px,Py,Phi)
    if strcmp(mode,'line')
        surf(Phi)
        colobar
        caxis([0,0.5])
    else
        surf(log10(max(Phi,1e-15)))
        colorbar
        if strcmp(mode,'lattice')
        caxis([-6,-1])
        else
        caxis([-6,1])
        end
    end
    xlabel('x')
    ylabel('y')
    title('hybrid')
    shading interp
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Coef = gauss_line_init(x_center,y_center,n_mu,n_deg)
%Generates initial data for line as normal-distribution
n_x=length(x_center);
n_y=length(y_center);
[X,Y]=meshgrid(x_center,y_center);
X=[X(:) Y(:)];
Coef = zeros(n_deg,n_mu,n_x,n_y);
f = @(x,y) exp(-((x-1.5).^2+(y-1.5).^2)/2/0.03^2)/2/pi/0.03^2; 
for i=1:n_mu
    Coef(1,i,:,:) = reshape(f(X(:,1),X(:,2)),[1,1,n_x,n_y]);
end
end  





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [sigma_s,sigma_a,q,Particles,x_left_bound,x_right_bound,y_left_bound,y_right_bound,t_f,N_Q] = Hohlraum_Init(delta_t)
%Sets up material parameters for Hohlraum

N_Q = 1000;
x_left_bound = 0;
x_right_bound = 1.3;
y_left_bound = 0;
y_right_bound = 1.3;
t_f =2.6;



blue = @(x,y) (x<=1.25)'*(y<0.05 | y >1.25)+(x>1.25)'*y.^0;
%
red = @(x,y) (x<0.05)'*(y>0.25 & y<1.05);
green = @(x,y) (x>0.45 & x<0.5)'*(y>0.25 & y<1.05)+(x<0.85 & x >0.5)'*((y<0.3 & y >0.25)|(y < 1.05 & y>1));
black = @(x,y) (x>0.5 & x<.85)'*(y>.3 & y<1);
white = @(x,y) ~black(x,y) & ~red(x,y) & ~green(x,y) & ~blue(x,y);
sigma_a = @(x,y) 1.0*(100*blue(x,y)+5*red(x,y)+10*green(x,y)+50*black(x,y));
sigma_s = @(x,y) 1.0*(0.1*white(x,y)+95*red(x,y)+90*green(x,y)+50*black(x,y));
q = @(om,x,y) 0*red(x,y);
Particles = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function C = hohl_sweep(mu,Q,sigma,matrix,n_deg)
%Sets up sweep for Hohlraum (regular sweep, but includes boundary conditions)
%Sweep solves linear system for coeficients for all angles starting at
%boundaries according to flow diretion.

n_mu = size(mu,1);

[n_x,n_y] = size(sigma);
C = zeros(n_deg,n_mu,n_x,n_y);


for i = 1:n_mu
    if mu(i,1)>0
        if mu(i,2)>0
           %bottom left 
           A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(1,1)*matrix.N+mu(i,1)*matrix.Bppx+mu(i,2)*matrix.Bppy);
           C(:,i,1,1) = A\(Q(:,i,1,1)+mu(i,1)*matrix.Bmpx*zeros(n_deg,1)+mu(i,2)*matrix.Bmpy*[1;zeros(n_deg-1,1)]);%[1/(mu(i,1)*c) ;zeros(n_deg-1,1)]);  %replace zeros with boundary
           for jx=2:n_x 
               %bottom row
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(jx,1)*matrix.N+mu(i,1)*matrix.Bppx+mu(i,2)*matrix.Bppy);
               C(:,i,jx,1) = A\(Q(:,i,jx,1)+mu(i,1)*matrix.Bmpx*C(:,i,jx-1,1)+mu(i,2)*matrix.Bmpy*zeros(n_deg,1));
           end
           for jy = 2:n_y
               %left colum
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(1,jy)*matrix.N+mu(i,1)*matrix.Bppx+mu(i,2)*matrix.Bppy);
               C(:,i,1,jy) = A\(Q(:,i,1,jy)+mu(i,1)*matrix.Bmpx*[1;zeros(n_deg-1,1)]...[1/(mu(i,1)*c) ;zeros(n_deg-1,1)]
                    +mu(i,2)*matrix.Bmpy*C(:,i,1,jy-1));
               for jx = 2:n_x
               %rest    
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(jx,jy)*matrix.N+mu(i,1)*matrix.Bppx+mu(i,2)*matrix.Bppy);
               C(:,i,jx,jy) = A\(Q(:,i,jx,jy)+mu(i,1)*matrix.Bmpx*C(:,i,jx-1,jy)+mu(i,2)*matrix.Bmpy*C(:,i,jx,jy-1));
               end
               
           end    
        end
        if mu(i,2)<0
            A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(1,end)*matrix.N+mu(i,1)*matrix.Bppx-mu(i,2)*matrix.Bmmy);
           C(:,i,1,end) = A\(Q(:,i,1,end)+mu(i,1)*matrix.Bmpx*zeros(n_deg,1)-mu(i,2)*matrix.Bpmy*[1;zeros(n_deg-1,1)]);%[1/(mu(i,1)*c) ;zeros(n_deg-1,1)]);  %replace zeros with boundary
           for jx=2:n_x 
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(jx,end)*matrix.N+mu(i,1)*matrix.Bppx-mu(i,2)*matrix.Bmmy);
               C(:,i,jx,end) = A\(Q(:,i,jx,end)+mu(i,1)*matrix.Bmpx*C(:,i,jx-1,1)-mu(i,2)*matrix.Bpmy*zeros(n_deg,1));
           end
           for jy = (n_y-1):-1:1
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(1,jy)*matrix.N+mu(i,1)*matrix.Bppx-mu(i,2)*matrix.Bmmy);
               C(:,i,1,jy) = A\(Q(:,i,1,jy)+mu(i,1)*matrix.Bmpx*[1;zeros(n_deg-1,1)]...[1/(mu(i,1)*c) ;zeros(n_deg-1,1)]
                   -mu(i,2)*matrix.Bpmy*C(:,i,1,jy+1));
               for jx = 2:n_x
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(jx,jy)*matrix.N+mu(i,1)*matrix.Bppx-mu(i,2)*matrix.Bmmy);
               C(:,i,jx,jy) = A\(Q(:,i,jx,jy)+mu(i,1)*matrix.Bmpx*C(:,i,jx-1,jy)-mu(i,2)*matrix.Bpmy*C(:,i,jx,jy+1));
               end
               
           end  
        end
    end    
        
    if mu(i,1)<0
        if mu(i,2)>0
        A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(end,1)*matrix.N-mu(i,1)*matrix.Bmmx+mu(i,2)*matrix.Bppy);
           C(:,i,end,1) = A\(Q(:,i,end,1)-mu(i,1)*matrix.Bpmx*zeros(n_deg,1)+mu(i,2)*matrix.Bmpy*zeros(n_deg,1));  %replace zeros with boundary
           for jx=(n_x-1):-1:1 
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(jx,1)*matrix.N-mu(i,1)*matrix.Bmmx+mu(i,2)*matrix.Bppy);
               C(:,i,jx,1) = A\(Q(:,i,jx,1)-mu(i,1)*matrix.Bpmx*C(:,i,jx+1,1)+mu(i,2)*matrix.Bmpy*zeros(n_deg,1));
           end
           for jy = 2:n_y
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(end,jy)*matrix.N-mu(i,1)*matrix.Bmmx+mu(i,2)*matrix.Bppy);
               C(:,i,end,jy) = A\(Q(:,i,end,1)-mu(i,1)*matrix.Bpmx*zeros(n_deg,1)+mu(i,2)*matrix.Bmpy*C(:,i,end,jy-1));
               for jx = (n_x-1):-1:1
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(jx,jy)*matrix.N-mu(i,1)*matrix.Bmmx+mu(i,2)*matrix.Bppy);
               C(:,i,jx,jy) = A\(Q(:,i,jx,jy)-mu(i,1)*matrix.Bpmx*C(:,i,jx+1,jy)+mu(i,2)*matrix.Bmpy*C(:,i,jx,jy-1));
               end               
           end      
        end
        
        
        
        if mu(i,2) <0
            A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(end,end)*matrix.N-mu(i,1)*matrix.Bmmx-mu(i,2)*matrix.Bmmy);
           C(:,i,end,end) = A\(Q(:,i,end,end)-mu(i,1)*matrix.Bpmx*zeros(n_deg,1)-mu(i,2)*matrix.Bpmy*zeros(n_deg,1));  %replace zeros with boundary
           for jx=(n_x-1):-1:1 
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(jx,end)*matrix.N-mu(i,1)*matrix.Bmmx-mu(i,2)*matrix.Bmmy);
               C(:,i,jx,end) = A\(Q(:,i,jx,end)-mu(i,1)*matrix.Bpmx*C(:,i,jx+1,end)-mu(i,2)*matrix.Bpmy*zeros(n_deg,1));
           end
           for jy = (n_y-1):-1:1
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(end,jy)*matrix.N-mu(i,1)*matrix.Bmmx-mu(i,2)*matrix.Bmmy);
               C(:,i,end,jy) = A\(Q(:,i,end,jy)-mu(i,1)*matrix.Bpmx*zeros(n_deg,1)-mu(i,2)*matrix.Bpmy*C(:,i,end,jy+1));
               for jx = (n_x-1):-1:1
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(jx,jy)*matrix.N-mu(i,1)*matrix.Bmmx-mu(i,2)*matrix.Bmmy);
               C(:,i,jx,jy) = A\(Q(:,i,jx,jy)-mu(i,1)*matrix.Bpmx*C(:,i,jx+1,jy)-mu(i,2)*matrix.Bpmy*C(:,i,jx,jy+1));
               end               
           end      
        end
    end

    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [sigma_s,sigma_a,q,Particles,x_left_bound,x_right_bound,y_left_bound,y_right_bound,t_f,N_Q] = Lattice_Init(delta_t)
%Initializes material parameters for Lattice-problem

N_Q = 1000;
x_left_bound = 0;
x_right_bound = 7;
y_left_bound = 0;
y_right_bound = 7;
t_f = 3.2;

grey = @(x,y) (~((floor(x)==5)'*(floor(y)==3))).*((x>1 & x<6)'*(y>1 & y<6)).*(~(mod(floor(x')+floor(y),2))).*~(((floor(x)==3)'*(floor(y)==3)));
%
red = @(x,y) ((floor(x)==3)'*(floor(y)==3));

sigma_a = @(x,y) 10*grey(x,y);
sigma_s = @(x,y) 1.0*~grey(x,y)+grey(x,y);
q = @(om,x,y) red(x,y);
Particles = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x,w]=lgwt(N,a,b)
%%Used to generate ordinates
% lgwt.m
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
%url: https://www.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes
N=N-1;
N1=N+1; N2=N+2;
xu=linspace(-1,1,N1)';
% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);
% Derivative of LGVM
Lp=zeros(N1,N2);
% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method
y0=2;
% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    
    
    L(:,1)=1;
    Lp(:,1)=0;
    
    L(:,2)=y;
    Lp(:,2)=1;
    
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
 
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    
    y0=y;
    y=y0-L(:,N2)./Lp;
    
end
% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;      
% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [sigma_s,sigma_a,q,Particles,x_left_bound,x_right_bound,y_left_bound,y_right_bound,t_f,N_Q]=Line_Init(delta_t)
%Sets up material parameters for Line-problem
t_f = 1;
x_left_bound = 0;
x_right_bound = t_f*3;
y_left_bound = 0;
y_right_bound = t_f*3;

sigma_s = @(x,y) x'.^0*y.^0*1;
sigma_a = @(x,y)  0*x'.^0*y.^0;
q = @(om,x,y) 0*x'*y;
Particles=[];
N_Q=0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Phi = make_phi(C,w_mu,n_high,dx,dy,plot_bool)
%Averages Coefficients along angles
[deg_tot,~,n_x,n_y]=size(C);
 %phi = zeros(sub*n_x,sub*n_y);     %phi(x)
for j = (1:n_x)
    for i = (1:n_y)
        C_bar(:,j,i) = (C(:,:,j,i)*w_mu);
    end
end

n_c = n_high/n_x;
xp = -dx/2 + 1/n_c*dx*(.5+(0:n_c-1));
[x,y] =  meshgrid(xp,xp);
cc = x.^0;
cl = x.^0.*y;
lc = x.*y.^0;
ll = x.*y;
Phi =  kron(reshape(C_bar(1,:,:),[n_x,n_y]),cc')... 
      + kron(reshape(C_bar(2,:,:),[n_x,n_y]),cl')...
      + kron(reshape(C_bar(3,:,:),[n_x,n_y]),lc')...
      + kron(reshape(C_bar(4,:,:),[n_x,n_y]),ll');

if plot_bool
figure

surf(Phi)
shading interp
colormap jet
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    function matrix = matrix_init(deg_x,deg_y,delta_x,delta_y)
%Sets up matrices resulting from Galerkin-method based on Lagrange
%polynomials


n_deg = (deg_x+1)*(deg_y+1);

%%N
Ny = diag(1./(2*((0:deg_y))+1));
N = Ny;
for i=1:deg_x
   N = blkdiag(N,1/(2*i+1)*Ny);
end
matrix.N=N*delta_y*delta_x;
%% M
[columns,rows] = meshgrid(0:(deg_y),0:(deg_y));
M =  (rows>columns).*(1-(-1).^(columns+rows));
My = M;
for i=1:deg_x
   My = blkdiag(My,1/(2*i+1)*M);
end
matrix.My = delta_x*My;

Mx = zeros(n_deg);



M = diag(1./(2*((0:deg_y))+1))*delta_y;
for i = 0:deg_x
       for j = 0:deg_x
       Mx((i*(deg_y+1)+1):(i*(deg_y+1)+deg_y+1),(j*(deg_y+1)+1):(j*(deg_y+1)+deg_y+1)) = (i>j)*(1-(-1)^(i+j))*M;    
       end
end
matrix.Mx=Mx;


Bppx = zeros(n_deg);
Bpmx = zeros(n_deg);
Bmpx = zeros(n_deg);
Bmmx = zeros(n_deg);
for i = 0:deg_x
       for j = 0:deg_x
       Bppx((i*(deg_y+1)+1):(i*(deg_y+1)+deg_y+1),(j*(deg_y+1)+1):(j*(deg_y+1)+deg_y+1)) = Ny*delta_y;   
       Bpmx((i*(deg_y+1)+1):(i*(deg_y+1)+deg_y+1),(j*(deg_y+1)+1):(j*(deg_y+1)+deg_y+1)) = (-1)^j*Ny*delta_y; 
       Bmpx((i*(deg_y+1)+1):(i*(deg_y+1)+deg_y+1),(j*(deg_y+1)+1):(j*(deg_y+1)+deg_y+1)) = (-1)^i*Ny*delta_y; 
       Bmmx((i*(deg_y+1)+1):(i*(deg_y+1)+deg_y+1),(j*(deg_y+1)+1):(j*(deg_y+1)+deg_y+1)) = (-1)^j*(-1)^i*Ny*delta_y; 
       end
end


App = 1.^columns.*1.^rows;
Apm = 1.^rows.*(-1).^columns;
Amp = (-1).^rows.*1.^columns;
Amm = (-1).^rows.*(-1).^columns;
Bppy = App*delta_x;
Bpmy = Apm*delta_x;
Bmpy = Amp*delta_x;
Bmmy = Amm*delta_x;

for i=1:deg_x
   Bppy = blkdiag(Bppy,1/(2*i+1)*App);
   Bpmy = blkdiag(Bpmy,1/(2*i+1)*Apm);
   Bmpy = blkdiag(Bmpy,1/(2*i+1)*Amp);
   Bmmy = blkdiag(Bmmy,1/(2*i+1)*Amm);
end
matrix.Bppy = Bppy;
matrix.Bpmy = Bpmy;
matrix.Bmpy = Bmpy;
matrix.Bmmy = Bmmy;
matrix.Bppx = Bppx;
matrix.Bpmx = Bpmx;
matrix.Bmpx = Bmpx;
matrix.Bmmx = Bmmx;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






function Q_t = re_source(C,Q,sigma_s,w_mu,N)
%Resets RHS for source iteration, old source+phi from current iteration.
[degree,n_mu,n_x,n_y] = size(Q);
Q_t=zeros(degree,n_mu,n_x,n_y);
C_bar = zeros(degree,n_x,n_y);             %coefficients averaged over mu
degree = degree-1;

 
for j = (1:n_x)
    for i = (1:n_y)
        C_bar(:,j,i) = (C(:,:,j,i)*w_mu);
    end
end
for i = (1:n_mu)
       for j=(1:n_x) 
           for k = (1:n_y)
           Q_t(:,i,j,k) =  Q(:,i,j,k) + sigma_s(j,k)*N*C_bar(:,j,k);
           end
       end
end  
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function [Q_t,Q2] = source(Q,C,sigma_s,w_mu,delta_t,delta_x,delta_y,N)
%Sets up initial RHS for equilibrium equation: RHS= Q+Psi^n/dt
[degree,n_mu,n_x,n_y] = size(Q);
Q_t=zeros(degree,n_mu,n_x,n_x);
Q2=Q_t;
C_bar = zeros(degree,n_x,n_y);             %coefficients averaged over mu
degree = degree-1;





    for j = (1:n_x)
        for i = (1:n_y)
        C_bar(:,j,i) = (C(:,:,j,i)*w_mu);
        end
    end
    D = diag(N);
    for d = 0:degree
        Q2(d+1,:,:,:) = delta_x*delta_y*Q(d+1,:,:,:) + 1/delta_t*C(d+1,:,:,:)*D(d+1);
    end
    for i = (1:n_mu)
       for j=(1:n_x) 
           for k = (1:n_y)
           Q_t(:,i,j,k) =  Q2(:,i,j,k) + sigma_s(j,k)*N*C_bar(:,j,k);
           end
       end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function C = sweep(mu,Q,sigma,matrix,n_deg)
%Solves linear system for coefficient for each cell and each angle
%Sweeps cells starting at boundary with direction best on angle

n_mu = size(mu,1);

[n_x,n_y] = size(sigma);
C = zeros(n_deg,n_mu,n_x,n_y);


for i = 1:n_mu
    if mu(i,1)>0
        if mu(i,2)>0
           %bottom left 
           A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(1,1)*matrix.N+mu(i,1)*matrix.Bppx+mu(i,2)*matrix.Bppy);
           C(:,i,1,1) = A\(Q(:,i,1,1)+mu(i,1)*matrix.Bmpx*zeros(n_deg,1)+mu(i,2)*matrix.Bmpy*zeros(n_deg,1));  %replace zeros with boundary
           for jx=2:n_x 
               %bottom row
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(jx,1)*matrix.N+mu(i,1)*matrix.Bppx+mu(i,2)*matrix.Bppy);
               C(:,i,jx,1) = A\(Q(:,i,jx,1)+mu(i,1)*matrix.Bmpx*C(:,i,jx-1,1)+mu(i,2)*matrix.Bmpy*zeros(n_deg,1));
           end
           for jy = 2:n_y
               %left colum
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(1,jy)*matrix.N+mu(i,1)*matrix.Bppx+mu(i,2)*matrix.Bppy);
               C(:,i,1,jy) = A\(Q(:,i,1,jy)+mu(i,1)*matrix.Bmpx*zeros(n_deg,1)+mu(i,2)*matrix.Bmpy*C(:,i,1,jy-1));
               for jx = 2:n_x
               %rest    
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(jx,jy)*matrix.N+mu(i,1)*matrix.Bppx+mu(i,2)*matrix.Bppy);
               C(:,i,jx,jy) = A\(Q(:,i,jx,jy)+mu(i,1)*matrix.Bmpx*C(:,i,jx-1,jy)+mu(i,2)*matrix.Bmpy*C(:,i,jx,jy-1));
               end
               
           end    
        end
        if mu(i,2)<0
            A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(1,end)*matrix.N+mu(i,1)*matrix.Bppx-mu(i,2)*matrix.Bmmy);
           C(:,i,1,end) = A\(Q(:,i,1,end)+mu(i,1)*matrix.Bmpx*zeros(n_deg,1)-mu(i,2)*matrix.Bpmy*zeros(n_deg,1));  %replace zeros with boundary
           for jx=2:n_x 
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(jx,end)*matrix.N+mu(i,1)*matrix.Bppx-mu(i,2)*matrix.Bmmy);
               C(:,i,jx,end) = A\(Q(:,i,jx,end)+mu(i,1)*matrix.Bmpx*C(:,i,jx-1,1)-mu(i,2)*matrix.Bpmy*zeros(n_deg,1));
           end
           for jy = (n_y-1):-1:1
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(1,jy)*matrix.N+mu(i,1)*matrix.Bppx-mu(i,2)*matrix.Bmmy);
               C(:,i,1,jy) = A\(Q(:,i,1,jy)+mu(i,1)*matrix.Bmpx*zeros(n_deg,1)-mu(i,2)*matrix.Bpmy*C(:,i,1,jy+1));
               for jx = 2:n_x
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(jx,jy)*matrix.N+mu(i,1)*matrix.Bppx-mu(i,2)*matrix.Bmmy);
               C(:,i,jx,jy) = A\(Q(:,i,jx,jy)+mu(i,1)*matrix.Bmpx*C(:,i,jx-1,jy)-mu(i,2)*matrix.Bpmy*C(:,i,jx,jy+1));
               end
               
           end  
        end
    end    
        
    if mu(i,1)<0
        if mu(i,2)>0
        A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(end,1)*matrix.N-mu(i,1)*matrix.Bmmx+mu(i,2)*matrix.Bppy);
           C(:,i,end,1) = A\(Q(:,i,end,1)-mu(i,1)*matrix.Bpmx*zeros(n_deg,1)+mu(i,2)*matrix.Bmpy*zeros(n_deg,1));  %replace zeros with boundary
           for jx=(n_x-1):-1:1 
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(jx,1)*matrix.N-mu(i,1)*matrix.Bmmx+mu(i,2)*matrix.Bppy);
               C(:,i,jx,1) = A\(Q(:,i,jx,1)-mu(i,1)*matrix.Bpmx*C(:,i,jx+1,1)+mu(i,2)*matrix.Bmpy*zeros(n_deg,1));
           end
           for jy = 2:n_y
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(end,jy)*matrix.N-mu(i,1)*matrix.Bmmx+mu(i,2)*matrix.Bppy);
               C(:,i,end,jy) = A\(Q(:,i,end,1)-mu(i,1)*matrix.Bpmx*zeros(n_deg,1)+mu(i,2)*matrix.Bmpy*C(:,i,end,jy-1));
               for jx = (n_x-1):-1:1
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(jx,jy)*matrix.N-mu(i,1)*matrix.Bmmx+mu(i,2)*matrix.Bppy);
               C(:,i,jx,jy) = A\(Q(:,i,jx,jy)-mu(i,1)*matrix.Bpmx*C(:,i,jx+1,jy)+mu(i,2)*matrix.Bmpy*C(:,i,jx,jy-1));
               end               
           end      
        end
        
        
        
        if mu(i,2) <0
            A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(end,end)*matrix.N-mu(i,1)*matrix.Bmmx-mu(i,2)*matrix.Bmmy);
           C(:,i,end,end) = A\(Q(:,i,end,end)-mu(i,1)*matrix.Bpmx*zeros(n_deg,1)-mu(i,2)*matrix.Bpmy*zeros(n_deg,1));  %replace zeros with boundary
           for jx=(n_x-1):-1:1 
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(jx,end)*matrix.N-mu(i,1)*matrix.Bmmx-mu(i,2)*matrix.Bmmy);
               C(:,i,jx,end) = A\(Q(:,i,jx,end)-mu(i,1)*matrix.Bpmx*C(:,i,jx+1,end)-mu(i,2)*matrix.Bpmy*zeros(n_deg,1));
           end
           for jy = (n_y-1):-1:1
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(end,jy)*matrix.N-mu(i,1)*matrix.Bmmx-mu(i,2)*matrix.Bmmy);
               C(:,i,end,jy) = A\(Q(:,i,end,jy)-mu(i,1)*matrix.Bpmx*zeros(n_deg,1)-mu(i,2)*matrix.Bpmy*C(:,i,end,jy+1));
               for jx = (n_x-1):-1:1
               A =  (-mu(i,1)*matrix.Mx-mu(i,2)*matrix.My+sigma(jx,jy)*matrix.N-mu(i,1)*matrix.Bmmx-mu(i,2)*matrix.Bmmy);
               C(:,i,jx,jy) = A\(Q(:,i,jx,jy)-mu(i,1)*matrix.Bpmx*C(:,i,jx+1,jy)-mu(i,2)*matrix.Bpmy*C(:,i,jx,jy+1));
               end               
           end      
        end
    end

    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [M2] = line_source_solution(n_x, plot_bool)
% Evaluates the true solution of the line source test case,
% using Ganapol's quadrature formulas.
t = 1.0;

resolution = 1000;
if nargin ~=2
    plot_res = 501;  % make this odd
    sigma = .03;
    plot_bool = true;
else
    plot_res = n_x;
    sigma = 0.03;
end

rho = [0:(t / resolution):t];
f = phi_l(rho,t);

if plot_bool
set(gca, 'FontName', 'Helvetica');
%plot(rho, f);
print('-depsc2', 'linesource1d.eps');
system('epstopdf linesource1d.eps');
end

dx = t*3.0 / plot_res;
axisLabels = dx * (0:(plot_res - 1))' + dx / 2 - t*3.0 / 2;

if plot_bool
figure;
end
M = zeros(plot_res,plot_res);
for i = 1:plot_res
    for j = 1:plot_res
        x = -3 * t / 2 + (i-1) * 3 * t/(plot_res-1);
        y = -3 * t / 2 + (j-1) * 3 * t/(plot_res-1);
        r = sqrt(x^2 + y^2);
        if r < t
            M(i,j) = interp1(rho, f, r);
        end
    end
end
if plot_bool
set(gca, 'FontName', 'Helvetica');
imagesc(axisLabels, axisLabels, M);
axis square;
caxis([0, 0.7]);
h = colorbar;
set(h, 'FontName', 'Helvetica');
print('-depsc2', 'linesource.eps');
system('epstopdf linesource.eps');
end
total_mass = sum(sum(M)) * (3 * t / (plot_res-1))^2;

if plot_bool
figure;
end
M2 = zeros(plot_res,plot_res);
for i = 1:plot_res
    for j = 1:plot_res
        x = -3 * t / 2 + (i-1) * 3 * t/(plot_res)+dx/2;
        y = -3 * t / 2 + (j-1) * 3 * t/(plot_res)+dx/2;
        r = sqrt(x^2 + y^2);
        value = 1 / (2 * pi * sigma^2) * exp(-r^2 / (2 * sigma^2)) * (3 * t / (plot_res-1))^2;
        iRange  = (max(i - (plot_res-1)/2, 1):min(i + (plot_res-1)/2, plot_res));
        jRange  = (max(j - (plot_res-1)/2, 1):min(j + (plot_res-1)/2, plot_res));
        iRange2 = (max((plot_res+3)/2 - i, 1):min((plot_res+1)/2 + plot_res - i, plot_res));
        jRange2 = (max((plot_res+3)/2 - j, 1):min((plot_res+1)/2 + plot_res - j, plot_res));
        M2(iRange, jRange) = M2(iRange, jRange) + value * M(iRange2, jRange2);
    end
end
if plot_bool
set(gca, 'FontName', 'Helvetica');
imagesc(axisLabels, axisLabels, M2);
axis square;
h = colorbar;
set(h, 'FontName', 'Helvetica');
print('-depsc2', 'gaussian.eps');
system('epstopdf gaussian.eps');
total_mass2 = sum(sum(M2)) * (3 * t / (plot_res-1))^2

figure;
set(gca, 'FontName', 'Helvetica');
plot(axisLabels, M2((plot_res+1)/2, :));
print('-depsc2', 'gaussian1d.eps');
system('epstopdf gaussian1d.eps');
end


end
 
function f = phi_l(rho,t)
eta = rho/t;
ind = find(eta<1);
f = rho*0; phi_l0 = f;
for k = ind
    f(k) = quadv(@(w) phi_pt(t*sqrt(eta(k).^2+w.^2),t),...
    0,sqrt(1-eta(k).^2),1e-5);
end
phi_l0(ind) = exp(-t)/(2*pi*t^2)./sqrt(1-eta(ind).^2);
f = phi_l0+(2*t)*f;
end
 
function f = phi_pt(r,t)
r = max(r,1e-10); % move zero value into small positive regime
eta = r/t;
ind = find(eta<1);
g = r*0;
% compute integral
for k = ind
    g(k) = quadv(@(u) sec(u/2).^2.*real(...
        (eta(k)+1i*tan(u/2)).*...
        xi(eta(k),u).^3.*exp(t/2*(1-eta(k).^2).*xi(eta(k),u)) ...
        ),0,pi,1e-3);
end
% final result
f = 1/(2*pi)*exp(-t)./(4*pi*r.*t.^2)*(t/2)^2.*(1-eta.^2).*g;
f(ind) = f(ind)+phi_pt1(r(ind),t);
end
 
function f = xi(eta,u)
q = (1+eta)./(1-eta);
q = min(max(q,-1e80),1e80); % remove singularity at eta = 1
f = (log(q)+1i*u)./(eta+1i*tan(u/2));
end
 
function y = phi_pt1(r,t)
eta = r/t;
q = (1+eta)./(1-eta);
q = min(max(q,-1e80),1e80); % remove singularity at eta = 1
y = exp(-t)./(4*pi*r*t^2)*t.*log(q);
end
