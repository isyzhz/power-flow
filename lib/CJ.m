function [V, converged, i] = CJ(Ybus, Sbus, V0,  ref, pv, pq, mpopt)

%% default arguments
if nargin < 7
    mpopt = mpoption;
end

%% options
tol         = mpopt.pf.tol;
max_it      = 50;
lin_solver  = mpopt.pf.nr.lin_solver;

%% initialize
converged = 0;
i = 0;
V = V0;
Va = angle(V);
Vm = abs(V);


%% set up indexing for updating V
npv = length(pv);
npq = length(pq);
j1 = 1;         j2 = npv;           %% j1:j2 - V angle of pv buses
j3 = j2 + 1;    j4 = j2 + npq;      %% j3:j4 - V angle of pq buses
j5 = j4 + 1;    j6 = j4 + npq;      %% j5:j6 - V mag of pq buses

%% evaluate Jacobian
j11 = imag(Ybus([pv;pq], [pv;pq])); 
j12 = -real(Ybus([pv;pq],pq));
j21 = real(Ybus(pq, [pv; pq]));
j22 = imag(Ybus(pq, pq));

%% reduce B matrices
 J = [   j11 j12;
         j21 j22;    ];

%% evaluate F(x0)
% pq,pv - indices
% Sbus: Voltage dependent ZIP load
mis = Sbus(Vm) - V .* conj(Ybus * V) ;
% F of size (nx1)
F_V =  [   real(mis([pv; pq])) ./ Vm([pv;pq]);
        imag(mis(pq)) ./ Vm(pq) ];

%% check tolerance
normF = norm(F_V, inf);

if mpopt.verbose > 1
    fprintf('\n it    max P & Q mismatch (p.u.)');
    fprintf('\n----  ---------------------------');
    fprintf('\n%3d        %10.3e', i, normF);
end
if normF < tol
    converged = 1;
    if mpopt.verbose > 1
        fprintf('\nConverged!\n');
    end
end

%% attempt to pick fastest linear solver, if not specified
if isempty(lin_solver)
    nx = length(F_V);
    if nx <= 10 || have_feature('octave')
        lin_solver = '\';       %% default \ operator
    else    %% MATLAB and nx > 10 or Octave and nx > 2000
        lin_solver = 'LU3';     %% LU decomp with 3 output args, AMD ordering
    end
end




%% do Newton iterations
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;

    %% compute update step
    dx = mplinsolve(J, -F_V, lin_solver);

    %% update voltage
    % j1 = 1, j2 = npv
    % j1:j2 - V angle of pv buses
    % j3:j4 - V angle of pq buses
    % j5:j6 - V mag of pq buses
    if npv
        da_pv = dx(j1:j2) ./ Vm(pv);
        Va(pv) = Va(pv) + da_pv;
    end
    if npq
        da_pq = dx(j3:j4) ./ Vm(pq);
        Va(pq) = Va(pq) + da_pq;
        Vm(pq) = Vm(pq) + dx(j5:j6);
    end
    V = Vm .* exp(1j * Va);
    Vm = abs(V);            %% update Vm and Va again in case
    Va = angle(V);          %% we wrapped around with a negative Vm

    %% evalute F(x)
    mis = Sbus(Vm) - V .* conj(Ybus * V) ;
    F_V =  [   real(mis([pv; pq])) ./ Vm([pv;pq]);
        imag(mis(pq)) ./ Vm(pq) ];


    %% check for convergence
    normF = norm(F_V, inf);
    if mpopt.verbose > 1
        fprintf('\n%3d        %10.3e', i, normF);
    end
    if normF < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nCJ converged in %d iterations.\n', i);
        end
    end
end

if mpopt.verbose
    if ~converged
        fprintf('\nCJ did not converge in %d iterations.\n', i);
    end
end
