function [V, converged, i] = CJ3_Vk(Ybus, Sbus, V0, ref, pv, pq, mpopt)

%% default arguments
if nargin < 7
    mpopt = mpoption;
end

%% options
tol         = mpopt.pf.tol;
max_it      = 50;


%% initialize
converged = 0;
i = 0;
V = V0;
Va = angle(V);
Vm = abs(V);


%% evaluate Jacobian
j11 = imag(Ybus([pv;pq], [pv;pq]));
j12 = -real(Ybus([pv;pq],pq));
j21 = real(Ybus(pq, [pv; pq]));
j22 = imag(Ybus(pq, pq));

H = j11;
N = j12;
M = j21;
L = j22; 

% M = Lx, x = (inv_L)*M = L \ M

L_inv_M = L \ M;
H_tilde = H - N*L_inv_M;


%% evaluate dP,dQ
% pq,pv - indices
% Sbus: Voltage dependent ZIP load. 33 bus: constant
mis = (Sbus(Vm) - V .* conj(Ybus * V)) ./ Vm;

P = real(mis([pv; pq]));
Q = imag(mis(pq));

%% check tolerance
normP = norm(P, inf);
normQ = norm(Q, inf);
if mpopt.verbose > 1
    fprintf('\niteration     max mismatch (p.u.)  ');
    fprintf('\ntype   #        P            Q     ');
    fprintf('\n---- ----  -----------  -----------');
    fprintf('\n  -  %3d   %10.3e   %10.3e', i, normP, normQ);
end
if normP < tol && normQ < tol
    converged = 1;
    if mpopt.verbose > 1
        fprintf('\nConverged!\n');
    end
end

%% do Newton iterations
% ~ = logic
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;
    %% -----  do Step-1 iteration, update V_temp  -----
    dVL = - L \ Q;
    V_temp = Vm;
    V_temp(pq) = V_temp(pq) + dVL; 
    
    %% calculate mismatch of Step-2
    Vp = V_temp .* exp(1j * Va);  
    mis_step2 = (Sbus(Vm) - Vp .* conj(Ybus * Vp)) ./ Vm ;
    P_step2 = real(mis_step2([pv; pq])) ; 

    %% -----  do Step-2 iteration, update Va  -----
    V_dVa = -H_tilde \ P_step2;
    dVa = V_dVa./Vm([pv; pq]);
    Va([pv; pq]) = Va([pv; pq]) + dVa;
    

    %% -----  do Step-3 iteration, update Vm  -----
    dVM = -L_inv_M *V_dVa;
    Vm(pq) = V_temp(pq) + dVM;
    V = Vm .* exp(1j * Va);

    %% evalute mismatch
    mis = (Sbus(Vm) - V .* conj(Ybus * V) ) ./ Vm;
    P = real(mis([pv; pq]));
    Q = imag(mis(pq));
    
    %% check tolerance
    normP = norm(P, inf);
    normQ = norm(Q, inf);

    if mpopt.verbose > 1
        fprintf('\n  Q  %3d   %10.3e   %10.3e', i, normP, normQ);
    end
    if normP < tol && normQ < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nCJ3-Vk converged in %d iterations.\n', i);
        end
        break;
    end
end

if mpopt.verbose
    if ~converged
        fprintf('\nCJ3-Vk did not converge in %d iterations.\n', i);
    end
end