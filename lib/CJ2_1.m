function [V, converged, i] = CJ2_1(Ybus, Sbus, V0, ref, pv, pq, mpopt)

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

%% Calculate H_tilde
% M = Lx, x = (inv_L)*M = L \ M
L_inv_M = L \ M;
H_tilde = H - N*L_inv_M;

%% evaluate initial mismatch
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


%% Step 1 from CJ3
dVL = - L \ Q;
Vm(pq) = Vm(pq) + dVL;

V =  Vm .* exp(1j * Va);  
mis = (Sbus(Vm) - V .* conj(Ybus * V)) ./ Vm;
P = real(mis([pv; pq]));
Q = imag(mis(pq));

%% do Newton iterations
% ~ = logic
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;
    %% -----  do Step-1 iteration, update Va  -----
    V_dVa = -H_tilde \ P;
    dVa = V_dVa ./ Vm([pv; pq]);

    %% update Va
    Va([pv; pq]) = Va([pv; pq]) + dVa;
    V = Vm .* exp(1j * Va);

    %% evalute mismatch
    mis = (Sbus(Vm) - V .* conj(Ybus * V)) ./ Vm;
    P = real(mis([pv; pq]));
    Q = imag(mis(pq));

    %% -----  do Step-2 iteration, update Vm  -----
    dVm = - L \ Q;

    %% update Vm
    Vm(pq) = Vm(pq) + dVm;
    V = Vm .* exp(1j * Va);

    %% evalute mismatch
    mis = (Sbus(Vm) - V .* conj(Ybus * V)) ./ Vm;
    P = real(mis([pv; pq]));
    Q = imag(mis(pq));
    
    %% check tolerance
    normP = norm(P, inf);
    normQ = norm(Q, inf);
    if mpopt.verbose > 1
        fprintf('\n  P  %3d   %10.3e   %10.3e', i, normP, normQ);
    end
    if normP < tol && normQ < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nCJ2-1 converged in %d iterations.\n', i);
        end
        break;
    end
end

if mpopt.verbose
    if ~converged
        fprintf('\nCJ2-1 did not converge in %d iterations.\n', i);
    end
end
