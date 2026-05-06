function [resonance_frequency]= determinant_solver_SL(tau_21,N_steady_2levels)

global W L d Num_of_state tau_32 tau_31 tau_sp beta tau_p tau_out Gamma neff e c0 ng vg sigma_32

% install Symbolic Math Toolbox from add-on
syms w

format short

Q_0 = N_steady_2levels(1);
S_0 = N_steady_2levels(2);
    
tau_3 = 1 / (1/tau_32 + 1/tau_31 + 1/tau_sp); % lifetime of upper level (should be 1.4ps)
    
    c_node = c0 / neff;  % Reduced speed of light in the medium
    V = W * L * d * Num_of_state;  % Cavity volume (m^3)

    Xi = tau_21/tau_32 + tau_21/tau_sp;

    diff_gain = sigma_32*(1-Xi);
    epsilon = sigma_32*vg*tau_21;

    G_N = Gamma*vg*diff_gain/V;
    K = Gamma* epsilon/V;

    G_0 = Q_0/(1+K*S_0);
    G_Q = 1/(1+K*S_0);
    G_S = -K*Q_0/((1+K*S_0)^2);

    neff_test = 3.318;
    Rout=-vg*2*log((neff_test-1)/(neff_test+1))/L;
    Rloss=vg*2000;
    tau_p = 1/(Rout+Rloss);
    % disp('tau_p:')
    % disp(tau_p)

% Define the matrix A(w) with symbolic variables
A = [1j*w + 1/tau_3 + G_N*G_Q*S_0,    G_N*G_0 + G_N*G_S*S_0;
    -(G_N*G_Q*S_0 + beta/tau_sp),     1j*w + 1/tau_p - G_N*G_0 - G_N*G_S*S_0];


% Calculate the determinant of the matrix A
det_A = det(A)*(tau_3*tau_sp*tau_p);   % normalization  

% output symbolic solutions
% Display the determinant expression (symbolic)
% disp('Determinant of A(w) is:');
% disp(det_A);

% Solve the determinant equation det(A) = 0 for w
solutions = vpasolve(det_A == 0, w); 

% Display the solutions for w (omega)
% disp('Solutions for w (omega) where det(A) = 0:');
% disp(solutions);

% output resonance
resonance_frequency = double(solutions); % Double To convert symbolic solutions into numerical values
resonance_frequency = resonance_frequency/(2*pi*1e9);

% ---------------------------
% monitoring result
disp('current tau21')
disp(tau_21*1e12)
disp('SL model')
disp('values for w (omega) in Ghz where det(A) = 0:');
disp(resonance_frequency);