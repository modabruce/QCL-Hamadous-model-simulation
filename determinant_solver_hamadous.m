function [resonance_frequency]= determinant_solver_hamadous(tau21,N_steady_values)

global W L d Num_of_state tau_32 tau_31 tau_sp beta tau_p tau_out Gamma neff e c0 ng vg sigma_32

% install Symbolic Math Toolbox from add-on
syms w

format short


N30 = N_steady_values(1);
N20 = N_steady_values(2);
% N10 = N_steady_values(3);
Nph0 = N_steady_values(4);

    tau3 = 1 / (1/tau_32 + 1/tau_31 + 1/tau_sp); % lifetime of upper level (should be 1.4ps)
    c_node = c0 / neff;  % Reduced speed of light in the medium
    V = W * L * d * Num_of_state;  % Cavity volume (m^3)

    K = Gamma * c_node * sigma_32 / V; % Coupling factor

% Define the matrix A(w) with symbolic variables
A = [1i*w + 1/tau3 + K*Nph0, -K*Nph0, 0, +K*(N30 - N20);
     -(1/tau_32 + 1/tau_sp) - K*Nph0, 1i*w + 1/tau21 + K*Nph0, 0, -K*(N30 - N20);
     -1/tau_31, -1/tau21, 1i*w + 1/tau_out, 0;
     -Num_of_state*K*Nph0 - Num_of_state*beta/tau_sp, Num_of_state*K*Nph0, 0, 1i*w + 1/tau_p - Num_of_state*K*(N30 - N20)];


% Calculate the determinant of the matrix A
det_A = det(A)*(tau3*tau21*tau_out*tau_p);   % normalization  !!!!

det_A = simplify(det_A);
% output symbolic solutions
% Display the determinant expression (symbolic)
% disp('Determinant of A(w) is:');
% disp(det_A);

% Solve the determinant equation det(A) = 0 for w
solutions = vpasolve(det_A == 0, w, 'Random', true); 

% Display the solutions for w (omega)
% disp('Solutions for w (omega) where det(A) = 0:');
% disp(solutions);

% output resonance
resonance_frequency = double(solutions); % Double To convert symbolic solutions into numerical values
resonance_frequency = resonance_frequency/(2*pi*1e9);


% ---------------------------
% monitoring result
disp('current tau21')
disp(tau21*1e12)
disp('hamadous model')
disp('values for w (omega) in Ghz where det(A) = 0:');
disp(resonance_frequency);