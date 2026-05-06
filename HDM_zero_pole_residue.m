function[pole_frequencies,zero_frequencies, transfer_function, R, P, K] = HDM_zero_pole_residue(freq,tau_21,N_steady_values)

global W L d Num_of_state tau_32 tau_31 tau_sp beta tau_p tau_out Gamma neff e c0 ng vg sigma_32

% Define symbolic variables
syms omega

Delta_I = 1e-3; % 1mA
Delta_J = Delta_I/(W*L); % Current density

tau_3 = 1 / (1/tau_32 + 1/tau_31 + 1/tau_sp); % lifetime of upper level
c_node = c0 / neff;  % Reduced speed of light in the medium
V = W * L * d * Num_of_state;  % Cavity volume (m^3)
K = Gamma * c_node * sigma_32 / V; % Coupling factor

% Extract the steady-state carrier and photon numbers
N_30 = N_steady_values(1);
N_20 = N_steady_values(2);
N_ph0 = N_steady_values(4);

% Frequency range

% Define matrix A(w)
A = [1i*omega + 1/tau_3 + K*N_ph0, -K*N_ph0, 0, +K*(N_30 - N_20);
     -(1/tau_32 + 1/tau_sp) - K*N_ph0, 1i*omega + 1/tau_21 + K*N_ph0, 0, -K*(N_30 - N_20);
     -1/tau_31, -1/tau_21, 1i*omega + 1/tau_out, 0;
     -Num_of_state*K*N_ph0 - Num_of_state*beta/tau_sp, Num_of_state*K*N_ph0, 0, 1i*omega + 1/tau_p - Num_of_state*K*(N_30 - N_20)];

% Right-hand side vector B
B = [W * L * Delta_J / e; 0; 0; 0];

% --- Step 1: Find Poles (roots of det(A) = 0)
det_A = det(A); % Compute determinant of A(w)
normalized_det_A = (det_A * tau_3 * tau_21 * tau_out * tau_p); % Normalize if needed


% with/without simplify ??
%normalized_det_A = (det_A * tau_3 * tau_21 * tau_out * tau_p); % Normalize if needed

pole_solutions = vpasolve(normalized_det_A == 0, omega); % Solve det(A) = 0
pole_frequencies = double(pole_solutions) ; % Convert poles to GHz

% --- Step 2: Find Zeros (roots of numerator of transfer function)
A_inv = A\ B;          % Solve A(w) * X = B
numerator = (A_inv(4)); % Extract and simplify the numerator of photon response
zero_solutions = vpasolve(numerator == 0, omega); % Solve numerator = 0
zero_frequencies = double(zero_solutions); % Convert zeros to GHz

%Display results
disp('hamadous Pole Frequencies (GHz):');
disp(pole_frequencies/ (2 * pi * 1e9));

disp('hamadous Zero Frequencies (GHz):');
disp(zero_frequencies/ (2 * pi * 1e9));

% Define symbolic transfer function
transfer_function = numerator / det_A;

% --- Step 3: Coefficients of numerator and denominator (polynomials in w)
[numerator_poly, numerator_den] = numden(numerator);
[det_A_poly, det_A_den] = numden(det_A);

% Extract coefficients
numerator_coeffs = coeffs(numerator_poly, omega, 'All');
denominator_coeffs = coeffs(det_A_poly, omega, 'All');

% --- Step 4: Perform partial fraction expansion
[R, P, K] = residue(double(numerator_coeffs), double(denominator_coeffs));

% Evaluate magnitude across a frequency range
%w_range = freq;
%T_magnitude = double(subs(abs(transfer_function), w, w_range));

% Display results
% disp('hamadous Poles (w):');
% disp(double(P)/(2 * pi * 1e9));  % Poles (frequencies)
% 
% disp('hamadous Amplitudes (R):');
% disp(double(R));  % Residues (amplitudes)
% 
% disp('hamadous Direct Terms (K):');
% disp(double(K));  % Direct terms (if any)

end