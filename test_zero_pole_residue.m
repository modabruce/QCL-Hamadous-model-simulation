clear;
clc;

% Define symbolic variables
syms w

freq = linspace(0.1e9, 100e9, 1280); % Frequency in Hz

% Physical parameters
W = 34e-6; % Width (m)
L = 1e-3;  % Length (m)
d = 45e-9; % Thickness (m)
Num_of_state = 48;  % Number of stages

tau_32 = 2.1e-12;  % Nonradiative lifetime (s)
tau_31 = 4.2e-12;  % Nonradiative lifetime (s)
tau_sp = 38e-9;    % Spontaneous lifetime (s)
beta = 2e-3;  % Proportion of spontaneous emission events
tau_p = 3.36e-12;  % Photon lifetime
tau_21 = 2.0e-12;

% Electron escape time
tau_out = 1e-12;  % (s)

Gamma = 0.32;      % Confinement factor

c0 = 2.9978e8;     % Exact speed of light (m/s)
neff = 3.27;       % Effective refractive index
ng = neff;         % Group refractive index
vg = c0 / ng;      % Group velocity

sigma_32 = 1.8e-18; % Stimulated emission cross-section (m^2)
e = 1.6e-19;        % Charge of an electron (C)

Delta_I = 1e-3; % 1mA
Delta_J = Delta_I/(W*L); % Current density

tau_3 = 1 / (1/tau_32 + 1/tau_31 + 1/tau_sp); % lifetime of upper level
c_node = c0 / neff;  % Reduced speed of light in the medium
V = W * L * d * Num_of_state;  % Cavity volume (m^3)
K = Gamma * c_node * sigma_32 / V; % Coupling factor

% Extract the steady-state carrier and photon numbers
N_30 = 2.1664e+08;
N_20 = 2.0802e+8;
N_ph0 = 1.3562e+08;

% Frequency range

% Define matrix A(w)
A = [1i*w + 1/tau_3 + K*N_ph0, -K*N_ph0, 0, +K*(N_30 - N_20);
     -(1/tau_32 + 1/tau_sp) - K*N_ph0, 1i*w + 1/tau_21 + K*N_ph0, 0, -K*(N_30 - N_20);
     -1/tau_31, -1/tau_21, 1i*w + 1/tau_out, 0;
     -Num_of_state*K*N_ph0 - Num_of_state*beta/tau_sp, Num_of_state*K*N_ph0, 0, 1i*w + 1/tau_p - Num_of_state*K*(N_30 - N_20)];

% Right-hand side vector B
B = [W * L * Delta_J / e; 0; 0; 0];

% --- Step 1: Find Poles (roots of det(A) = 0)
det_A = det(A); % Compute determinant of A(w)
normalized_det_A = simplify(det_A * tau_3 * tau_21 * tau_out * tau_p); % Normalize if needed


% with/without simplify ??
%normalized_det_A = (det_A * tau_3 * tau_21 * tau_out * tau_p); % Normalize if needed



pole_solutions = vpasolve(normalized_det_A == 0, w); % Solve det(A) = 0
pole_frequencies = double(pole_solutions) / (2 * pi * 1e9); % Convert poles to GHz

% --- Step 2: Find Zeros (roots of numerator of transfer function)
A_inv = A\ B;          % Solve A(w) * X = B
numerator = (A_inv(4)); % Extract and simplify the numerator of photon response
zero_solutions = vpasolve(numerator == 0, w); % Solve numerator = 0
zero_frequencies = double(zero_solutions) / (2 * pi * 1e9); % Convert zeros to GHz

% Display results
disp('Pole Frequencies (GHz):');
disp(pole_frequencies);

disp('Zero Frequencies (GHz):');
disp(zero_frequencies);

% Define symbolic transfer function
transfer_function = numerator / det_A;

% --- Step 3: Coefficients of numerator and denominator (polynomials in w)
[numerator_poly, numerator_den] = numden(numerator);
[det_A_poly, det_A_den] = numden(det_A);

% Extract coefficients
numerator_coeffs = coeffs(numerator_poly, w, 'All');
denominator_coeffs = coeffs(det_A_poly, w, 'All');

% --- Step 4: Perform partial fraction expansion
[R, P, K] = residue(double(numerator_coeffs), double(denominator_coeffs));

% Display results
disp('Poles (w):');
disp(double(P)/(2 * pi * 1e9));  % Poles (frequencies)

disp('Amplitudes (R):');
disp(double(R));  % Residues (amplitudes)

disp('Direct Terms (K):');
disp(double(K));  % Direct terms (if any)

% Evaluate magnitude across a frequency range
w_range = freq;
T_magnitude = double(subs(abs(transfer_function), w, w_range));

% Plot magnitude response
figure(1)
set(gca, 'XScale', 'log');
plot(freq / (1e9), 20*log10(abs(T_magnitude)));
set(gca, 'XScale', 'log');
xlabel('Frequency (GHz)');
ylabel('|T(w)|');
title('Transfer Function Magnitude Response');
grid on;