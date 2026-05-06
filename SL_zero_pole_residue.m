function[pole_frequencies,zero_frequencies, transfer_function, R, P, K] = SL_zero_pole_residue(freq,tau_21,N_steady_values)

global W L d Num_of_state tau_32 tau_31 tau_sp beta tau_p tau_out Gamma neff e c0 ng vg sigma_32


% Define symbolic variables
syms w

Delta_I = 1e-3; % 1mA
Delta_J = Delta_I/(W*L); % Current density

tau_3 = 1 / (1/tau_32 + 1/tau_31 + 1/tau_sp); % lifetime of upper level
V = W * L * d * Num_of_state;  % Cavity volume (m^3)\

Q_0 = N_steady_values(1);
S_0 = N_steady_values(2);

Xi = tau_21/tau_32 + tau_21/tau_sp;

diff_gain = sigma_32*(1-Xi);
epsilon = sigma_32*vg*tau_21;

G_N = Gamma*vg*diff_gain/V;
K = Gamma* epsilon/V;

G_0 = Q_0/(1+K*S_0);
G_Q = 1/(1+K*S_0);
G_S = -K*Q_0/((1+K*S_0)^2);

%beta = G_N *tau_sp;

% Define the matrix A(w) with symbolic variables
A = [1j*w + 1/tau_3 + G_N*G_Q*S_0,    G_N*G_0 + G_N*G_S*S_0;
    -(G_N*G_Q*S_0 + Num_of_state*beta/tau_sp),     1j*w + 1/tau_p - G_N*G_0 - G_N*G_S*S_0];


% Right-hand side vector B
B = [W * L * Delta_J / e; 0;];

% --- Step 1: Find Poles (roots of det(A) = 0)
det_A = det(A); % Compute determinant of A(w)
normalized_det_A = simplify(det_A * tau_3 * tau_sp * tau_p); % Normalize if needed

pole_solutions = vpasolve(normalized_det_A == 0, w); % Solve det(A) = 0
pole_frequencies = double(pole_solutions); % 

% --- Step 2: Find Zeros (roots of numerator of transfer function)
A_inv = A\ B;          % Solve A(w) * X = B
numerator = (A_inv(2)); % Extract and simplify the numerator of photon response


zero_solutions = vpasolve(numerator == 0, w); % Solve numerator = 0

% Check if solutions exist
if isempty(zero_solutions)
    zero_frequencies = []; % Return an empty array if no zeros are found
    disp('No zeros found for the SL transfer function.');
else
    zero_frequencies = double(zero_solutions) / (2 * pi * 1e9); % Convert zeros to GHz
end

% Display results
disp('SL Pole Frequencies (GHz):');
disp(pole_frequencies/ (2 * pi * 1e9));

disp('SL Zero Frequencies (GHz):');
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
disp('SL Poles (w):');
disp(double(P)/(2 * pi * 1e9));  % Poles (frequencies)

disp('SL Amplitudes (R):');
disp(double(R));  % Residues (amplitudes)

disp('SL Direct Terms (K):');
disp(double(K));  % Direct terms (if any)

% Evaluate magnitude across a frequency range
% w_range = freq;
% T_magnitude = double(subs(abs(transfer_function), w, w_range));