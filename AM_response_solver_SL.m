%% Solve small signal modulation with varying steady values and tau_21

function [S] = AM_response_solver_SL(freq, tau_21, N_steady_2level)

global W L d Num_of_state tau_32 tau_31 tau_sp beta tau_p tau_out Gamma neff e c0 ng vg sigma_32

    % check value of steady states
    % -------------------------------
    % monitor
    % disp('N_steady:');
    % disp(N_steady_results);  % Print the steady-state values to check
    
    Delta_I = 1e-3; %1mA
    Delta_J = Delta_I/(W*L); % density
    
    e = 1.602e-19; % Charge of an electron (in Coulombs)
    
    % Nonradiative scattering times between 32, 31, 21
    tau_3 = 1 / (1/tau_32 + 1/tau_31 + 1/tau_sp); % lifetime of upper level (should be 1.4ps)
    
    %c_node = c / neff;  % Reduced speed of light in the medium
    V = W * L * d * Num_of_state;  % Cavity volume (m^3)

Xi = tau_21/tau_32 + tau_21/tau_sp;

diff_gain = sigma_32*(1-Xi);
epsilon = sigma_32*vg*tau_21;

G_N = Gamma*vg*diff_gain/V;
K = Gamma* epsilon/V;

% carrier number on steady states from hamadous_original_reproduce.m
Q_0 = N_steady_2level(1);
S_0 = N_steady_2level(2);

G_0 = Q_0/(1+K*S_0);
G_Q = 1/(1+K*S_0);
G_S = -K*Q_0/((1+K*S_0)^2);

% neff_test = 3.318;
% Rout=-vg*2*log((neff_test-1)/(neff_test+1))/L;
% Rloss=vg*2000;
% tau_p = 1/(Rout+Rloss);
% disp('tau_p:')
% disp(tau_p)

A_sol = zeros(2, length(freq)); % Preallocate A matrix

% Loop over each frequency to compute A for each w freq
for idx = 1:length(freq)
    w = 2 * pi * freq(idx); % Angular frequency 

    H = [1j*w + 1/tau_3 + G_N*G_Q*S_0,    G_N*G_0 + G_N*G_S*S_0;
    -(G_N*G_Q*S_0 + beta/tau_sp),     1j*w + 1/tau_p - G_N*G_0 - G_N*G_S*S_0];

    % Define vector B for current frequency
    B = [ ...
        W * L * Delta_J / e; ...
        0; ...
    ];

    % Solve for A (A = H \ B)
    A = H \ B;  % More efficient than inv(H) * B ???

    % Store the solution A for each frequency
    A_sol(:, idx) = A;   % 48 times fitted
end

% Output the photon number response (N_ph)
    S = A_sol(2, :);  % N_ph is the 4th element of the solution

end