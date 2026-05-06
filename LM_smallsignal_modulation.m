clear;
clc;

% Define the parameters   
W = 34e-6; % lateral widths  34 micrometer = 34e-6m 
L = 1e-3; % lateral length 1mm=1e-3m
d = 45e-9; % 45nm = 45e-9m

Delta_I = 1e-3; %1mA
Delta_J = Delta_I/(W*L); % density

e = 1.6e-19; % Charge of an electron (in Coulombs)
% nonradiative scattering times between 32 31 21
tau_32 = 2.1e-12;   % 2.1ps
tau_31 = 4.2e-12;   % 4.2ps
tau_21 = 0.3e-12;   % 0.3ps
tau_sp = 38e-9;   % 38ns
tau_3 = 1/(1/tau_32+1/tau_31+1/tau_sp);    %  life time of upper level: should be 1.4ps
% radiative spontaneous relaxation time between level 3 and 2

tau_p = 3.36e-12;    % photon life time in cavity 3.36ps
N = 48;       % number of stages
beta = 2e-3;  % ref17 in hamadous original
% proportion of spontaneous emission events that emita photon into the cavity mode

tau_out = 1e-12;  % electron escape time

Gamma = 0.32; %
c = 3e8;  % unit? 30 cm/ps = 3e8 m/s
neff = 3.27;
c_node = c/neff; % c' = c/neff
sigma_32 = 1.8e-18; % 1.8e-14 cm^2 = 1.8e-18 m^2
% stimulated emission crosssection between the upper and lower levels
V = W*L*d*N; % WLd?  m^3

ng = neff;
vg = c/ng;

% K = Gamma*c'*sigma32/V

Xi = tau_21/tau_32 + tau_21/tau_sp;

diff_gain = sigma_32*(1-Xi);
epsilon = sigma_32*vg*tau_21;

G_N = Gamma*vg*diff_gain/V;
K = Gamma* epsilon/V;

% carrier number on steady states from hamadous_original_reproduce.m
% for J/Jth = 2.5
Q_0 = 5.0138e+08;
S_0 = 7.0616e+07;


G_0 = Q_0/(1+K*S_0);
G_Q = 1/(1+K*S_0);
G_S = -K*Q_0/((1+K*S_0)^2);


% Frequency range from 0.1 GHz to 30 GHz with 128 points
freq = linspace(0.1e9, 100e9, 1280); % Frequency in Hz
A_sol = zeros(2, length(freq)); % Preallocate A matrix

% Loop over each frequency to compute A for each w freq
for idx = 1:length(freq)
    w = 2 * pi * freq(idx); % Angular frequency  ??

    % jw term with imaginary part
    jw = 1j * w;
    
    % replace all 0 to 1e-6

    H = [1j*w + 1/tau_3 + G_N*G_Q*S_0,    G_N*G_0 + G_N*G_S*S_0;
    -(G_N*G_Q*S_0 + N*beta/tau_sp),     1j*w + 1/tau_p - G_N*G_0 - G_N*G_S*S_0];

    % Define vector B for current frequency
    B = [ ...
        W * L * Delta_J / e; ...
        0; ...
    ];

    % Solve for A (A = H \ B)
    A = H \ B;  % More efficient than inv(H) * B ???

    % Store the solution A for each frequency
    A_sol(:, idx) = A;
end


figure(1)
plot(freq / 1e9, 20 * log10(abs(A_sol(2, :))), 'LineWidth', 1.5);  % Plot in GHz and dB scale
set(gca, 'XScale', 'log');  % Set x-axis to logarithmic scale
title('Absolute Value of \Delta N_{ph} vs Frequency');  % Title with subscript formatting
xlabel('Frequency (GHz)');  % Label x-axis
ylabel('Magnitude (dB)');  % Label y-axis in dB
grid on;  % Add grid for easier reading