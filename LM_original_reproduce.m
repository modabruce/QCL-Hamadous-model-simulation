%% 4 to 2 rate equations groups
%% Q1: V = dLW or V=dLW*N   Answer   V = dLWN
%% Q2: J_input = N*r*Jth    correct
% SL laser ratye equation
% reproduce by fourth order Runge–Kutta method
% change all parameter unit to second and meter
clear;
clc;

r = 1.25; % ratio J/Jth

% Unify unit
% Define the parameters   
W = 34e-6; % lateral widths  34 micrometer = 34e-6m 
L = 1e-3; % lateral length 1mm=1e-3m
d = 45e-9; % 45nm = 45e-9m
e = 1.6e-19; % Charge of an electron (in Coulombs)

% nonradiative scattering times between 32 31 21
tau_32 = 2.1e-12;   % 2.1ps
tau_31 = 4.2e-12;   % 4.2ps
tau_21 = 0.3e-12;   % 0.3ps  1e-15
tau_sp = 38e-9;   % 38ns
tau_3 = 1/(1/tau_32+1/tau_31+1/tau_sp);    %  life time of upper level: should be 1.4ps
% radiative spontaneous relaxation time between level 3 and 2
tau_p = 3.36e-12;    % photon life time in cavity 3.36ps

Num_of_state = 48;       % number of stages

% proportion of spontaneous emission events that emita photon into the cavity mode
% beta = 2e-3;  % ref17 in hamadous original
tau_out = 1e-12;  % electron escape time

Gamma = 0.32; % comfinement factor
c = 3e8;  % unit? 30 cm/ps = 3e8 m/s
neff = 3.27;
c_node = c/neff; % c' = c/neff
sigma_32 = 1.8e-18; % 1.8e-14 cm^2 = 1.8e-18 m^2

% stimulated emission crosssection between the upper and lower levels
V = W * L* d* Num_of_state; % WLd?  m^3

neff = 3.27;
c0 = 2.9978e8;
ng= neff;
vg = c0/ng;

Xi = tau_21/tau_32 + tau_21/tau_sp;

diff_gain = sigma_32*(1-Xi);
epsilon = sigma_32*vg*tau_21;

Ith = e*V/(Gamma*vg*sigma_32*tau_p*tau_3*Num_of_state*(1-Xi));  % threshold current   % with ./N
Jth = Ith/(W*L); % threshold current density
J_in = r*Jth;
I_in = r*Ith;
I = Num_of_state*I_in;
J = Num_of_state*J_in;  % in LM, J = J_input times num of stages(48)

% K = Gamma*c_node*sigma_32/V;  % old K values for 4 rate equations 

disp('diff gain(cm^2)')
disp(diff_gain*1e4)
disp('epsilon(cm^3)')
disp(epsilon*1e6)

A_nr = 1/tau_3;
disp('Anr')
disp(A_nr)

disp('Hom.th.current')
disp(I)

% diff_gain = diff_gain*1e4;  % change unit from cm to m?
% epsilon = epsilon*1e6;

G_N = Gamma*vg*diff_gain/V;
K = Gamma* epsilon/V;

beta = G_N *tau_sp;
%%
% Initial conditions [N, P]
% Nc is carrier number; P is photon number

% unit of diff gain and epsilon
N_initial = [0; 0];  % Initial values
Ntr = 0;
R_sp = 0;  %sp

% Create a time vector 
t0 = 0.1e-12;
h = 1.0e-15; % time step
t_end = 300e-12; % time span 0.1-300ps
t = t0:h:t_end;

% Preallocate the solution vectors
Nc = zeros(1, length(t));
photon = zeros(1, length(t));

% Set initial values
Nc(1) = N_initial(1);
photon(1) = N_initial(2);

% Define the system of RATE equations as a function
fN = @(N3, S) W*L*J/e - A_nr*N3 - G_N *(N3 - Ntr) *S /(1+K*S); 
fP = @(N3, S) G_N *(N3 - Ntr) * S/(1+K*S)+ Num_of_state*beta*N3/tau_sp - S/tau_p+ R_sp;

% Runge-Kutta 4th Order Method
for i = 1:(length(t) - 1)
    % k1 values
    k1_N3 = fN(Nc(i), photon(i));
    k1_N0 = fP(Nc(i), photon(i));
    
    % k2 values
    k2_N3 = fN(Nc(i) + 0.5 * h * k1_N3, photon(i) + 0.5 * h * k1_N0);
    k2_N0 = fP(Nc(i) + 0.5 * h * k1_N3, photon(i) + 0.5 * h * k1_N0);
    
    % k3 values
    k3_N3 = fN(Nc(i) + 0.5 * h * k2_N3, photon(i) + 0.5 * h * k2_N0);
    k3_N0 = fP(Nc(i) + 0.5 * h * k2_N3, photon(i) + 0.5 * h * k2_N0);
    
    % k4 values
    k4_N3 = fN(Nc(i) + h * k3_N3, photon(i) + h * k3_N0);
    k4_N0 = fP(Nc(i) + h * k3_N3, photon(i) + h * k3_N0);
    
    % Update the values for the next step
    Nc(i + 1) = Nc(i) + (h / 6) * (k1_N3 + 2 * k2_N3 + 2 * k3_N3 + k4_N3);
    photon(i + 1) = photon(i) + (h / 6) * (k1_N0 + 2 * k2_N0 + 2 * k3_N0 + k4_N0);
end


% Plot the results in 2×2 plot
figure(1);
subplot(2, 1, 1);
plot(t, Nc/Num_of_state, '-o');
xlabel('Time (t)');
ylabel('Carrier population Q');
title('Q vs Time');
grid on;

subplot(2, 1, 2);
plot(t, photon, '-o');
xlabel('Time (t)');
ylabel('Photon number S');
title('S vs Time');
grid on;

% exportgraphics(gcf, fullFilePath, 'Resolution', 300); % 300 DPI resolution
% test carrier number and photon number
disp('carrier numbers')
disp(Nc(end))
disp('photon numbers')
disp(photon(end))

%%
% Nph_sat = 9.16e8; 
% delta_Nth = 8.5e6;
% % population inversion between upper level and lower level
% 
% delta_N = (W*L*J/e)*tau_3*(1-tau_21/tau_32-tau_21/tau_sp)./(1+Photon./Nph_sat);
% % start value of J is zero
% 
% %delta_N = N3-N2;% definition
% 
% figure(3)
% plot(t*1e12, Photon./Nph_sat, '-o', 'DisplayName', 'norm photon');
% hold on;
% plot(t*1e12, delta_N./delta_Nth, '-s', 'DisplayName', 'norm population inversion');
% xlabel('Time (ps)');
% legend show;
% grid on;