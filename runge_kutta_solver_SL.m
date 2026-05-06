%% function to solve Hamadous original rate equations
% Get input current from get_Jth, r= J/Jth
function [Nc, photon ,N_steady] = runge_kutta_solver_SL(r, tau_21, Jth , time_step, t_final)

global W L d Num_of_state tau_32 tau_31 tau_sp beta tau_p tau_out Gamma neff e c0 ng vg sigma_32

    % test values of current r
    % disp(['Current r in SL: ', num2str(r)]);  % Print the current value of r
    % disp(['Current tau21(ps) in SL: ', num2str(tau_21*1e12)]);  % print current value of tau_21
    
    % Jth = 3.38e7; % threshold current density (A/m^2)
    J = r * Jth;  % current density (A/m^2)
    
    e = 1.6e-19; % Charge of an electron (C)
    
    % Other constants
    Gamma = 0.32;  % comfinement factor
    c = 3e8;  % Speed of light (m/s)
    neff = 3.27;  % Effective refractive index
    c_node = c / neff;  % Reduced speed of light
    sigma_32 = 1.8e-18;  % Stimulated emission cross-section (m^2)
    V = W * L * d * Num_of_state;  % Volume (m^3)

    % Nonradiative scattering times
    tau_3 = 1 / (1/tau_32 + 1/tau_31 + 1/tau_sp); % Lifetime of upper level (s)
    
    Xi = tau_21/tau_32 + tau_21/tau_sp;

    diff_gain = sigma_32*(1-Xi);
    epsilon = sigma_32*vg*tau_21;

    J = Num_of_state*J;  % in LM, J = J_input times num of stages(48)

    A_nr = 1/tau_3;
    G_N = Gamma*vg*diff_gain/V;
    K = Gamma* epsilon/V;

%beta = G_N *tau_sp;

% neff_test = 3.318;
% Rout=-vg*2*log((neff_test-1)/(neff_test+1))/L;
% Rloss=vg*2000;
% tau_p = 1/(Rout+Rloss);
% disp('tau_p:')
% disp(tau_p)

% Initial conditions [N, P]
% Nc is carrier number; P is photon number

% unit of diff gain and epsilon
N_initial = [0; 0];  % Initial values
Ntr = 0;
R_sp = 0;  %sp

% Create a time vector 
t0 = 0.1e-12;
h = time_step; % time step
t_end = t_final; % time span 0.1-300ps
t = t0:h:t_end;

% Preallocate the solution vectors
Nc = zeros(1, length(t));
photon = zeros(1, length(t));

% Set initial values
Nc(1) = N_initial(1);
photon(1) = N_initial(2);

% Define the system of RATE equations as a function
fN = @(N3, S) W*L*J/e - A_nr*N3 - G_N *(N3 - Ntr) *S /(1+K*S); 
fP = @(N3, S) G_N *(N3 - Ntr) * S/(1+K*S)+ 1*beta*N3/tau_sp - S/tau_p+ R_sp;

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

    % Steady-state values
    N_steady = [Nc(end), photon(end)];

    % disp('N_steady in SL:');
    % disp(N_steady);
end
