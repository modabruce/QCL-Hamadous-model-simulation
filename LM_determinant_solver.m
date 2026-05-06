% function [resonance_frequency]= determinant_solver(tau21,N_steady_values)


% install Symbolic Math Toolbox from add-on
syms w

format short

% Q_0 = N_steady_values(1);
% S_0 = N_steady_values(2);

Q_0 = 4.4072e+08;
S_0 = 1.9555e+08;


% Define the parameters
    W = 34e-6; % lateral widths  34 micrometer = 34e-6m 
    L = 1e-3; % lateral length 1mm=1e-3m
    d = 45e-9; % 45nm = 45e-9m
    
    % Delta_I = 1e-3; %1mA
    % Delta_J = Delta_I/(W*L); % density
    
    % e = 1.602e-19; % Charge of an electron (in Coulombs)
    
    % Nonradiative scattering times between 32, 31, 21
    tau_32 = 2.1e-12;   % 2.1ps
    tau_31 = 4.2e-12;   % 4.2ps
    tau_21 = 0.1e-12;   % 0.3ps   smallest e.g. 1e-15
    tau_sp = 38e-9;     % 38ns
    tau_3 = 1 / (1/tau_32 + 1/tau_31 + 1/tau_sp); % lifetime of upper level (should be 1.4ps)
    
    % Radiative and photon lifetimes
    tau_p = 3.36e-12;   % Photon lifetime in cavity (3.36ps)
    N = 48;             % Number of stages
    beta = 2e-3;        % Proportion of spontaneous emission events into the cavity mode
    
    tau_out = 1e-12;    % Electron escape time

    Gamma = 0.32;       % Confinement factor
    c = 3e8;            % Speed of light in vacuum (m/s)
    neff = 3.27;        % Effective refractive index
    c_node = c / neff;  % Reduced speed of light in the medium
    sigma_32 = 1.8e-18; % Stimulated emission cross-section (m^2)
    V = W * L * d * N;  % Cavity volume (m^3)

    %K = Gamma * c_node * sigma_32 / V; % Coupling factor

    ng = neff;
    vg = c/ng;

    % K = Gamma*c'*sigma32/V

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
    -(G_N*G_Q*S_0 + N*beta/tau_sp),     1j*w + 1/tau_p - G_N*G_0 - G_N*G_S*S_0];


% Calculate the determinant of the matrix A
det_A = det(A)*(tau_3*tau_sp*tau_p);   % normalization  !!!!

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
disp('values for w (omega) in Ghz where det(A) = 0:');
disp(resonance_frequency);