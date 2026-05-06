function [epsilon,dgdN] = LM_parameter(tau_21)

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

    %tau_21 = linspace(0.1e-12,0.3e-12,3);   % 0.3ps   smallest e.g. 1e-15

    tau_sp = 38e-9;     % 38ns
    tau3 = 1 / (1/tau_32 + 1/tau_31 + 1/tau_sp); % lifetime of upper level (should be 1.4ps)
    
    % Radiative and photon lifetimes
    tau_p = 3.36e-12;   % Photon lifetime in cavity (3.36ps)
    N = 48;             % Number of stages
    beta = 2e-3;        % Proportion of spontaneous emission events into the cavity mode
    
    tau_out = 1e-12;    % Electron escape time

    Gamma = 0.32;       % Confinement factor
    c = 3e8;            % Speed of light in vacuum (m/s)
    neff = 3.27;        % Effective refractive index
    c_node = c / neff;  % Reduced speed of light in the medium
    sigma_32 = 1.80e-18; % Stimulated emission cross-section (m^2)
    V = W * L * d * N;  % Cavity volume (m^3)

    K = Gamma * c_node * sigma_32 / V; % Coupling factor


    neff = 3.27;
    c0 = 2.9978e8;
    ng= neff;
    vg = c0/ng;  % group velocity

    % G = 0; % effective gain optical co-efficient
    % Xi*tauN = tau_21;
    % epsilon = G*Xi_times_tauN*V/Gamma;

    % in hamadous model guidance pdf
    epsilon = sigma_32*vg*tau_21; % in m^2
    epsilon = epsilon*1e6; % in cm^3
    disp('first value of epsilon')
    disp(epsilon(1))
    disp('final value of epsilon')
    disp(epsilon(length(epsilon)))
    
    Xi = tau_21/tau_32;
    dgdN = sigma_32*(1-Xi)*1e4;  % cm^2
    disp('first value of a')
    disp(dgdN(1))
    disp('final value of a')
    disp(dgdN(length(dgdN)))
