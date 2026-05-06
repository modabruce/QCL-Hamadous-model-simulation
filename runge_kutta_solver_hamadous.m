%% function to solve Hamadous original rate equations
% Get input current from get_Jth, r= J/Jth
function [N3, N2, N1, N0, N_steady] = runge_kutta_solver_hamadous(r, tau_21, Jth, time_step, t_final)

global W L d Num_of_state tau_32 tau_31 tau_sp beta tau_p tau_out Gamma neff e c0 ng vg sigma_32

    % test values of current r
    % disp(['Current r: ', num2str(r)]);  % Print the current value of r
    % disp(['Current tau21(ps): ', num2str(tau_21*1e12)]);  % print current value of tau_21
    
    % Jth = 3.38e7; % threshold current density (A/m^2)
    J = r * Jth;  % current density (A/m^2)
    
    V = W * L * d * Num_of_state;  % Volume (m^3)

    c_node = c0/neff;
    K = Gamma * c_node * sigma_32 / V;  % Constant
    
    tau_3 = 1 / (1/tau_32 + 1/tau_31 + 1/tau_sp); % Lifetime of upper level (s)
    
    % Initial conditions [N3, N2, N1, N0]
    N_initial = [0; 0; 0; 0];  % Initial values
    % when [0; 0; 0; 1e8] beta is not essential
    % when [0; 0; 0; 0] beta is essential

    % Time setup
    t0 = 0.1e-12;  % Initial time (s)
    h = time_step;   % Time step (s)
    t_end = t_final; % End time (s)   longer t_end or shorter h
    t = t0:h:t_end;  % Time vector
    
    % Preallocate solution vectors
    N3 = zeros(1, length(t));
    N2 = zeros(1, length(t));
    N1 = zeros(1, length(t));
    N0 = zeros(1, length(t));

    % Set initial values
    N3(1) = N_initial(1);
    N2(1) = N_initial(2);
    N1(1) = N_initial(3);
    N0(1) = N_initial(4);

    % Define the system of rate equations
    fN3 = @(N3, N2, N0) W * L * J / e - N3 / tau_3 - K * (N3 - N2) * N0;
    fN2 = @(N3, N2, N0) (1/tau_32 + 1/tau_sp) * N3 - N2 / tau_21 + K * (N3 - N2) * N0;
    fN1 = @(N3, N2, N1) N3 / tau_31 + N2 / tau_21 - N1 / tau_out;
    fN0 = @(N3, N2, N0) Num_of_state * K * (N3 - N2) * N0 + Num_of_state * beta * N3 / tau_sp - N0 / tau_p;

    % Runge-Kutta 4th Order Method
    for i = 1:(length(t) - 1)
        % k1 values
        k1_N3 = fN3(N3(i), N2(i), N0(i));
        k1_N2 = fN2(N3(i), N2(i), N0(i));
        k1_N1 = fN1(N3(i), N2(i), N1(i));
        k1_N0 = fN0(N3(i), N2(i), N0(i));

        % k2 values
        k2_N3 = fN3(N3(i) + 0.5 * h * k1_N3, N2(i) + 0.5 * h * k1_N2, N0(i) + 0.5 * h * k1_N0);
        k2_N2 = fN2(N3(i) + 0.5 * h * k1_N3, N2(i) + 0.5 * h * k1_N2, N0(i) + 0.5 * h * k1_N0);
        k2_N1 = fN1(N3(i) + 0.5 * h * k1_N3, N2(i) + 0.5 * h * k1_N2, N1(i) + 0.5 * h * k1_N1);
        k2_N0 = fN0(N3(i) + 0.5 * h * k1_N3, N2(i) + 0.5 * h * k1_N2, N0(i) + 0.5 * h * k1_N0);

        % k3 values
        k3_N3 = fN3(N3(i) + 0.5 * h * k2_N3, N2(i) + 0.5 * h * k2_N2, N0(i) + 0.5 * h * k2_N0);
        k3_N2 = fN2(N3(i) + 0.5 * h * k2_N3, N2(i) + 0.5 * h * k2_N2, N0(i) + 0.5 * h * k2_N0);
        k3_N1 = fN1(N3(i) + 0.5 * h * k2_N3, N2(i) + 0.5 * h * k2_N2, N1(i) + 0.5 * h * k2_N1);
        k3_N0 = fN0(N3(i) + 0.5 * h * k2_N3, N2(i) + 0.5 * h * k2_N2, N0(i) + 0.5 * h * k2_N0);

        % k4 values
        k4_N3 = fN3(N3(i) + h * k3_N3, N2(i) + h * k3_N2, N0(i) + h * k3_N0);
        k4_N2 = fN2(N3(i) + h * k3_N3, N2(i) + h * k3_N2, N0(i) + h * k3_N0);
        k4_N1 = fN1(N3(i) + h * k3_N3, N2(i) + h * k3_N2, N1(i) + h * k3_N1);
        k4_N0 = fN0(N3(i) + h * k3_N3, N2(i) + h * k3_N2, N0(i) + h * k3_N0);

        % Update values for the next step
        N3(i + 1) = N3(i) + (h / 6) * (k1_N3 + 2 * k2_N3 + 2 * k3_N3 + k4_N3);
        N2(i + 1) = N2(i) + (h / 6) * (k1_N2 + 2 * k2_N2 + 2 * k3_N2 + k4_N2);
        N1(i + 1) = N1(i) + (h / 6) * (k1_N1 + 2 * k2_N1 + 2 * k3_N1 + k4_N1);
        N0(i + 1) = N0(i) + (h / 6) * (k1_N0 + 2 * k2_N0 + 2 * k3_N0 + k4_N0);
    end

    % Steady-state values
    N_steady = [N3(end), N2(end), N1(end), N0(end)];

    % disp('N3_steady and Nph_steady in hamadous:');
    % disp(N_steady(1));
    % disp(N_steady(4));
    % disp('48 times fitted carrier numbers in hamadous')
    % disp(48*N_steady(1));
end
