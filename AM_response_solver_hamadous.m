%% Solve small signal modulation with varying steady values and tau_21
function [N_ph] = AM_response_solver_hamadous(freq, tau_21, N_steady_results)

global W L d Num_of_state tau_32 tau_31 tau_sp beta tau_p tau_out Gamma neff e c0 sigma_32

    % check value of steady states
    % -------------------------------
    % monitor
    % disp('N_steady:');
    % disp(N_steady_results);  % Print the steady-state values to check

    Delta_I = 1e-3; %1mA
    Delta_J = Delta_I/(W*L); % density
    
    tau_3 = 1 / (1/tau_32 + 1/tau_31 + 1/tau_sp); % lifetime of upper level (should be 1.4ps)
    
    c_node = c0 / neff;  % Reduced speed of light in the medium
    V = W * L * d * Num_of_state;  % Cavity volume (m^3)

    K = Gamma * c_node * sigma_32 / V; % Coupling factor

    % Extract the steady-state carrier and photon numbers from N_steady_results
    N_30 = double(N_steady_results(1));
    N_20 = double(N_steady_results(2));
    N_ph0 = double(N_steady_results(4));

    % Frequency range from 0.1 GHz to 30 GHz with 128 points
    A_sol = zeros(4, length(freq));    % Preallocate matrix for solutions

    disp('current tau21 in solving AM response')
    disp(tau_21*1e12)
    
    % Loop over each frequency to compute the response
    for idx = 1:length(freq)
        w = 2 * pi * freq(idx);        % Angular frequency
        
        % jw term with imaginary part
        jw = 1j * w;

        % Define matrix H(omega)
        H = [ ...
            jw + 1/tau_3 + K * N_ph0, -K * N_ph0, 1e-6, +K * (N_30 - N_20); ...
            -(1/tau_32 + 1/tau_sp) - K * N_ph0, jw + 1/tau_21 + K * N_ph0, 1e-6, -K * (N_30 - N_20); ...
            -1/tau_31, -1/tau_21, jw + 1/tau_out, 1e-6; ...
            -Num_of_state * K * N_ph0 - Num_of_state * beta/tau_sp , Num_of_state * K * N_ph0, 0, jw + 1/tau_p - Num_of_state * K * (N_30 - N_20)
        ];

        % Replace 0 to 1e-6

        % Define vector B for the current frequency
        B = [ ...
            W * L * Delta_J / e; ...
            0; ...
            0; ...
            0
        ];
        
        % Solve for A (A = H \ B)
        A = H \ B;  % More efficient than inv(H) * B

        % Ensure A is numeric (convert from symbolic if necessary)
        % A = double(A);  % Explicitly convert to numeric if needed

        % Store the solution A for the current frequency
        A_sol(:, idx) = A;
    end

    % Output the photon number response (N_ph)
    N_ph = A_sol(4, :);  % N_ph is the 4th element of the solution

end
