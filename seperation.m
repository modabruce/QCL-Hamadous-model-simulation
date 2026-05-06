close all
clear
clc
clf

%% set values of xi and J
% r: ratio of J/Jth
r_values = linspace(1.25, 1.25, 1);  % reproduce results in LaserMatrix
% xi: ratio of tau21/tau32
xi = linspace(0.05, 0.9, 18);
%xi = [0.1, 0.15, 0.25, 0.5, 0.6, 0.8, 0.9];
%xi = [0.8667,0.8857,0.9048];
%xi = 1/7;
Tau32 = 2.1e-12; 
tau_21_values = Tau32*xi; 
freq = linspace(0.1e9, 1000e9, 1280*6); % Frequency in Hz 

% time step and total time for runge_kutta method
time_step = 5.0e-15;
t_final = 2000e-12;
% freq = linspace(0.1e9, 3000e9, 1280);
% freq = 1j*freq;

% Define global parameters
global W L d Num_of_state tau_32 tau_31 tau_sp beta tau_p tau_out Gamma neff e c0 ng vg sigma_32
W = 34e-6; % Width (m)
L = 1e-3;  % Length (m)
d = 45e-9; % Thickness (m)
Num_of_state = 48;  % Number of stages

tau_32 = 2.1e-12;  % Nonradiative lifetime (s)
tau_31 = 4.2e-12;  % Nonradiative lifetime (s)
tau_sp = 38e-9;    % Spontaneous lifetime (s)
beta = 2e-3;  % Proportion of spontaneous emission events
tau_p = 3.36e-12;  % Photon lifetime

% Electron escape time
tau_out = 1e-12;  % (s)

Gamma = 0.32;      % Confinement factor

c0 = 2.9978e8;     % Exact speed of light (m/s)
neff = 3.27;       % Effective refractive index
ng = neff;         % Group refractive index
vg = c0 / ng;      % Group velocity

sigma_32 = 1.8e-18; % Stimulated emission cross-section (m^2)
e = 1.6e-19;        % Charge of an electron (C)

%% Save figures
% %folder = 'Your folder';
% filename_J = 'your_name';
% fullFilePath = fullfile(folder, filename_J);
% exportgraphics(gcf, fullFilePath, 'Resolution', 1000); % 300 DPI resolution

%% get f_parameter 1/2pi*tau_21
% f_parameter = 1./(2*pi*tau_21_values);
% f_parameter = f_parameter/1e9;
% disp('f_paramter in Ghz')
% disp(f_parameter)

%% get LM parameter epsilon and dgdN("a" in LaserMatrix guide) from tau_21 array
% in LaserMatrix, epsilon as par 1 adn dgdN as par 2
[epsilon, dgdN] = LM_parameter(tau_21_values);

%% get value of Jth with varing tau_21
Jth_value = zeros(1,length(tau_21_values));

for j = 1:length(tau_21_values)
    Jth_value(j) = get_Jth(tau_21_values(j));
end

% Preallocate arrays to store results for each r
N3_results_Hamadous = cell(length(tau_21_values), length(r_values));
N2_results_Hamadous = cell(length(tau_21_values), length(r_values));
N1_results_Hamadous = cell(length(tau_21_values), length(r_values));
N0_results_Hamadous = cell(length(tau_21_values), length(r_values));
N_steady_Hamadous = cell(length(tau_21_values), length(r_values));  

Q_results_SL = cell(length(tau_21_values), length(r_values));
S_results_SL = cell(length(tau_21_values), length(r_values));
N_steady_SL = cell(length(tau_21_values), length(r_values));   % Store steady-state values

%% solve Hamadous runge_kutta and SL runge_kutta
% Loop over r values and run the solver
for j = 1:length(tau_21_values)
    tau_21 = tau_21_values(j);
    Jth = Jth_value(j);

for i = 1:length(r_values)
    r = r_values(i);  % Current value of r
    
    % Call the runge_kutta_solver for each r of hamadous 4 rate equations
    [N3, N2, N1, N0, N_steady_HD] = runge_kutta_solver_hamadous(r,tau_21, Jth, time_step, t_final);
    
    % Store the results
    N3_results_Hamadous{j,i} = N3;
    N2_results_Hamadous{j,i} = N2;
    N1_results_Hamadous{j,i} = N1;
    N0_results_Hamadous{j,i} = N0;
    N_steady_Hamadous{j,i} = N_steady_HD;  % Store steady-state values for small signal modulation


    % Call the runge_kutta_solver for each r of SL 2 rate equations
    [carrier_num, photon_num, N_steady_2_level] = runge_kutta_solver_SL(r, tau_21, Jth, time_step, t_final);

    Q_results_SL{j,i} = carrier_num;
    S_results_SL{j,i} = photon_num;
    N_steady_SL{j,i} = N_steady_2_level;  % Store steady-state values for small signal modulation
end

end

%% test steady value accuracy 
N3_test = zeros(1,length(tau_21_values));
Nph_test = zeros(1,length(tau_21_values));
SL_Q_test = zeros(1,length(tau_21_values));
SL_S_test = zeros(1,length(tau_21_values));

for i = 1:length(tau_21_values)
    N3_test(1,i) = 48*N_steady_Hamadous{j,1}(1);  % 48 times fitted MAT to LM
    Nph_test(1,i) = N_steady_Hamadous{j,1}(4);
    SL_Q_test(1,i) = N_steady_SL{j,1}(1);
    SL_S_test(1,i) = N_steady_SL{j,1}(2);
end

Carrier_relative_error = abs(N3_test - SL_Q_test)./N3_test;
Carrier_relative_error = Carrier_relative_error*100; % 100%
Photon_relative_error = abs(Nph_test - SL_S_test)./Nph_test;
Photon_relative_error = 100*Photon_relative_error; % 100%

%test value
disp('tau_values:');
disp(tau_21_values);
disp('Carrier of hamadous:');
disp(N3_test);
disp('Carrier of SL:');
disp(SL_Q_test);
disp('photon of hamadous:');
disp(Nph_test);
disp('photon of SL:');
disp(SL_S_test);
disp('difference of carrier')
disp(Carrier_relative_error)
disp('difference of photon')
disp(Photon_relative_error)

%% Solve determinant to get omega array
% omega array can be shown from eterminant_solver function
omega_array_Hamadous = cell(length(tau_21_values), length(r_values));
omega_array_SL = cell(length(tau_21_values), length(r_values));

pole_frequencies_hamadous = cell(length(tau_21_values), length(r_values));
zero_frequencies_hamadous = cell(length(tau_21_values), length(r_values));
T_magnitude_hamadous = cell(length(tau_21_values), length(r_values));
am_hamadous = cell(length(tau_21_values), length(r_values));   %P  amplitude
R_hamadous = cell(length(tau_21_values), length(r_values));  %R  resonance poles
K_hamadous = cell(length(tau_21_values), length(r_values));  %K

pole_frequencies_SL = cell(length(tau_21_values), length(r_values));
zero_frequencies_SL = cell(length(tau_21_values), length(r_values));
T_magnitude_SL = cell(length(tau_21_values), length(r_values));
P_SL = cell(length(tau_21_values), length(r_values));   %P  amplitude
R_SL = cell(length(tau_21_values), length(r_values));  %R  resonance poles
K_SL = cell(length(tau_21_values), length(r_values));  %K

for j = 1:length(tau_21_values)
    for i = 1:length(r_values)
        % Get the steady-state values for current parameters
        N_steady_Hamadous_current = N_steady_Hamadous{j,i};
        N_steady_SL_current = N_steady_SL{j,i};
        
        % Solve determinants for omega arrays
        omega_array_Hamadous{j,i} = determinant_solver_hamadous(tau_21_values(j), N_steady_Hamadous_current);
        omega_array_SL{j,i} = determinant_solver_SL(tau_21_values(j), N_steady_SL_current);
        
        % Call the zero_pole_residue function with appropriate parameters
        [pf_hamadous, zf_hamadous, Tm_hamadous, res_R_hamadous, res_P_hamadous, res_K_hamadous] = HDM_zero_pole_residue(freq, tau_21_values(j), N_steady_Hamadous_current);
        [pf_SL, zf_SL, Tm_SL, res_R_SL, res_P_SL, res_K_SL] = SL_zero_pole_residue(freq, tau_21_values(j), N_steady_SL_current);
        
        % Store the results in the corresponding cells
        pole_frequencies_hamadous{j,i} = pf_hamadous;    pole_frequencies_SL{j,i} = pf_SL;
        zero_frequencies_hamadous{j,i} = zf_hamadous;    zero_frequencies_SL{j,i} = zf_SL;
        T_magnitude_hamadous{j,i} = Tm_hamadous;         T_magnitude_SL{j,i} = Tm_SL;
        R_hamadous{j,i} = res_R_hamadous;                R_SL{j,i} = res_R_SL;
        am_hamadous{j,i} = res_P_hamadous;                P_SL{j,i} = res_P_SL;
        K_hamadous{j,i} = res_K_hamadous;                K_SL{j,i} = res_K_SL;
    end
end
%% Get resonance frequency && get relax frequency (omega_1 and omega_2)
pole_zero_save = zeros(length(tau_21_values),length(r_values));
omega1_hamadous = pole_zero_save;
omega2_hamadous = pole_zero_save;
omega1_SL = pole_zero_save;
omega2_SL = pole_zero_save;
zeros_hamadous = pole_zero_save;
zeros_SL = pole_zero_save;

for j = 1:length(tau_21_values)
for i = 1:length(r_values)
    omega1_hamadous(j,i) = omega_array_Hamadous{j,i}(1);
    omega2_hamadous(j,i) = omega_array_Hamadous{j,i}(2);
    omega1_SL(j,i) = omega_array_SL{j,i}(1);
    omega2_SL(j,i) = omega_array_SL{j,i}(2);
    zeros_hamadous(j,i) = zero_frequencies_hamadous{j,i}(1);

    % Check if zero_frequencies_SL{j,i} is empty
        if isempty(zero_frequencies_SL{j, i})
            zeros_SL(j, i) = NaN; % Assign NaN if empty
        else
            % Assign the first zero frequency (or any other logic if multiple zeros exist)
            zeros_SL(j, i) = zero_frequencies_SL{j, i}(1); 
        end
end
end

omega1_hamadous = abs(omega1_hamadous);
omega2_hamadous = abs(omega2_hamadous);
omega1_SL = abs(omega1_SL);
omega2_SL = abs(omega2_SL);
zeros_hamadous = abs(zeros_hamadous);
zeros_SL = abs(zeros_SL);

% Test reson_freq and relax_freq
% figure
% plot(tau_21_values*1e12,abs(omega1_MAT))  % abs or imag
% %plot(tau_21_values*1e12,abs(relax_freq)/1e9)

%% Solve AM-response after normalization
% N_steady can be shown from AM_response_solver fucntion
number_of_photons_hamadous = cell(length(tau_21_values), length(r_values));
number_of_photons_SL = cell(length(tau_21_values), length(r_values));
Nor_number_of_photons_hamadous = cell(length(tau_21_values), length(r_values));
Nor_number_of_photons_SL{j,i} = cell(length(tau_21_values), length(r_values));

for j = 1:length(tau_21_values)
    tau_21 = tau_21_values(j);

for i = 1:length(r_values)
    number_of_photons_hamadous{j,i} = AM_response_solver_hamadous(freq, tau_21, N_steady_Hamadous{j,i});
    Nor_number_of_photons_hamadous{j,i} = number_of_photons_hamadous{j,i}./number_of_photons_hamadous{j,i}(1);
    number_of_photons_SL{j,i} = 48*AM_response_solver_SL(freq, tau_21, N_steady_SL{j,i});
    Nor_number_of_photons_SL{j,i} = number_of_photons_SL{j,i}./number_of_photons_SL{j,i}(1);
end

end
%% compare LM and MATLAB first pole in frequency
figure(1)

plot(xi, omega1_hamadous, 'o-', 'MarkerFaceColor', 'blue', 'DisplayName', 'Hamadous'); % Circle markers
hold on
plot(xi, omega1_SL, 's--', 'MarkerFaceColor', 'red', 'DisplayName', 'SL'); % Square markers
xlabel('\xi');
ylabel('First pole frequency (Ghz)');
title('\omega_{1} for different \xi Values');
legend show; 
grid on;
hold off;

% Adjust figure size 8x5 inches)
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 3]); % [left, bottom, width, height]
% Set renderer to ensure quality
set(gcf, 'Renderer', 'painters');
% 
% folder = 'C:\Users\bruce\Desktop\QCL figures';
% filename = 'Pole vs zero.png';
% fullFilePath = fullfile(folder, filename);
% exportgraphics(gcf, fullFilePath, 'Resolution', 1200); % 1200 DPI resolution
%% plot hamadous pole and zero in same figure for r=1.25, tau_21 from 1.8 to 2.0

figure(2)

plot(xi, omega2_hamadous, 'o-', 'MarkerFaceColor', 'blue', 'DisplayName', 'Hamadous'); % Circle markers
hold on
plot(xi, omega2_SL, 's--', 'MarkerFaceColor', 'red', 'DisplayName', 'SL'); % Square markers
xlabel('\xi');
ylabel('Second pole frequency (Ghz)');
title('\omega_{2} for different \xi Values');
legend ('show','FontSize', 12) 
grid on;
hold off;

% Adjust figure size 8x5 inches)
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 3]); % [left, bottom, width, height]
% Set renderer to ensure quality
set(gcf, 'Renderer', 'painters');



%% log and linear error 
log_mean_squared_difference = zeros(1,length(tau_21_values));
lin_mean_squared_difference = zeros(1,length(tau_21_values));
log_meandB = cell(1,length(tau_21_values));
lin_meandB = cell(1,length(tau_21_values));
 
for j = 1:length(tau_21_values)
    tau_21 = tau_21_values(j);
    
    % Extract the number_of_photons for the current tau_21 (for i=1 since r_values has only one value)
    photon_response_hamadous = abs(number_of_photons_hamadous{j,1});
    photon_response_SL = abs(number_of_photons_SL{j,1});
    first_value_hamadous = photon_response_hamadous(1);

    % option 1 log error
    log_squared_difference = abs((20*log10(photon_response_hamadous)-20*log10(photon_response_SL)));
    log_squared_diff_plot = log_squared_difference ;
    % Average the squared difference to get the desired value
    log_mean_squared_difference(j) = mean(log_squared_difference);
    log_meandB{1,j} = log_squared_difference;

    % option 2 linear error
    lin_squared_difference = abs((photon_response_hamadous) - (photon_response_SL))./abs(first_value_hamadous);
    lin_squared_diff_plot = 20*log10(lin_squared_difference);
    % Average the squared difference to get the desired value
    lin_mean_squared_difference(j) = mean(lin_squared_difference);
    lin_meandB{1,j} = 10*log10(lin_squared_difference);

end

%%
% Plot for log_mean_squared_difference
figure(3)
set(gca, 'XScale', 'log');
plot(xi, (log_mean_squared_difference), 'o-', 'DisplayName', 'LogError mean squared difference', 'LineWidth', 2);
xlabel('\xi', 'FontSize', 12);
ylabel('Average error (dB)', 'FontSize', 12);
title('Log Mean Squared Error', 'FontSize', 12);
grid on;
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 3]); % [left, bottom, width, height]
set(gcf, 'Renderer', 'painters');

% Plot for lin_mean_squared_difference
figure(4)
set(gca, 'XScale', 'log');
plot(xi, (lin_mean_squared_difference), 'o-', 'DisplayName', 'LinError mean squared difference', 'LineWidth', 2);
xlabel('\xi', 'FontSize', 12);
ylabel('Average error (dB)', 'FontSize', 12);
title('Lin Mean Squared Error', 'FontSize', 12);
grid on;
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 3]); % [left, bottom, width, height]
set(gcf, 'Renderer', 'painters');









