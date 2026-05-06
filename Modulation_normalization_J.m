close all
clear
clc
clf

%% set values of xi and J
% r: ratio of J/Jth
r_values = linspace(1.25, 10, 4);  % reproduce results in LaserMatrix
% xi: ratio of tau21/tau32
xi = linspace(1/7, 1/7, 1);
Tau32 = 2.1e-12;
tau_21_values = Tau32*xi; 
freq = linspace(0.1e9, 1000e9, 1280*8); % Frequency in Hz 

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
N3_test = zeros(1,length(r_values));
Nph_test = zeros(1,length(r_values));
SL_Q_test = zeros(1,length(r_values));
SL_S_test = zeros(1,length(r_values));

for i = 1:length(r_values)
    N3_test(1,i) = 48*N_steady_Hamadous{1,i}(1);  % 48 times fitted MAT to LM
    Nph_test(1,i) = N_steady_Hamadous{1,i}(4);
    SL_Q_test(1,i) = N_steady_SL{1,i}(1);
    SL_S_test(1,i) = N_steady_SL{1,i}(2);
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
disp('difference of carrier in %')
disp(Carrier_relative_error)
disp('difference of photon in %')
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

plot(r_values, omega1_hamadous, ':', 'Color', 'blue', 'LineWidth', 2, 'MarkerFaceColor', 'blue', 'DisplayName', 'Hamadous'); % Blue dotted line with circle markers
hold on
plot(r_values, omega1_SL, 'o', 'Color', 'red', 'MarkerFaceColor', 'red', 'LineWidth', 1, 'DisplayName', 'SL'); % Red line with circle markers

xlabel('$J/J_{th}$', 'Interpreter', 'latex');
ylabel('First pole frequency (Ghz)');
title('\omega_{1} for different bias ratio');
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

plot(r_values, omega2_hamadous, 'o-', 'MarkerFaceColor', 'blue', 'DisplayName', 'Hamadous'); % Circle markers
hold on
plot(r_values, omega2_SL, 's--', 'MarkerFaceColor', 'red', 'DisplayName', 'Hamadous'); % Square markers
xlabel('$J/J_{th}$', 'Interpreter', 'latex');
ylabel('Frequency (Ghz)');
title('\omega_{2} frequency for different bias ratio');
legend ('show','FontSize', 12) 
grid on;
hold off;

% Adjust figure size 8x5 inches)
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 3]); % [left, bottom, width, height]
% Set renderer to ensure quality
set(gcf, 'Renderer', 'painters');

% folder = 'C:\Users\bruce\Desktop\QCL figures';
% filename = 'Pole vs zero.png';
% fullFilePath = fullfile(folder, filename);
% exportgraphics(gcf, fullFilePath, 'Resolution', 1200); % 1200 DPI resolution

%% plot zeros in same figure
figure(3)
plot(r_values, zeros_hamadous/(2 * pi * 1e9), 'o-', 'MarkerFaceColor', 'blue', 'DisplayName', 'Hamadous'); % Circle markers
hold on
plot(r_values, zeros_SL/(2 * pi * 1e9), 's--', 'MarkerFaceColor', 'red', 'DisplayName', 'SL'); % Square markers
xlabel('$J/J_{th}$', 'Interpreter', 'latex');ylabel('Zero point frequency');
title('Zero point for different J/J_{th} Values');
legend show; 
grid on;
hold off;


%% plot AM response in same figure(photon number)
figure(4)
hold on; % Keep the plots on the same figure
set(gca, 'XScale', 'log');
grid on

% Define a color palette to distinguish lines for each tau_21 value
color_palette = lines(length(r_values)); % Generates distinct colors

% Loop through the values of r
for i = 1:length(r_values)
    r_current = r_values(i);
    % Extract the number_of_photons for the current tau_21 (for j=1 since xi has only one value)
    photon_response_hamadous = number_of_photons_hamadous{1,i};
    photon_response_SL = number_of_photons_SL{1,i};
    % Define line styles and markers
    line_color = color_palette(i, :); % Assign a unique color for each tau_21

    % Plot photon_response_hamadous with solid lines and markers
    plot(abs(freq) / 1e9, ...
         20 * log10(abs(photon_response_hamadous)), ...
         'LineWidth', 2, ...
         'LineStyle', '-', ...
         'MarkerIndices', 1:10:length(freq), ... % Add markers at intervals
         'Color', line_color, ...
         'DisplayName',['$Hamadous:\ J/J_{th} = ', num2str(r_current), '$']);

    % Plot photon_response_SL with dashed lines and markers
    plot(abs(freq) / 1e9, ...
         20 * log10(abs(photon_response_SL)), ...
         'LineWidth', 2, ...
         'LineStyle', '--', ...
         'MarkerIndices', 1:10:length(freq), ... % Add markers at intervals
         'Color', line_color, ...
         'DisplayName',['$SL:\ J/J_{th} = ', num2str(r_current), '$']);

end

% Customize the plot
legend('Interpreter', 'latex');
xlabel('Modulation Frequency (GHz)', 'FontSize', 10);
ylabel('Photon response (dBe)', 'FontSize', 10);
%title('Photon Response Comparison', 'FontSize', 14);
legend('show', 'Location', 'southwest', 'FontSize', 6); 
hold off;

% Adjust figure size 8x5 inches)
set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 5]); % [left, bottom, width, height]
% Set renderer to ensure quality
set(gcf, 'Renderer', 'painters');

folder = 'C:\Users\bruce\Desktop\QCL figures';
filename = 'HDM_SL_AM.png';
fullFilePath = fullfile(folder, filename);
exportgraphics(gcf, fullFilePath, 'Resolution', 1200); % 1200 DPI resolution

%% plot AM response in same figure(normalization)
figure(5)
hold on; % Keep the plots on the same figure
set(gca, 'XScale', 'log');
grid on

% Define a color palette to distinguish lines for each tau_21 value
color_palette = lines(length(r_values)); % Generates distinct colors

% Loop through the values of r
for i = 1:length(r_values)
    r_current = r_values(i);

    % Extract the number_of_photons for the current tau_21 (for j=1 since xi has only one value)
    photon_response_hamadous = Nor_number_of_photons_hamadous{1,i};
    photon_response_SL = Nor_number_of_photons_SL{1,i};
    % Define line styles and markers
    line_color = color_palette(i, :); % Assign a unique color for each tau_21

    % Plot photon_response_hamadous with solid lines and markers
    plot(abs(freq) / 1e9, ...
         20 * log10(abs(photon_response_hamadous)), ...
         'LineWidth', 2, ...
         'LineStyle', '-', ...
         'MarkerIndices', 1:10:length(freq), ... % Add markers at intervals
         'Color', line_color, ...
         'DisplayName', ['$Hamadous:\ J/J_{th} = ', num2str(r_current), '$']);

    % Plot photon_response_SL with dashed lines and markers
    plot(abs(freq) / 1e9, ...
         20 * log10(abs(photon_response_SL)), ...
         'LineWidth', 2, ...
         'LineStyle', '--', ...
         'MarkerIndices', 1:10:length(freq), ... % Add markers at intervals
         'Color', line_color, ...
         'DisplayName', ['$SL:\ J/J_{th} = ', num2str(r_current), '$']);

end

% Customize the plot
legend('Interpreter', 'latex');
xlabel('Modulation Frequency (GHz)', 'FontSize', 10);
ylabel('Photon response (dBe)', 'FontSize', 10);
%title('Photon Response Comparison', 'FontSize', 14);
legend('show', 'Location', 'southwest', 'FontSize', 7.5); 
hold off;

% Adjust figure size 8x5 inches)
set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 3]); % [left, bottom, width, height]
% Set renderer to ensure quality
set(gcf, 'Renderer', 'painters');

folder = 'C:\Users\bruce\Desktop\QCL figures';
filename = 'HDM_SL_AM_N.png';
fullFilePath = fullfile(folder, filename);
exportgraphics(gcf, fullFilePath, 'Resolution', 1200); % 1200 DPI resolution

%% log and linear error 
log_epsilon = zeros(1,length(r_values));
lin_epsilon = zeros(1,length(r_values));
logE_freq = cell(1,length(r_values));
linE_freq = cell(1,length(r_values));
 
for i = 1:length(r_values)
    r_current = r_values(i);
    
    % Extract the number_of_photons for the current tau_21 (for i=1 since r_values has only one value)
    photon_response_hamadous = abs(number_of_photons_hamadous{1,i});
    photon_response_SL = abs(number_of_photons_SL{1,i});
    first_value_hamadous = photon_response_hamadous(1);

    % option 1 log error
    log_Ef = abs(2*log10(photon_response_hamadous)-2*log10(photon_response_SL));
    
    logE_freq{1,i} = log_Ef;

    % Average the squared difference to get the desired value
    log_epsilon(i) = 10* mean(log_Ef);
    

    % option 2 linear error
    lin_Ef = abs((photon_response_hamadous) - (photon_response_SL))./abs(first_value_hamadous);
    lin_Ef = lin_Ef.^2;
    linE_freq{1,i} = lin_Ef;
    % Average the squared difference to get the desired value
    lin_epsilon(i) = 10*log10(mean(lin_Ef));

end

%%
figure(6)
% Define the color palette
color_palette = lines(length(r_values));

% First subplot: logE_freq
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
hold on;
set(gca, 'XScale', 'log'); % Set x-axis to logarithmic scale
grid on;

for i = 1:length(r_values)
    r_current = r_values(i);
    line_color = color_palette(i, :);

    % Plot logE_freq
    plot(abs(freq)/1e9, logE_freq{i}, 'Color', line_color, ...
        'LineStyle', '-', 'DisplayName', ['$\mathrm{log\ error:\ } J/J_{\mathrm{th}} = ', num2str(r_current), '$'], ...
        'LineWidth', 2);
end

ylabel('Squared difference (dB)', 'FontSize', 12, 'Interpreter', 'latex'); % LaTeX y-label
title('Log Error Response', 'FontSize', 12, 'Interpreter', 'latex'); % LaTeX title
legend('show', 'Location', 'southwest', 'FontSize', 8, 'Interpreter', 'latex'); % Set legend interpreter to LaTeX
hold off;

% Second subplot: linE_freq
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
hold on;
set(gca, 'XScale', 'log'); % Set x-axis to logarithmic scale
grid on;

for i = 1:length(r_values)
    r_current = r_values(i);
    line_color = color_palette(i, :);

    % Plot linE_freq
    plot(abs(freq)/1e9, linE_freq{i}, 'Color', line_color, ...
        'LineStyle', '--', 'DisplayName', ['$\mathrm{lin\ error:\ } J/J_{\mathrm{th}} = ', num2str(r_current), '$'], ...
        'LineWidth', 2);
end

xlabel('Modulation frequency (GHz)', 'FontSize', 12, 'Interpreter', 'latex'); % LaTeX x-label
ylabel('Squared difference', 'FontSize', 12, 'Interpreter', 'latex'); % LaTeX y-label
title('Linear Error Response', 'FontSize', 12, 'Interpreter', 'latex'); % LaTeX title
legend('show', 'Location', 'southwest', 'FontSize', 8, 'Interpreter', 'latex'); % Set legend interpreter to LaTeX
hold off;

% Adjust figure size
set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 5]); % [left, bottom, width, height] with increased height for subplots
set(gcf, 'Renderer', 'painters');


%% Mean squared error

figure(7)

% Define the color palette for consistency
color_palette = lines(2); % Two colors for log_epsilon and lin_epsilon

% First subplot: log_epsilon
subplot(2, 1, 1); % 2 rows, 1 column, first subplot
hold on;
grid on;

% Plot log_epsilon
plot(r_values, log_epsilon, 'o-', 'Color', color_palette(1, :), ...
    'MarkerFaceColor', color_palette(1, :), 'LineWidth', 2, ...
    'DisplayName', '$\epsilon_{\log}$'); % Use LaTeX syntax for legend
legend('show', 'Location', 'southwest', 'FontSize', 12, 'Interpreter', 'latex'); % LaTeX legend
ylabel('Average error (dB)', 'FontSize', 12, 'Interpreter', 'latex'); % LaTeX y-label
title('Log Error Mean Squared Difference', 'FontSize', 12, 'Interpreter', 'latex'); % LaTeX title
hold off;

% Second subplot: lin_epsilon
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
hold on;
grid on;

% Plot lin_epsilon
plot(r_values, lin_epsilon, 'o-', 'Color', color_palette(2, :), ...
    'MarkerFaceColor', color_palette(2, :), 'LineWidth', 2, ...
    'DisplayName', '$\epsilon_{\mathrm{lin}}$'); % Use LaTeX syntax for legend
legend('show', 'Location', 'southwest', 'FontSize', 12, 'Interpreter', 'latex'); % LaTeX legend
xlabel('$J/J_{\mathrm{th}}$', 'FontSize', 12, 'Interpreter', 'latex'); % LaTeX x-label
ylabel('Average error (dB)', 'FontSize', 12, 'Interpreter', 'latex'); % LaTeX y-label
title('Linear Error Mean Squared Difference', 'FontSize', 12, 'Interpreter', 'latex'); % LaTeX title
hold off;

% Adjust the figure size
set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 5]); % [left, bottom, width, height] with increased height for subplots
set(gcf, 'Renderer', 'painters');



%% response with zeros and poles for both hamadous and SL with normalization
figure(9)
hold on; % Keep the plots on the same figure
set(gca, 'XScale', 'log');  % Set X-axis to logarithmic scale
grid on
color_palette = lines(length(r_values));
% Loop through the values of tau_21
for i = 1:length(r_values)
    r_current = r_values(i);
   
    % Extract the number_of_photons for the current tau_21 (for i=1 since r_values has only one value)
    photon_response_hamadous = Nor_number_of_photons_hamadous{1,i};
    photon_response_SL = Nor_number_of_photons_SL{1,i};
    line_color = color_palette(i, :); % Assign a unique color for each tau_21

    % Extract zero and pole frequencies for current tau_21
    zeros_freq_hm = zero_frequencies_hamadous{1,i}./(2*pi);  % zero frequency
    poles_freq_hm = pole_frequencies_hamadous{1,i}(1)./(2*pi);  % all pole freq
    poles_freq_SL = pole_frequencies_SL{1,i}./(2*pi);  % all pole freq

    zero_point_hamadous = AM_response_solver_hamadous(abs(zeros_freq_hm), tau_21_values, N_steady_Hamadous{1,i});
    zero_point_hamadous = zero_point_hamadous/number_of_photons_hamadous{1,i}(1);
    pole_point_hamadous = AM_response_solver_hamadous(abs(poles_freq_hm), tau_21_values, N_steady_Hamadous{1,i});
    pole_point_hamadous = pole_point_hamadous/number_of_photons_hamadous{j,1}(1);

    pole_point_SL = 48*AM_response_solver_SL(abs(poles_freq_SL), tau_21_values, N_steady_SL{1,i});
    pole_point_SL = pole_point_SL/number_of_photons_SL{1,i}(1);

   % Plot photon response with solid line
    plot(abs(freq) / 1e9, ...
         20 * log10(abs(photon_response_hamadous)), ...
         'LineWidth', 1.5, ...
         'Color', line_color, ...
         'LineStyle', '-', ...
         'DisplayName', ['Hamadous Response: $\ J/J_{th} = ', num2str(r_current), '$']);

     plot(abs(freq) / 1e9, ...
         20 * log10(abs(photon_response_SL)), ...
         'LineWidth', 1.5, ...
         'Color', line_color, ...
         'LineStyle', '--', ...
         'DisplayName', ['SL Response: $\ J/J_{th} = ', num2str(r_current), '$']);

    % Plot zeros as circles
    plot(abs(zeros_freq_hm) / 1e9, ...
         20 * log10(abs(zero_point_hamadous)), ...
         'o', ...
         'MarkerSize', 10, ...
         'LineWidth', 3, ...
         'MarkerFaceColor', 'none', ...
         'MarkerEdgeColor', line_color, ...
         'DisplayName', ['Hamadous zeros: $\ J/J_{th} = ', num2str(r_current), '$']);

end

hold off;

% Add labels, title, and grid
legend('Interpreter', 'latex');
xlabel('Modulation Frequency (GHz)', 'FontSize', 10);
ylabel('Photon response (dBe)', 'FontSize', 10);
%title('Photon Response with poles and zeros frequencies', 'FontSize', 12);
legend('show', 'Location', 'southwest', 'FontSize', 7);
grid on;

% Adjust figure size 8x5 inches)
set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 5]); % [left, bottom, width, height]
% Set renderer to ensure quality
set(gcf, 'Renderer', 'painters');

folder = 'C:\Users\bruce\Desktop\QCL figures';
filename = '2.png';
fullFilePath = fullfile(folder, filename);
exportgraphics(gcf, fullFilePath, 'Resolution', 1000); % 300 DPI resolution

%% response with zeros and poles for both hamadous and SL without normalization
figure(10)
hold on; % Keep the plots on the same figure
set(gca, 'XScale', 'log');  % Set X-axis to logarithmic scale
grid on
color_palette = lines(length(r_values));
% Loop through the values of tau_21
for i = 1:length(r_values)
    r_current = r_values(i);
   
    % Extract the number_of_photons for the current tau_21 (for i=1 since r_values has only one value)
    photon_response_hamadous = number_of_photons_hamadous{1,i};
    photon_response_SL = number_of_photons_SL{1,i};
    line_color = color_palette(i, :); % Assign a unique color for each tau_21

    % Extract zero and pole frequencies for current tau_21
    zeros_freq_hm = zero_frequencies_hamadous{1,i}./(2*pi);  % zero frequency
    poles_freq_hm = pole_frequencies_hamadous{1,i}(1)./(2*pi);  % all pole freq
    poles_freq_SL = pole_frequencies_SL{1,i}./(2*pi);  % all pole freq

    zero_point_hamadous = AM_response_solver_hamadous(abs(zeros_freq_hm), tau_21_values, N_steady_Hamadous{1,i});
    pole_point_hamadous = AM_response_solver_hamadous(abs(poles_freq_hm), tau_21_values, N_steady_Hamadous{1,i});
    pole_point_SL = 48*AM_response_solver_SL(abs(poles_freq_SL), tau_21_values, N_steady_SL{1,i});
    

   % Plot photon response with solid line
    plot(abs(freq) / 1e9, ...
         20 * log10(abs(photon_response_hamadous)), ...
         'LineWidth', 1.5, ...
         'Color', line_color, ...
         'LineStyle', '-', ...
         'DisplayName', ['Hamadous Response: $\ J/J_{th} = ', num2str(r_current), '$']);

     plot(abs(freq) / 1e9, ...
         20 * log10(abs(photon_response_SL)), ...
         'LineWidth', 1.5, ...
         'Color', line_color, ...
         'LineStyle', '--', ...
         'DisplayName', ['SL Response: $\ J/J_{th} = ', num2str(r_current), '$']);

    % Plot zeros as circles
    plot(abs(zeros_freq_hm) / 1e9, ...
         20 * log10(abs(zero_point_hamadous)), ...
         'o', ...
         'MarkerSize', 10, ...
         'LineWidth', 3, ...
         'MarkerFaceColor', 'none', ...
         'MarkerEdgeColor', line_color, ...
         'DisplayName', ['Hamadous zeros: $\ J/J_{th} = ', num2str(r_current), '$']);

end

hold off;

% Add labels, title, and grid
legend('Interpreter', 'latex');
xlabel('Modulation Frequency (GHz)', 'FontSize', 10);
ylabel('Photon response (dBe)', 'FontSize', 10);
%title('Photon Response with poles and zeros frequencies', 'FontSize', 12);
legend('show', 'Location', 'southwest', 'FontSize', 7);
grid on;

% Adjust figure size 8x5 inches)
set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 5]); % [left, bottom, width, height]
% Set renderer to ensure quality
set(gcf, 'Renderer', 'painters');

folder = 'C:\Users\bruce\Desktop\QCL figures';
filename = '2.png';
fullFilePath = fullfile(folder, filename);
exportgraphics(gcf, fullFilePath, 'Resolution', 1000); % 300 DPI resolution