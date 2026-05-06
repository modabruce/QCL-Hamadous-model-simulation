close all
clf
clear all
clc

%% set values of xi and J
% r: ratio of J/Jth
r_values = [1.25,1.5,2,2.5];  % reproduce results in LaserMatrix
% xi: ratio of tau21/tau32
xi = linspace(1/7, 1/7, 1);
Tau32 = 2.1e-12;
tau_21_values = Tau32*xi; 
freq = linspace(0.1e9, 300e9, 1280*3); % Frequency in Hz 

% time step and total time for runge_kutta method
time_step = 5.0e-15;
t_final = 300e-12;
% freq = linspace(0.1e9, 3000e9, 1280);
% freq = 1j*freq;

t0 = 0.1e-12;
t = t0:time_step:t_final;

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



%% 
figure(1)
hold on; % Keep the plots on the same figure
grid on

% Define a set of colors (you can use a predefined colormap or specify them manually)
colors = lines(length(r_values)); % 'lines' colormap will generate unique colors for each iteration

% Loop through the values of r
for i = 1:length(r_values)
    r_current = r_values(i);

    % Get the color for the current iteration
    color = colors(i, :);  % Fetch the i-th color for the loop

    % Plot photon_response_hamadous with solid lines and the selected color
    plot(t*1e12, ...
         N3_results_Hamadous{j,i}, ...
         'LineWidth', 2, ...
         'LineStyle', '-', ...
         'Color', color, ... % Set color here
         'DisplayName',['$Hamadous:\ J/J_{th} = ', num2str(r_current), '$']);

    % Plot Q_results_SL with dashed lines and the same color
    plot(t*1e12, ...
         Q_results_SL{j,i}/48, ...
         'LineWidth', 2, ...
         'LineStyle', '--', ...
         'Color', color, ... % Set color here too
         'DisplayName',['$SL:\ J/J_{th} = ', num2str(r_current), '$']);
end

% Customize the plot
legend('Interpreter', 'latex');
xlabel('Time(ps)', 'FontSize', 12);
ylabel('Number of carriers', 'FontSize', 12);
legend('show', 'Location', 'southwest', 'FontSize', 8); 
hold off;

% Adjust figure size (8x5 inches)
set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 5]); % [left, bottom, width, height]
% Set renderer to ensure quality
set(gcf, 'Renderer', 'painters');



%% 
figure(2)
hold on; % Keep the plots on the same figure
grid on

% Define a set of colors (you can use a predefined colormap or specify them manually)
colors = lines(length(r_values)); % 'lines' colormap will generate unique colors for each iteration

% Loop through the values of r
for i = 1:length(r_values)
    r_current = r_values(i);

    % Get the color for the current iteration
    color = colors(i, :);  % Fetch the i-th color for the loop

    % Plot photon_response_hamadous with solid lines and the selected color
    plot(t*1e12, ...
         N0_results_Hamadous{j,i}, ...
         'LineWidth', 2, ...
         'LineStyle', '-', ...
         'Color', color, ... % Set color here
         'DisplayName',['$Hamadous:\ J/J_{th} = ', num2str(r_current), '$']);

    % Plot S_results_SL with dashed lines and the same color
    plot(t*1e12, ...
         S_results_SL{j,i}, ...
         'LineWidth', 2, ...
         'LineStyle', '--', ...
         'Color', color, ... % Set color here too
         'DisplayName',['$SL:\ J/J_{th} = ', num2str(r_current), '$']);
end

% Customize the plot
legend('Interpreter', 'latex');
xlabel('Time(ps)', 'FontSize', 12);
ylabel('Number of photons', 'FontSize', 12);
legend('show', 'Location', 'southwest', 'FontSize', 5); 
hold off;

% Adjust figure size (8x5 inches)
set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 5]); % [left, bottom, width, height]
% Set renderer to ensure quality
set(gcf, 'Renderer', 'painters');
