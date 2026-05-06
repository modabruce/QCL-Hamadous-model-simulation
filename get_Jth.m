%% Get Jth with varying tau_21
function [Jth] = get_Jth(tau_21)

W = 34e-6; % lateral widths  34 micrometer = 34e-6m 
L = 1e-3; % lateral length 1mm=1e-3m
d = 45e-9; % 45nm = 45e-9m

neff = 3.27;
c0 = 2.9978e8;
ng= neff;
vg = c0/ng;  % group velocity
 
N = 48;   % number of states
q = 1.602e-19;
V = W*L*d*N; % m^3

sigma_32 = 1.8e-18; % Stimulated emission cross-section (m^2)
gamma = 0.32; % confinement factor

tau_32 = 2.1e-12;   % 2.1ps
tau_31 = 4.2e-12;   % 4.2ps
%tau_21 = 0.3e-12;  % 0.3ps  1e-15
tau_sp = 38e-9;     % 38ns 
% Only 2 levels in LaserMatrix. When make tau_21 as small as possible
% convert 3 levels to 2 levels

tau_3 = 1/(1/tau_32+1/tau_31+1/tau_sp);    %  life time of upper level: should be 1.4ps

tau_p = 3.36e-12;    % photon life time in cavity 3.36ps

Xi = tau_21/tau_32 + tau_21/tau_sp;

Ith = q*V/(gamma*vg*sigma_32*tau_p*tau_3*N*(1-Xi));  % threshold current
Jth = Ith/(W*L); % threshold current density