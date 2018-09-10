function [shape] = hrf_shape(rho_j, K_j, gamma_j, alpha_j, tau_j, Fs)
% simulates 30s of hemodynamic response to an impulse stimulus, 
% given haemodynamic parameters and sampling resolution Fs [Hz]

global rho_i
global K_i
global gamma_i
global alpha_i
global tau_i

rho_i = rho_j;
K_i = K_j;
gamma_i = gamma_j;
alpha_i = alpha_j;
tau_i = tau_j;

V_0   = 0.02; 
init  = [0 1 1 1];
T     = 30; %[s]

if numel(Fs)==0
    Fs = 200; %[Hz]
end

t = (1/Fs):(1/Fs):T;
% create impulse signal:
u1 = zeros(length(t), 1); u1(1:40) = 1; 

% evaluate HRF response to this impuls signal:
fun1 = @(xc) interp1(t, u1, xc);
options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6 1e-6]);
[T1,Y1] = ode45(@(t,y1) hrf_eqi(t,y1,fun1), [10/Fs T], init, options);
q1 = Y1(:,4);
v1 = Y1(:,3);
k11 = 7*rho_i;
k12 = 2;
k13 = 2*rho_i - 0.2;
y1_obs = V_0*(k11*(1 - q1) + k12*(1 - q1./v1) + k13*(1 - v1));
shape = interp1(T1, y1_obs, t);
