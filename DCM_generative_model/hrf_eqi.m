function dy = hrf_eqi(t,y,fun)
% just equations for Balloon-Windkessel model
global rho_i;
global K_i;
global gamma_i;
global alpha_i;
global tau_i;

dy = zeros(4,1);    
E = 1 - mpower(1 - rho_i, 1/y(2));                                       %auxiliary variable
dy(1) = fun(t) - K_i*y(1) - gamma_i*(y(2) - 1);                          %s
dy(2) = y(1);                                                            %f
dy(3) = (1/tau_i)*(y(2) - mpower(y(3), 1/alpha_i));                      %v
dy(4) = (1/tau_i)*(y(2)*E/rho_i - mpower(y(3), 1/alpha_i)*y(4)/y(3));    %q
