function dy = hrf_eq(t,y,fun)
% just equations for Balloon-Windkessel model
global rho;
global K;
global gamma;
global alpha;
global tau;
global i;

dy = zeros(4,1);    
E = 1 - mpower(1 - rho(i), 1/y(2));                                        %auxiliary variable
dy(1) = fun(t) - K(i)*y(1) - gamma(i)*(y(2) - 1);                          %s
dy(2) = y(1);                                                              %f
dy(3) = (1/tau(i))*(y(2) - mpower(y(3), 1/alpha(i)));                      %v
dy(4) = (1/tau(i))*(y(2)*E/rho(i) - mpower(y(3), 1/alpha(i))*y(4)/y(3));   %q
