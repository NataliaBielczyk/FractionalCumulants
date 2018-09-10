function [time_series, BOLD, binary_inputs, inputs, rho, K, gamma, alpha, tau] = forward_model(A, T, Fs, TR, noise_level, version, noise_type)
% reproduces a computational study by Smith et al. (2011):
% http://www.ncbi.nlm.nih.gov/pubmed/20817103
% runs a DCM forward model on the basis of Friston et al. (2003):
% http://www.ncbi.nlm.nih.gov/pubmed/12948688
% written by Natalia Bielczyk, 2017

% Inputs:
%           A:  the adjacency matrix
%               note that the adjency matrix must have self inhibition on the diagonal 
%               + eigenvalues with negative real parts!!
%           T:  length of the simulation (in minutes, 10min as default)
%           Fs: the sampling frequency (Hz, 200 as a default)
%           TR: TR for subsampling the BOLD response (is seconds, 3s as a default)
%           noise_level: STD of the noise (as a fraction of the value of high input)
%           version: whether you want to use the integration with the full Balloon-Windkessel model ('bw') for the haemodynamics, 
%               or simply convolute with the HRF function ('conv'), 'conv' as default
%               beware: the data generated with use of full BW model ('bw') will be 12 seconds shorter than the data geneated with convolution
%           noise_type: can be chosen from 'white', 'pink', 'red', 'violet', 'blue' ('white' as default), look here: https://en.wikipedia.org/wiki/Colors_of_noise

% Outputs:
%           time_series: neuronal time series (fast dynamics) in each node, nodes in the rows
%           BOLD: structure containing BOLD response in each node (variable 'full')
%               and BOLD subsampled every 3s (variable 'subsampled')
%           binary inputs: the record of the randomly generated on- and off- Poissonian trains
%           inputs: the record of the full inputs (Poissonian trains + noise)
%           rho, K, gamma, alpha, tau: record of the randomly picked sets of hemodynamic parameters

global rho K gamma alpha tau i
warning('off','all');

%% seed the random number generator based on the current time:
rng('shuffle');
display('preparing the simulation:...')
pause(3*rand);  %make sure each instance starts at different time in qsubcellfun

%% set simulation parameters:
lag = round((1/Fs)/20);                         % lag in information transfer between different brain regions, after Smith et al. (2012): lag = 10 frames, which is 50ms for Fs = 200
T = T*60;                                       % express time in seconds
V_0 = 0.02;                                     % haemodynamic parameter
thermal_noise = 0.0;                            % add thermal noise to the resulting BOLD (as the fraction of the mean)
low_input = 0;                                  % value of the stimulus-off state
high_input = 0.05;                              % value of the stimulus-on state
prob_of_switch_on = (1/(7.5*Fs));               % expected length of silence: 7.5s
prob_of_switch_off = (1/(2.5*Fs));              % expected length of the input signal: 2.5s
N = length(A);

%% set the defaults:
if numel(T) == 0
    T = 10;
end
if numel(Fs) == 0
    Fs = 200;
end
if numel(TR) == 0
    TR = 3;
end
if numel(version) == 0
    version = 'conv';
end
if numel(noise_type) == 0
    noise_type = 'white';
end

%% sample a random set of haemodynamic parameters from the dictionary:
indexes = randperm(1000);
par = load('parameters.mat');
N_hrf = size(par.K, 2);
p = randperm(N_hrf);
indexes = p(1:N);
for i = 1:N
    rho(i) = par.rho(indexes(i));
    K(i) = par.K(indexes(i));
    gamma(i) = par.gamma(indexes(i));
    alpha(i) = par.alpha(indexes(i));
    tau(i) = par.tau(indexes(i));
end

%% input connections:
C = ones(N,1);

%% initial condition for inputs (random assignment to low and high states):
kickoff = high_input*round(rand(N,1));
u = kickoff*low_input + (ones(N,1)-kickoff)*high_input;

%% prepare noise trains:
noise_trains = [];
for i = 1:N
    if strcmp(noise_type, 'white') == 1        
        noise_trains = [noise_trains; randn(1, Fs*T)];
    elseif strcmp(noise_type, 'pink') == 1
        noise_trains = [noise_trains; pinknoise(Fs*T)];
    elseif strcmp(noise_type, 'red') == 1
        noise_trains = [noise_trains; rednoise(Fs*T)];
    elseif strcmp(noise_type, 'blue') == 1
        noise_trains = [noise_trains; bluenoise(Fs*T)];
    elseif strcmp(noise_type, 'violet') == 1
        noise_trains = [noise_trains; violetnoise(Fs*T)];   
    end
end
noise_trains = noise_level*(high_input-low_input)*noise_trains;
    
%% initiate neuronal time series, and observation variable:
time_series = zeros(N,round(T*Fs));            %neuronal time series for the connected network
BOLD = struct('full',{},'subsampled',{});      %observation (full 'obs' sampled as frequently as time_series, and BOLD 'subsampled' every TR)    %observation (full 'obs' sampled as frequently as time_series, and BOLD 'subsampled' every TR), flippe
times = struct([]);                            %auxiliary variable for the ode45 solver

%% initiate haemodynamic variables: s = 0, f = v = q = 1;
init = [0 1 1 1];

%% calculate the whole time series:
display('calculating the neuronal time series:...')
time_series(:,1:lag) = zeros(N, lag);          %initial condition, time series defined on the lag period

for i = lag+1:round(T*Fs)
    for j = 1:N
        %define current input:
        if u(j) == low_input
            u(j) = u(j) + (high_input - low_input)*poissrnd(prob_of_switch_on);
        else
            u(j) = high_input;
            u(j) = u(j) - (high_input - low_input)*poissrnd(prob_of_switch_off);
        end
        u_curr = u(j) + normrnd(0, noise_level*(high_input-low_input));
        noisy_inputs(j,i) = u_curr;
        binary_inputs(j,i) = u(j);
        %calculate current activity:
        time_series(j,i) = time_series(j,i-lag) + sum(A(:,j).*time_series(:,i-lag)) + C(j)*u_curr;    
    end
end

display('calculating the BOLD:...')
%% calculate the haemodynamics:
switch version
    case 'bw'
        for i = 1:N % apply the Balloon-Windkessel haemodynamic model to each region apart:
            t = (1/Fs):(1/Fs):T;
            u1 = time_series(i,:);
            fun = @(xc) interp1(t, u1, xc);
            options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6 1e-6 1e-6]);
            [T1,Y1] = ode45(@(t,y) hrf_eq(t,y,fun), [10/Fs T], init, options);
            q = Y1(:,4);
            v = Y1(:,3);
            k11 = 7*rho(i); k12 = 2; k13 = 2*rho(i) - 0.2;
            BOLD(i).obs = V_0*(k11*(1 - q1) + k12*(1 - q1./v1) + k13*(1 - v1));
            times(i).time = T1;
        end
        
        %% add thermal noise to the resulting BOLD if you wish:
        for i = 1:N
            thermal_noise_per_node = normrnd(0, thermal_noise, 1, length(y(i).obs));
            BOLD(i).obs = BOLD(i).obs + thermal_noise_per_node.*repmat(mean(mean(BOLD(i).obs)), 1, length(BOLD(i).obs));
        end
        
        %% subsample BOLD every TR and subtract first 4 frames from the readout BOLD:
        t_vec = TR:TR:TR*floor(T/TR);
        for i = 1:N
            subsampled = zeros(1,length(t_vec));
            subsampled = interp1(times(i).time, y(i).obs, t_vec);
            BOLD(i).subsampled = subsampled(5:end);
        end
        
        %% subsample BOLD as dense as the time series and cut first 12s:
        t_vec = (1/Fs):(1/Fs):T;
        for i = 1:N
            full = zeros(1,length(t_vec));
            full = interp1(times(i).time, BOLD(i).obs, t_vec);
            BOLD(i).full = full(Fs*12+1:end);
        end
    case 'conv'
        for j = 1:N
            rho_i = rho(j); K_i = K(j); gamma_i = gamma(j); alpha_i = alpha(j); tau_i = tau(j);
            shape = hrf_shape(rho_i, K_i, gamma_i, alpha_i, tau_i, Fs);
            shape_cut = shape(21:end);
            beginning = linspace(0.0, shape(20), 21);
            shape_interp = [beginning(1:end-1), shape_cut];
            BOLD(j).full = conv(time_series(j, :), shape_interp, 'same');
            thermal_noise_per_node = normrnd(0, thermal_noise, 1, length(BOLD(j).full));
            BOLD(j).full = BOLD(j).full + thermal_noise_per_node.*repmat(mean(mean(BOLD(j).full)), 1, length(BOLD(j).full));
            t_vec = TR:TR:TR*(floor(T/TR));
            BOLD(j).subsampled = interp1(1/Fs:1/Fs:T, BOLD(j).full, t_vec);
        end
end

inputs          = noisy_inputs;
binary_inputs   = compress_inputs(binary_inputs);