clear all;

%Ornstein-Uhlenbeck generator
% OU is  dX = -lambda*X dt + mu dW,   X(0) = Xzero
% variane = mu^2/(2*lamda);

% #### INPUT Parameters ####
d0 = 1.0;                   % expected mean value of the signal
delw = 5;                   % noise parameter: mostly affect correlation length
q = 0.0655;                 % noise parameter: variance of the OU signal
T1 = 0:0.0001:150;    % time for the O-U signal


randn('state',100)  % if one wish to generate truly random data, then comment this line
lamda = delw; mu = sqrt(2*delw*q); 
Xzero = 0;    % problem parameters. Theoritecally this value should belong to a normal distribution of mean 0 and variance (mu^2/(2*lamda))
T = max(T1); N = 2^18; dt = T/N;         
dW = sqrt(dt)*randn(1,N);         % Brownian increments

R = 4; Dt = R*dt; L = N/R;        % L EM steps of size Dt = R*dt
Xem = zeros(1,L);                 % preallocate for efficiency
Xtemp = Xzero;
for j = 1:L
   Winc = sum(dW(R*(j-1)+1:R*j)); 
   Xtemp = Xtemp - lamda*Xtemp*Dt + mu*Winc;
   Xem(j) = Xtemp;
end

t=0:Dt:T;
d=d0+[Xzero,Xem];
u_noisy=interp1(t,d,T1);
u_noisy=u_noisy';

% ### NOTES ###

% Discretized Brownian path over [0,1] has dt = 2^(-20). N = 2^20; dt = T/N
% Euler-Maruyama uses timestep R*dt. ##From Higham##

% if you change max(T1), then you need to change the power of 2 in N =
% 2^20 as well; this will eventually modify dt.

% if you want to regenerate some old data then use the exact same d0, delw, q,
% T1, N, R, randn state seed and Xzero as was used in the old run