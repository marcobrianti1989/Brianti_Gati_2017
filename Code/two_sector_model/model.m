% MODEL_LINEAR - Linearized version of the simple RBC model example from 
% EC860.
%
% usage:
%
% [fy, fx, fyp, fxp] = model(param)
%
% where
%
% param = a parameters object, created in parmeters.m
%
% NOTES: This program, which generates the model matrices in cannonical form,
%        requires the MATLAB symbolic toolbox! The algorithm used to solve 
%        the model, however, is numerical and requires no m-files
%        beyond those that appear in this directory.
%
% Code by Ryan Chahrour, Boston College, 2014
 

function [fyn, fxn, fypn, fxpn] = model(param)

%Steady State
[ss, param] = model_ss(param);

%Declare parameters
bet      = param.bet;  %Never call anything beta...it's a matlab function
nu       = param.nu; % Frisch
alph1    = param.alph1;  %capital-share of final good PF 
alph2    = param.alph2 ; % labor-share of final good PF
gam      = param.gam; % spillover of IT capital
phi      = param.phi; % capital share of IT good PF
del1     = param.del1; % depreciation rate of capital
del2     = param.del2; % depreciation rate of IT capital 
%b        = param.b; % AR-coefficient of IT productivity
%c        = param.c; % AR-coefficient of common component of TFP
psi      = param.psi; % a final-good-specific technology component (should be a process, abstract from that for now)
lam      = param.lam; % IT productivity param
eta      = param.eta; % common tech param

%Declare Needed Symbols
syms K K1 K2 S
syms K_p K1_p K2_p S_p
syms YY C L L1 L2 R W I IT P
syms YY_p C_p L_p L1_p L2_p R_p W_p I_p IT_p P_p

%Declare X and Y vectors
X  = [K K1 K2 S]; % vector of state variables  
XP = [K_p K1_p K2_p S_p]; % p signifies t+1

Y  = [YY C L L1 L2 R W I IT P]; % vector of controls
YP = [YY_p C_p L_p L1_p L2_p R_p W_p I_p IT_p P_p] ;


% Model Equations (my notes p. 152):
f(1)     = bet*C/C_p*(R_p + 1 - del1) -1; % eq 1
f(end+1) = L^(nu-1)/W - 1/C;              % eq 2
f(end+1) = alph1*YY/K1 - R;               % eq 3
f(end+1) = alph2*YY/L1 - W;               % eq 4
f(end+1) = (1-alph1-alph2)*YY/S -P;       % eq 5
f(end+1) = lam*(1-phi)*IT/L2*P - W;       % eq 6
f(end+1) = lam*phi*IT/K2*P -R;            % eq 7
f(end+1) = -K + K1 +K2;                   % eq 8
f(end+1) = -L + L1 +L2;                   % eq 9
f(end+1) = (1-del2)*S + IT - S_p;         % eq 10
f(end+1) = (1-del1)*K + I - K_p;          % eq 11
f(end+1) = W*L + R*K + P*S - C -I;        % eq 12
f(end+1) = -YY +K1^(alph1)*L1^(alph2)*...
    S^(1-alph1-alph2)*eta*psi*S^gam;      % eq 13 % *S^gam causes it not to solve right now (explodes)
f(end+1) = -IT*lam*eta*K2^phi*L2^(1-phi); % eq 14

%Check Computation of Steady-State Numerically
fnum = double(subs(f, [Y X YP XP], [ss, ss]));
disp('Checking steady-state equations:')
disp(fnum);

%Log-linear approx
log_var = [];
f = subs(f, log_var, exp(log_var)); 
   
%Differentiate
fx  = jacobian(f,X);
fy  = jacobian(f,Y);
fxp = jacobian(f,XP); 
fyp = jacobian(f,YP);

%Plug back into levels
fx =  subs(fx , log_var, log(log_var));
fy =  subs(fy , log_var, log(log_var));
fxp = subs(fxp, log_var, log(log_var));
fyp = subs(fyp, log_var, log(log_var));

%Compute numerical values
fxn =  double(subs(fx , [Y X YP XP], [ss, ss]));
fyn =  double(subs(fy , [Y X YP XP], [ss, ss]));
fxpn = double(subs(fxp, [Y X YP XP], [ss, ss]));
fypn = double(subs(fyp, [Y X YP XP], [ss, ss]));
