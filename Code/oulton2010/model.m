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
 
% This is OULTON 2010 WITH GROWTH

function [fyn, fxn, fypn, fxpn] = model(param)

%Steady State
[ss, param] = model_ss(param);

%Declare parameters
bet      = param.bet;  %Never call anything beta...it's a matlab function
a        = param.a; % share of standard capital
b        = param.b; % share of IT capital
gami     = param.gami; % growth rate of equipment productivity
gamc     = param.gamc;    % growth rate of neutral technology
di       = param.di; % depreciation rate of capital IT
dc       = param.dc; % depreciation rate of capital standard
chi      = param.chi; % preference parameter.

%Declare Needed Symbols
syms KC KI KC1 KC2 KI1 KI2 % 6 vars
syms KC_p KI_p KC1_p KC2_p KI1_p KI2_p 
syms GAMC GAMI G GI GP % 5 vars
syms GAMC_p GAMI_p G_p GI_p GP_p
syms YC YI C IC IT W RC RI H H1 H2 % 11 vars
syms YC_p YI_p C_p IC_p IT_p W_p RC_p RI_p H_p H1_p H2_p  % -> 22 vars

%Declare X and Y vectors
% KC KI GAMC GAMI
X  = [KC KI GAMC GAMI]; % vector of state variables  
XP = [KC_p KI_p GAMC_p GAMI_p]; % p signifies t+1 

% YC YI C IC IT W RC RI H H1 H2 KC1 KC2 KI1 KI2 G GI GP
Y  = [YC YI C IC IT W RC RI H H1 H2 KC1 KC2 KI1 KI2 G GI GP]; % vector of controls
YP = [YC_p YI_p C_p IC_p IT_p W_p RC_p RI_p H_p H1_p H2_p KC1_p KC2_p KI1_p KI2_p G_p GI_p GP_p] ;

% Model Equations (my notes p. 228):
f(1)    = -YC + (1+GAMC)*H1^(1-a-b)*KC1^(a)*KI1^(b);
f(end+1)= -YI + (1+GAMI)*H2^(1-a-b)*KC2^(a)*KI2^(b); 
f(end+1)= -KC + KC1 + KC2;
f(end+1)= -KI + KI1 + KI2;
f(end+1)= -H + H1 + H2;
f(end+1)= -KC_p*(1+G) + (1-dc)*KC + IC;
f(end+1)= -KI_p*(1+GI)+ (1-di)*KI + IT; 
% f(end+1)= -KI_p*(1+GI)*(1+GP_p)+ (1-di)*KI + IT; % brown attempt 3
f(end+1)= -YC + C + IC;
f(end+1)= -YI + IT;
f(end+1)= -W/C + chi;
f(end+1)= -1 + bet*C/C_p*1/(1+G)*(RC_p + 1-dc);
f(end+1)= -1 + bet*C/C_p*1/(1+G)*(RI_p/(1+GP_p) + 1-di); 
% f(end+1)= -1 + bet*C/C_p*1/(1+G)*(RI_p + 1-di); % Brown option 2
f(end+1)= -W + (1-a-b)*(1+GAMC)*H1^(-a-b)*KC1^(a)*KI1^(b); 
f(end+1)= -RC + a*(1+GAMC)*H1^(1-a-b)*KC1^(a-1)*KI1^(b);
f(end+1)= -RI + b*(1+GAMC)*H1^(1-a-b)*KC1^(a)*KI1^(b-1);
f(end+1)= -W + (1-a-b)*(1+GAMI)*(1+GP)*H2^(-a-b)*KC2^(a)*KI2^(b); %%% <-
f(end+1)= -RC + a*(1+GAMI)*(1+GP)*H2^(1-a-b)*KC2^(a-1)*KI2^(b);
f(end+1)= -RI + b*(1+GAMI)*(1+GP)*H2^(1-a-b)*KC2^(a)*KI2^(b-1);
f(end+1)= log(GAMC_p/gamc) - .95*log(GAMC/gamc); %taken directly from Ryan's example code.
f(end+1)= log(GAMI_p/gami) - .95*log(GAMI/gami); %taken directly from Ryan's example code.
% f(end+1)= -GP + GAMC -GAMI; % eq. I p. 222
f(end+1)= -G + ((1-b)*GAMC + b*GAMI)/(1-a-b);  % eq. IIs
f(end+1)= -GI + ((a)*GAMC + (1-a)*GAMI)/(1-a-b);  % eq. III

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
