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

function [fyn, fxn, fypn, fxpn] = model_spillover(param)

%Steady State
[ss, param] = model_spillover_ss(param);

%Declare parameters
bet      = param.bet;  %Never call anything beta...it's a matlab function
a        = param.a; % share of standard capital
b        = param.b; % share of IT capital
biggami  = param.biggami; % growth rate of equipment productivity
biggamc  = param.biggamc;    % growth rate of neutral technology
di       = param.di; % depreciation rate of capital IT
dc       = param.dc; % depreciation rate of capital standard
chi      = param.chi; % preference parameter.
gam      = param.gam; % spillover elasticity

%Declare Needed Symbols
syms KC KI KC1 KC2 KI1 KI2 % 6 vars
syms KC_p KI_p KC1_p KC2_p KI1_p KI2_p 
syms BIGGAMC BIGGAMI % 2 vars
syms BIGGAMC_p BIGGAMI_p 
syms YC YI C IC IT W RC RI H H1 H2 % 11 vars
syms YC_p YI_p C_p IC_p IT_p W_p RC_p RI_p H_p H1_p H2_p  
syms P P_p % 1 var  % -> 20 vars
syms CL CL_p KIL KIL_p KCL KCL_p YCL YCL_p YIL YIL_p HL HL_p PL PL_p KC2L KC2L_p KI2L KI2L_p
syms GAMC GAMC_p GAMKI GAMKI_p GAMKC GAMKC_p GAMYC GAMYC_p GAMYI GAMYI_p GAMH GAMH_p GAMP GAMP_p GAMKC2 GAMKC2_p GAMKI2 GAMKI2_p
syms BIGGAMCL BIGGAMCL_p BIGGAMIL BIGGAMIL_p
syms EXPGC EXPGC_p EXPGI EXPGI_p

%Declare X and Y vectors
X  = [KC KI BIGGAMC BIGGAMI]; % vector of state variables  
XP = [KC_p KI_p BIGGAMC_p BIGGAMI_p]; % p signifies t+1 

Y  = [YC YI C IC IT W RC RI H H1 H2 KC1 KC2 KI1 KI2 P EXPGC EXPGI]; % vector of controls
YP = [YC_p YI_p C_p IC_p IT_p W_p RC_p RI_p H_p H1_p H2_p KC1_p KC2_p KI1_p KI2_p P_p EXPGC_p EXPGI_p] ;

%Make index variables for future use
make_index([Y,X])

% Model Equations 
f(1)    = -YC + KI^gam * (BIGGAMC)*H1^(1-a-b)*KC1^(a)*KI1^(b); 
f(end+1)= -YI + KI^gam *(BIGGAMI)*H2^(1-a-b)*KC2^(a)*KI2^(b); 
f(end+1)= -KC + KC1 + KC2;
f(end+1)= -KI + KI1 + KI2; %
f(end+1)= -H + H1 + H2;
f(end+1)= -EXPGC + BIGGAMC^((1-b-gam)/(1-a-b-gam))*BIGGAMI^((b+gam)/(1-a-b-gam));
f(end+1)= -EXPGI + BIGGAMC^(a/(1-a-b-gam))*BIGGAMI^((1-a)/(1-a-b-gam));
f(end+1)= -KC_p*(EXPGC) + (1-dc)*KC + IC;
f(end+1)= -KI_p*(EXPGI)+ (1-di)*KI + IT; 
f(end+1)= -YC + C + IC;
f(end+1)= -YI + IT;
f(end+1)= -W/C + chi;
f(end+1)= -1 + bet*C/C_p*1/EXPGC*(RC_p + 1-dc);
f(end+1)= -1 + bet*C/C_p*1/EXPGI*(RI_p/P_p + 1-di);
f(end+1)= -W + KI^gam * (1-a-b)*(BIGGAMC)*H1^(-a-b)*KC1^(a)*KI1^(b); %
f(end+1)= -RC + KI^gam * a*(BIGGAMC)*H1^(1-a-b)*KC1^(a-1)*KI1^(b); %
f(end+1)= -RI + KI^gam * b*(BIGGAMC)*H1^(1-a-b)*KC1^(a)*KI1^(b-1); %
f(end+1)= -W + KI^gam*(1-a-b)*(BIGGAMI)*H2^(-a-b)*KC2^(a)*KI2^(b)*P;%
f(end+1)= -RC + KI^gam*a*(BIGGAMI)*H2^(1-a-b)*KC2^(a-1)*KI2^(b)*P;%
f(end+1)= -RI + KI^gam*b*(BIGGAMI)*H2^(1-a-b)*KC2^(a)*KI2^(b-1)*P;%
f(end+1)= log(BIGGAMC_p/biggamc) - .8*log(BIGGAMC/biggamc); %taken directly from Ryan's example code.
f(end+1)= log(BIGGAMI_p/biggami) - .8*log(BIGGAMI/biggami); %taken directly from Ryan's example code.


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
