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
bet = param.bet;  %Never call anything beta...it's a matlab function
kapp = param.kapp;  % output-gap response parameter of NKPC
phi = param.phi;  % inflation-response parameter of Taylor-rule
rho = param.rho;  % coefficient of AR(1) tech process


%Declare Needed Symbols
syms XX XX_p ETA ETA_p
syms YY YY_p PI PI_p I I_p A A_p   

%Declare X and Y vectors
X  = [XX ETA]; % vector of state variables  
XP = [XX_p ETA_p]; % p signifies t+1

Y  = [YY PI I A]; % vector of controls
YP = [YY_p PI_p I_p A_p] ;


% Lorenzoni Model Equations (Full Information Model):
f(1) = -YY + YY_p + PI_p - I;
f(2) = -PI + kapp*(YY-A) +bet*PI_p;
f(3) = -I + phi*PI;
f(4) = XX_p - rho*XX;
f(5) = ETA_p;
f(6) = A - (XX + ETA);

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
