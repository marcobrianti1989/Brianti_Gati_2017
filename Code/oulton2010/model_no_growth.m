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
 
% This is OULTON 2010 W/O GROWTH

function [fyn, fxn, fypn, fxpn] = model_no_growth(param)

%Steady State
[ss, param] = model_ss_no_growth(param);

%Declare parameters
bet      = param.bet;  %Never call anything beta...it's a matlab function
a        = param.a; % share of standard capital
b        = param.b; % share of IT capital
gamq     = param.gamq; % growth rate of equipment productivity
gamz     = param.gamz;    % growth rate of neutral technology
g        = param.g; % growth rate of the rest of the model
di       = param.di; % depreciation rate of capital IT
dc       = param.dc; % depreciation rate of capital standard
chi      = param.chi; % preference parameter.
bc       = param.bc; % treating final goods technology as a param
bi       = param.bi; % treating IT goods technology as a param

%Declare Needed Symbols
syms KC KI KC1 KC2 KI1 KI2 BC BI % 8 vars
syms KC_p KI_p KC1_p KC2_p KI1_p KI2_p BC_p BI_p
syms YC YI C IC IT W RC RI P H H1 H2 % 12 vars
syms YC_p YI_p C_p IC_p IT_p W_p RC_p RI_p P_p H_p H1_p H2_p % 12 vars -> 20 vars -2 = 18 vars

%Declare X and Y vectors
% KC KI 
X  = [KC KI]; % vector of state variables  
XP = [KC_p KI_p]; % p signifies t+1 

% YC YI C IC IT W RC RI P H H1 H2 KC1 KC2 KI1 KI2
Y  = [YC YI C IC IT W RC RI P H H1 H2 KC1 KC2 KI1 KI2]; % vector of controls
YP = [YC_p YI_p C_p IC_p IT_p W_p RC_p RI_p P_p H_p H1_p H2_p KC1_p KC2_p KI1_p KI2_p] ;


% Model Equations (my notes p. 192):
f(1)     = -YC + bc * H1^(1-a-b) * KC1^(a) * KI1^(b);                       % eq 1  %
f(end+1) = -YI + bi * H2^(1-a-b) * KC2^(a) * KI2^(b);                       % eq 2  %
f(end+1) = -KC + KC1 + KC2;                                                 % eq 3
f(end+1) = -KI + KI1 + KI2;                                                 % eq 4
f(end+1) = -H + H1 + H2;                                                    % eq 5
f(end+1) = -KC_p + (1-dc)*KC + IC;                                          % eq 6
f(end+1) = -KI_p + (1-di)*KI + IT;                                          % eq 7  
f(end+1) = -YC + C + IC;                                                    % eq 8
f(end+1) = -YI + IT;                                                        % eq 9
% skip eq 10, 11 and 12 b/c we take p as given now and the technologies too
f(end+1) = -W/C + chi;                                                      % eq 13 
f(end+1) = -1 + bet*C/C_p*(RC_p + 1 - dc);                                  % eq.14 %
f(end+1) = -1 + bet*C/C_p*(RI_p/P_p + 1-di);                                % eq 15 %
f(end+1) = -W + (1-a-b)*bc* H1^(-a-b) * KC1^(a) * KI1^b;                    % eq 16
f(end+1) = -RC + a*bc* H1^(1-a-b) * KC1^(a-1) * KI1^b;                      % eq 17
f(end+1) = -RI + b*bc* H1^(1-a-b) * KC1^(a) * KI1^(b-1);                    % eq 18
f(end+1) = -W + (1-a-b)*bi* H2^(-a-b) * KC2^(a) * KI2^b * P;                % eq 19
f(end+1) = -RC + a*bi* H2^(1-a-b) * KC2^(a-1) * KI2^(b) * P;                  % eq 20  %
f(end+1) = -RI + b*bi* H2^(1-a-b) * KC2^(a) * KI2^(b-1) * P;                % eq 21  -3 = 18 equations %
%f(end+1) = P -1; % might have to add this if we don't have enough eqs.

% Stuff maybe needed for growth version later:
% f(end+1) = log(GAMQ_p/gamq) - .95*log(GAMQ/gamq); %taken directly from Ryan's example code.
% f(end+1) = log(GAMZ_p/gamz) - .95*log(GAMZ/gamz); %taken directly from Ryan's example code.
% f(end+1) = -G + GAMZ^(1/(1-alphe-alphs)) + GAMQ^(alphe/(1-alphe-alphs)); % the relation between growth rates / we'll need something like this   

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
