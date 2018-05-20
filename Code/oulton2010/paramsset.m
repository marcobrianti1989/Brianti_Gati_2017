%% PARAMSSET -- define which parameters are to optimized over, and which ones are set (fix)

function [param,set] = paramsset() %GHK params if anything
set.bet      = 0.95;  % 0.95 Never call anything beta...it's a matlab function
set.a        = 0.13;  %0.13 share of standard capital
set.b        = 0.17; %0.17 share of IT capital
set.biggami  = 1.032; %1.032; % growth rate of IT productivity
set.biggamc  = 1.00328; %1.00328;    % growth rate of cons. productivity
set.di       = 0.124; %0.124; % depreciation rate of IT capital (Some work suggests 0.1830)
set.dc       = 0.056; %0.056 depreciation rate of standard capital
set.chi      = 1; % preference parameter of preference for leisure. (From Ryan)
param.gam      = 0.2; %0.2 spillover elasticity (only used in model_spillover.m and model_spillover_ss.m)
set.siggami  = 1; % variance of BIGGAMI
set.sige     = 2; % variance of signal on BIGGAMI

%Always keep these two settings last
set.adiff = 0;       %Anlytical differentiation of model (no used here)
set.me_eq = [];      %Measurement Error to the 5th observation equation
set.approx_deg = 1;  %Degree of Approximation