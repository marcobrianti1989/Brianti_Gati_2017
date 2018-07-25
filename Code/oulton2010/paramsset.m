%% PARAMSSET -- define which parameters are to optimized over, and which ones are set (fix)

function [param,set] = paramsset() %GHK params if anything
set.bet        = 0.9924;  % 0.95 Never call anything beta...it's a matlab function
set.a          = 0.297;  %0.13 share of standard capital (should be param eventually)
set.b          = 0.031; %0.17 share of IT capital (should be param eventually)
set.biggami    = 1.016; %1.032; % growth rate of IT productivity
set.biggamc    = 1.0034; %1.00328;    % growth rate of cons. productivity
set.di         = 0.0398; %0.124; % depreciation rate of IT capital (Some work suggests 0.1830)
set.dc         = 0.0206; %0.056 depreciation rate of standard capital
set.chi        = 1; % preference parameter of preference for leisure. (From Ryan)
param.gam      = 0.2; %0.2 spillover elasticity (only used in model_spillover.m and model_spillover_ss.m)
set.siggami    = 1; % variance of BIGGAMI
set.rhogami    = 0.8;
set.sige       = 1; %2 variance of signal on BIGGAMI
param.sigitlev = 1; % variance of the IT prod. shock (level)
param.rhoitlev = 0.8; % persistence of the IT prod shock


%Always keep these two settings last
set.adiff = 0;       %Anlytical differentiation of model (no used here)
set.me_eq = [];      %Measurement Error to the 5th observation equation
set.approx_deg = 1;  %Degree of Approximation