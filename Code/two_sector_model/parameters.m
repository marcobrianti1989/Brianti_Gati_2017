% PARAMETETRS - This function returns a parameter structure to use in the model solution.
% Object is used for easier passsing of many input values to function. 
%
% usage
%
% param = parameters
%
% (No input arguments, file is just for storing values and creating object)

function param = parameters()


param.bet      = 0.98;  %Never call anything beta...it's a matlab function
param.eta      = 0.25;  % Frisch-elasticity
param.alpha_f  = 0.3;  %capital-share of final good PF
param.beta_f   = 0.3; % labor-share of final good PF
param.gam      = 1.1; % spillover of IT capital
param.phi      = 0.3; % capital share of IT good PF
param.delta_k  = 0.1; % depreciation rate of capital
param.delta_it = 0.1; % depreciation rate of IT capital 
param.b        = 0.9; % AR-coefficient of IT productivity
param.c        = 0.9; % AR-coefficient of common component of TFP

