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
param.nu       = 0.25;  % Frisch-elasticity
param.alph1    = 0.3;  %capital-share of final good PF
param.alph2    = 0.6; % labor-share of final good PF
param.gam      = 0; % param.alph2; % spillover of IT capital
param.phi      = 0.3; % capital share of IT good PF
param.del1     = 0.3; % depreciation rate of capital
param.del2     = 0.3; % depreciation rate of IT capital
% not implementing growth right now, treat tech variables as parameters
%param.b        = 0.9; % AR-coefficient of IT productivity (lambda)
%param.c        = 0.9; % AR-coefficient of common component of TFP (eta)
param.psi      = 1; % a final-good-specific technology component (should be a process, abstract from that for now)
param.lam      = 1; % IT-specific productivity param
param.eta      = 1; % common component tech param