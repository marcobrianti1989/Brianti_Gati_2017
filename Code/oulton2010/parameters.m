% PARAMETETRS - This function returns a parameter structure to use in the model solution.
% Object is used for easier passsing of many input values to function. 
%
% usage
%
% param = parameters
%
% (No input arguments, file is just for storing values and creating object)

function param = parameters() %GHK params if anything
param.bet      = 0.95;  %Never call anything beta...it's a matlab function
param.a        = 0.13;  % share of standard capital
param.b        = 0.17; % share of IT capital
param.gami     = 1.032; %1.032; % growth rate of IT productivity
param.gamc     = 1.00328; %1.00328;    % growth rate of cons. productivity
param.di       = 0.124; %0.124; % depreciation rate of IT capital (Some work suggests 0.1830)
param.dc       = 0.056; %0.056 depreciation rate of standard capital
param.chi      = 1; % preference parameter in labor disutility. (From Ryan)
