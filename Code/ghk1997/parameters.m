% PARAMETETRS - This function returns a parameter structure to use in the model solution.
% Object is used for easier passsing of many input values to function. 
%
% usage
%
% param = parameters
%
% (No input arguments, file is just for storing values and creating object)

function param = parameters()
param.bet      = 0.95;  %Never call anything beta...it's a matlab function
param.alphe    = 0.17;  % share of equipment
param.alphs    = 0.13; % share of structures
param.gamq     = 1.032; % growth rate of equipment productivity
param.gamz     = 1.00328;    % growth rate of neutral technology
param.g        = 1.0124;%1.0124; % growth rate of the rest of the model
param.dele     = 0.124; % depreciation rate of capital equipment
param.dels     = 0.056; % depreciation rate of capital structures
param.thet     = 0.4; % preference parameter.
param.d        = 1; %1 the steady state ratio of z*q^alphe / (y^(1-alphe-alphs)). Open question.