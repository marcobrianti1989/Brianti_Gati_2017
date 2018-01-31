% MODEL_SS - Return the steady state of the model (computed analytically)
%
% usage:
% 
% [ss, param] = model_ss(param)


function [ss,param] = model_ss(param)

%Parameters from param object
bet      = param.bet;  %Never call anything beta...it's a matlab function
eta      = param.eta; % Frisch
alpha_f  = param.alpha_k;  %capital-share of final good PF 
beta_f   = param.beta_f ; % labor-share of final good PF
gam      = param.gam; % spillover of IT capital
phi      = param.phi; % capital share of IT good PF
delta_k  = param.delta_k; % depreciation rate of capital
delta_it = param.delta_it; % depreciation rate of IT capital 
b        = param.b; % AR-coefficient of IT productivity
c        = param.c; % AR-coefficient of common component of TFP

%Use closed form expressions for the ss values.
R = 1/bet -(1-delta_f);

pi = 0; % shouldn't zero st.st. inflation imply pi = 1?
y = 0;
i = 0;
a = 0; 
x = 0;
eta = 0;

%Put the ss values in a vector consistent with Y and X vectors in model.m
yy  = [y pi i a];
xx  = [x eta];
ss  = [yy xx];
