% MODEL_SS - Return the steady state of the model (computed analytically)
%
% usage:
% 
% [ss, param] = model_ss(param)
% STEADY STATE OF OULTON 2010 WITH GROWTH

function [ss,param] = model_exog_stat_ss(param)

%Parameters from param object
bet      = param.bet;  %Never call anything beta...it's a matlab function
a        = param.a; % share of standard capital
b        = param.b; % share of IT capital
gami     = param.gami; % growth rate of equipment productivity
gamc     = param.gamc;    % growth rate of neutral technology
di       = param.di; % depreciation rate of capital IT
dc       = param.dc; % depreciation rate of capital standard
chi      = param.chi; % preference parameter.

%Use closed form expressions for the ss values. 
% My notes p. 289-292
p  = (1-gamc)/(1+gami); % p. 291
rc = 1/(bet* ((1/(1+gamc))^((1-b)/(1-a-b)) *(1/(1+gami))^((b)/(1-a-b))) ) - (1-dc); % p. 289
ri = (1/(bet* ((1/(1+gamc))^((1-b)/(1-a-b)) *(1/(1+gami))^((b)/(1-a-b))) ) -(1-di))*p; 

% p. 292
kc_h = (a/rc)^((1-b)/(1-a-b)) * (1+gamc)^(1/(1-a-b)) * (b/ri)^(b/(1-a-b));%  eq.1
ki_h = (b*(1+gamc)/ri)^(1/(1-b)) * kc_h^(a/(1-b)); % eq.2
w    = (1-a-b)*(1+gamc)*kc_h^a * ki_h^b; % eq.3
h    = w/(chi*(1+gamc)*kc_h^a*ki_h^b - chi*kc_h*( ...
    ((1/(1+gamc))^((1-b)/(1-a-b)) *(1/(1+gami))^((b)/(1-a-b))) - (1-dc))); %eq.4

kc = kc_h*h;
ki = ki_h*h;
h1 = h; % See Oulton 2010, p. 10.
h2 = h;
kc1= kc;
kc2= kc;
ki1= ki;
ki2= ki;

c      = w/chi;
ic     = (  ((1/(1+gamc))^((1-b)/(1-a-b)) *(1/(1+gami))^((b)/(1-a-b)))   -(1-dc))*kc; % p. 289
it     = (  ((1+gamc)^((a)/(1-a-b)) *(1+gami)^((1-a)/(1-a-b)))   -(1-di))*ki;
yc     = c+ic;
yi     = it;

%Put the ss values in a vector consistent with Y and X vectors in model.m
% KC KI GAMC GAMI
xx  = [kc ki gamc gami]; 
yy  = [yc yi c ic it w rc ri h h1 h2 kc1 kc2 ki1 ki2 p];
% YC YI C IC IT W RC RI H H1 H2 KC1 KC2 KI1 KI2 P

ss  = [yy xx];
