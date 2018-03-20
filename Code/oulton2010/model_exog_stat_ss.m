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
p  = (1+gamc)/(1+gami); % p. 291
rc = 1/(bet* ((1/(1+gamc))^((1-b)/(1-a-b)) *(1/(1+gami))^((b)/(1-a-b))) ) - (1-dc); % p. 289
ri = (1/(bet* ((1/(1+gamc))^((1-b)/(1-a-b)) *(1/(1+gami))^((b)/(1-a-b))) ) -(1-di))*p; 

% p. 292
kc1_h = (a/rc)^((1-b)/(1-a-b)) * (1+gamc)^(1/(1-a-b)) * (b/ri)^(b/(1-a-b));%  eq.1 and p. 293 comment
ki1_h = (b*(1+gamc)/ri)^(1/(1-b)) * kc1_h^(a/(1-b)); % eq.2 and p. 293 comment
w    = (1-a-b)*(1+gamc)*kc1_h^a * ki1_h^b; % eq.3 and p. 293 comment
h1    = w/(chi*(1+gamc)*kc1_h^a*ki1_h^b - chi*kc1_h*( ...
    ((1/(1+gamc))^((1-b)/(1-a-b)) *(1/(1+gami))^((b)/(1-a-b))) - (1-dc))); %eq.4 and p. 293 comment

kc1 = kc1_h*h1;
ki1 = ki1_h*h1;
%h1 = h1; % See Oulton 2010, p. 10.
h2 = h1;
%kc1= kc1;
kc2= kc1;
%ki1= ki1;
ki2= ki1;

% p. 293
kc = kc1+kc2;
ki = ki1+ki2;
h  = h1+h2;

c      = w/chi;
ic     = (  ((1+gamc)^((1-b)/(1-a-b)) *(1+gami)^((b)/(1-a-b)))   -(1-dc))*kc1; % p. 289
it     = (  ((1+gamc)^((a)/(1-a-b)) *(1+gami)^((1-a)/(1-a-b)))   -(1-di))*ki1;
yc     = c+ic;
yi     = it;

%Put the ss values in a vector consistent with Y and X vectors in model.m
% KC KI GAMC GAMI
xx  = [kc ki gamc gami]; 
yy  = [yc yi c ic it w rc ri h h1 h2 kc1 kc2 ki1 ki2 p];
% YC YI C IC IT W RC RI H H1 H2 KC1 KC2 KI1 KI2 P

ss  = [yy xx];
