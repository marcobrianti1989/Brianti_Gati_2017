% MODEL_SS - Return the steady state of the model (computed analytically)
%
% usage:
% 
% [ss, param] = model_ss(param)
% STEADY STATE OF OULTON 2010 WITH GROWTH AND WITH NEWS SHOCKS

function [ss,param] = model_news_ss(param)

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
v1 = 0;
v2 = 0;

% starting at my notes p. 231, green equations.
gp = gamc - gami; % eq. I p. 222
g = ((1-b)*gamc + b*gami)/(1-a-b);  % eq. II
gi = ((a)*gamc + (1-a)*gami)/(1-a-b);  % eq. III

rc = (1+g)/bet - (1-dc); % green 1
% ri = ((1+g)/bet - (1-di))*(1+gp); % green 2
ri = ((1+g)/bet -(1-di)*(1+gp)); % yellow 2

% Purple suggestion instead of orange one
ki1_h1  = (b/ri)^((1-a)/(1-a-b)) *(a/rc)^(a/(1-a-b)) * (1+gamc)^(1/(1-a-b)); %purple 1
kc1_h1  = (a*(1+gamc)/rc)^(1/(1-a)) *ki1_h1^(b/(1-a)); % purple 2
w       = (1-a-b)*(1+gamc)*kc1_h1^(a) * ki1_h1^(b); % purple 3
kc2_h2  = (a/rc)^((1-b)/(1-a-b)) * (b/ri)^(b/(1-a-b)) *((1+gami)*(1+gp))^(1/(1-a-b)); % purple 4 = orange 1
ki2_h2  = (b/ri * (1+gami)*(1+gp))^(1/(1-b))*(kc2_h2)^(a/(1-b)); % purple 5 = orange 2
ki1_h2  = ((1+gami)/(1+gi -(1-di)))*kc2_h2^a*ki2_h2^b - ki2_h2; % purple 6 alternative
h1_h2   = ki1_h2 / ki1_h1; % purple 7
kc2_h1  = kc2_h2/ h1_h2;  % purple 7
h1      = (w/chi)/((1+gamc)*(kc1_h1)^a*(ki1_h1)^b -(1+g-(1-dc))*(kc1_h1+kc2_h1)); % purple 8
h2      = h1/h1_h2;
ki1     = ki1_h1*h1;
kc1     = kc1_h1*h1;
kc2     = kc2_h2*h2;
ki2     = ki2_h2*h2;

h      = h1 + h2;
kc     = kc1 + kc2;
ki     = ki1 + ki2;
c      = w/chi;
ic     = (1+g-(1-dc))*kc;
it     = (1+gi-(1-di))*ki;
yc     = c+ic;
yi     = it;

%Put the ss values in a vector consistent with Y and X vectors in model.m
% KC KI GAMC GAMI
xx  = [kc ki gamc gami v1 v2]; 
yy  = [yc yi c ic it w rc ri h h1 h2 kc1 kc2 ki1 ki2 g gi gp];
% YC YI C IC IT W RC RI H H1 H2 KC1 KC2 KI1 KI2 G GI GP

ss  = [yy xx];
