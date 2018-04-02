% MODEL_SS - Return the steady state of the model (computed analytically)
%
% usage:
% 
% [ss, param] = model_ss(param)
% STEADY STATE OF OULTON 2010 WITH GROWTH

function [ss,param] = model_spillover_ss(param)

%Parameters from param object
bet      = param.bet;  %Never call anything beta...it's a matlab function
a        = param.a; % share of standard capital
b        = param.b; % share of IT capital
biggami  = param.biggami; % growth rate of equipment productivity
biggamc  = param.biggamc;    % growth rate of neutral technology
di       = param.di; % depreciation rate of capital IT
dc       = param.dc; % depreciation rate of capital standard
chi      = param.chi; % preference parameter.
gam      = param.gam; % spillover elasticity

%Use closed form expressions for the ss values. 
% My notes p. 289-292
p  = (biggamc)/(biggami); % p. 291
expgc = biggamc^((1-b-gam)/(1-a-b-gam))*biggami^((b+gam)/(1-a-b-gam));
expgi = biggamc^(a/(1-a-b-gam))*biggami^((1-a)/(1-a-b-gam));
gc  = expgc - (1-dc); 
gi = expgi - (1-di);
rc = 1/bet * expgc - (1-dc); % p. 289
ri = (1/bet* expgi -(1-di))*p; 

% What we do know: from Marco's notes and from FOCs simply.
ki_kc = b/a*rc/ri;
w_ki = (1-a-b)/b*ri;


% p. 292 Now these ratios are wrong, but they're used to get the correct
% h1/h2.
kc1_h1 = (a/rc)^((1-b)/(1-a-b)) * (biggamc)^(1/(1-a-b)) * (b/ri)^(b/(1-a-b));%  eq.1 and p. 293 comment
ki1_h1 = (b*(biggamc)/ri)^(1/(1-b)) * kc1_h1^(a/(1-b)); % eq.2 and p. 293 comment
w      = (1-a-b)*(biggamc)*kc1_h1^a * ki1_h1^b; % eq.3 and p. 293 comment


h1_h2 = (biggami*kc1_h1^a*ki1_h1^b - ki1_h1*gi)*1/(ki1_h1*gi); % Hopefully this still holds w/ spillover.
h2 = ((biggamc/biggami*h1_h2 - (ki_kc)^(-1)*(gc/gi)) *chi/w_ki *(1+h1_h2)*gi)^(-1); % This comes from "Guess Marco"
h1 = h1_h2*h2;
h = h1+h2;
yc_yi = h1_h2*biggamc/biggami;
c_ii = yc_yi - ki_kc^(-1)*gc/gi;
w = (1-a-b)*( biggamc* h^gam *(b/ri)^(gam+b) * (a/rc)^a )^1/(1-gam-a-b);
ki1 = b/(1-a-b)*w/ri*h1;
kc1 = a/(1-a-b)*w/rc*h1;
ki2 = ki1/h1*h2;
kc2 = kc1/h1*h2;


% p. 293
kc = kc1+kc2;
ki = ki1+ki2;
h  = h1+h2;

c      = w/chi;
ic     = (  gc)*kc; % p. 289
it     = (  gi)*ki;
yc     = c+ic;
yi     = it;

%Put the ss values in a vector consistent with Y and X vectors in model.m
% KC KI BIGGAMC BIGGAMI
xx  = [kc ki biggamc biggami]; 
yy  = [yc yi c ic it w rc ri h h1 h2 kc1 kc2 ki1 ki2 p expgc expgi];
% YC YI C IC IT W RC RI H H1 H2 KC1 KC2 KI1 KI2 P

ss  = [yy xx];
