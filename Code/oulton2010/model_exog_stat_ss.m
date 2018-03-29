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
biggami  = param.biggami; % growth rate of equipment productivity
biggamc  = param.biggamc;    % growth rate of neutral technology
di       = param.di; % depreciation rate of capital IT
dc       = param.dc; % depreciation rate of capital standard
chi      = param.chi; % preference parameter.

%Use closed form expressions for the ss values. 
% My notes p. 289-292
p  = (biggamc)/(biggami); % p. 291
rc = 1/bet * (biggamc)^((1-b)/(1-a-b)) * (biggami)^((b)/(1-a-b)) - (1-dc); % p. 289
ri = (1/bet* ((biggamc)^((1-b)/(1-a-b)) *(biggami)^((b)/(1-a-b)))*biggami/biggamc -(1-di))*p; 

% p. 292
kc1_h1 = (a/rc)^((1-b)/(1-a-b)) * (biggamc)^(1/(1-a-b)) * (b/ri)^(b/(1-a-b));%  eq.1 and p. 293 comment
ki1_h1 = (b*(biggamc)/ri)^(1/(1-b)) * kc1_h1^(a/(1-b)); % eq.2 and p. 293 comment
w      = (1-a-b)*(biggamc)*kc1_h1^a * ki1_h1^b; % eq.3 and p. 293 comment
% h1     = w/(chi*(biggamc)*kc1_h1^a*ki1_h1^b - chi*2*kc1_h1*(((biggamc)^((1-b)/(1-a-b)) *(biggami)^((b)/(1-a-b))) - (1-dc))); %eq.4 and p. 293 comment
expg = ((biggamc)^((1-b)/(1-a-b)) *(biggami)^((b)/(1-a-b)));
expgi = ((biggamc)^((a)/(1-a-b)) *(biggami)^((1-a)/(1-a-b)));
g  = expg - (1-dc); 
gi = expgi - (1-di);

h1_h2 = (biggami*kc1_h1^a*ki1_h1^b - ki1_h1*gi)*1/(ki1_h1*gi); % Note: kc2_h2 = kc1_h1
h1    = w/chi *1/(biggamc*kc1_h1^a*ki1_h1^b - (kc1_h1 + kc1_h1/h1_h2)*g);
h2 = h1/h1_h2;

% A check:
% RHS = biggamc/biggami*a/b *ri/p * 1/rc * (rc/(biggamc*a))^(1/b);
% kc1_h1_alt = RHS^(b/(b+a-1));
% ki1_h1_alt = (rc/(biggamc*a))^(1/b)*kc1_h1_alt^((1-a)/b);
% w_alt = (1-a-b)*biggamc*kc1_h1_alt^a*ki1_h1_alt^b;

kc1 = kc1_h1*h1;
ki1 = ki1_h1*h1;

% % Marco's notes:
% LHS = w/(chi*h1) + kc1_h1*(expg - (1-dc)) - biggamc/biggami*ki1_h1*(expgi - (1-di));
% S1  = biggamc/biggami*ki1_h1*(expgi - (1-di)); % gi = expgi - (1-di) in Marco's notes
% S2  = kc1_h1*(expg - (1-dc)); % g = expg - (1-dc) in Marco's notes
% h1_h2 = (LHS + sqrt(LHS^2 + 4*S1*S2))/(2*S1); % the positive root
% h2    = h1/h1_h2;

% h2  = h1;  % See Oulton 2010, p. 10.
% kc2 = kc1;
% ki2 = ki1;
kc2 = kc1_h1*h2; % Because kc2_h2 = kc1_h1
ki2 = ki1_h1*h2; % - || - 

% p. 293
kc = kc1+kc2;
ki = ki1+ki2;
h  = h1+h2;

c      = w/chi;
ic     = (  ((biggamc)^((1-b)/(1-a-b)) *(biggami)^((b)/(1-a-b)))   -(1-dc))*kc; % p. 289
it     = (  ((biggamc)^((a)/(1-a-b)) *(biggami)^((1-a)/(1-a-b)))   -(1-di))*ki;
yc     = c+ic;
yi     = it;

%Put the ss values in a vector consistent with Y and X vectors in model.m
% KC KI BIGGAMC BIGGAMI
xx  = [kc ki biggamc biggami]; 
yy  = [yc yi c ic it w rc ri h h1 h2 kc1 kc2 ki1 ki2 p];
% YC YI C IC IT W RC RI H H1 H2 KC1 KC2 KI1 KI2 P

ss  = [yy xx];
