% MODEL_SS - Return the steady state of the model (computed analytically)
%
% usage:
% 
% [ss, param] = model_ss(param)


function [ss,param] = model_ss_no_growth(param)

%Parameters from param object
bet      = param.bet;  %Never call anything beta...it's a matlab function
a        = param.a; % share of standard capital
b        = param.b; % share of IT capital
gamq     = param.gamq; % growth rate of equipment productivity
gamz     = param.gamz;    % growth rate of neutral technology
g        = param.g; % growth rate of the rest of the model
di       = param.di; % depreciation rate of capital IT
dc       = param.dc; % depreciation rate of capital standard
chi      = param.chi; % preference parameter.
bc       = param.bc; % treating final goods technology as a param
bi       = param.bi; % treating IT goods technology as a param

%Use closed form expressions for the ss values. 

% % Version 1. Assume P given. (copied from steadystate_no_growth.nb, based on my notes p. 196.)
% p = 1;
% rc     = 1/bet - (1 - dc);
% ri     = (1/bet - (1 - di))*p;
% ki1ki2 = (bc/bi)^((1 - 2*a*b)/b) * p^((1 - 2*(a + b) - 2*a*(a + b))/(a + b));
% kc1kc2 = p^(-1/b)*(bc/bi)^(1/b)*(ki1ki2)^(-1);
% h1h2   = (1/p * bc/bi * kc1kc2^a * ki1ki2^b)^(1/(a+b));
% h1     = ((1-a-b)/(chi)* rc/a) / (rc/a - dc - dc*kc1kc2^-1); % p. 204
% h2     = h1/h1h2;
% kc2ki2 = (di*ki1ki2 + di) / (bi *rc/(a*bi*p));
% ki2    = (rc/(a*bc)* kc1kc2^(1-a)* ki1ki2^(-b) * h1^(a+b-1) * kc2ki2^(a-1) )^(1/(1-a-b));
% kc2    = kc2ki2 * ki2;
% ki1    = ki1ki2 * ki2;
% kc1    = kc1kc2 * kc2;
% w      = (1-a-b)*bc* h1^(-a-b) * kc1^(a) * ki1^b;
% h      = h1 + h2;
% kc     = kc1 + kc2;
% ki     = ki1 + ki2;
% c      = w/chi;
% ic     = dc*kc;
% it     = di*ki;
% yc     = c+ic;
% yi     = it;

% Version 2. Don't assume anything about P. Eq numbers refer to blue
% numbers frm my notes, starting p. 209. Red eq numbers refer to results of my rechecking starting on p. 213.
% The errors are greater than in the previous formulation. 
p      = bc/bi; %  red eq 2
rc     = 1/bet - (1 - dc); %  red eq 1
ri     = (1/bet - (1 - di))*p; % red eq 3
ki1_ki2= (ri/p)/(bet*di); % blue eq 3
kc1_kc2= ki1_ki2; % From red eq. 2
h1_h2  = ki1_ki2; % From red eq. 2
w_kc1  = chi*(rc/a - dc*(kc1_kc2)^(-1) - dc); % red eq 4
h1     = ((1-a-b)/chi * rc/a)/(rc/a -dc -dc*kc1_kc2^-1); % blue eq 5 = red eq 5
h2     = kc1_kc2*h1;  % blue = red eq 6

ki1_h1= ((b*bc/ri)*(a*bc/rc)^a)^((1-a)/(1-a-b)); % red eq 7
ki1   = ki1_h1*h1;
ki2   = ki1/ki1_ki2;
kc1_h1= (a*bc/rc)^(1/(1-a))*ki1_h1^(b/(1-a)); % red eq 10
kc1   = kc1_h1*h1;
kc2   = kc1/kc1_kc2;
w      =(1-a-b)*bc*h1^(-a-b)*kc1^a*ki1^b; % blue eq14

h      = h1 + h2;
kc     = kc1 + kc2;
ki     = ki1 + ki2;
c      = w/chi;
ic     = dc*kc;
it     = di*ki;
yc     = c+ic;
yi     = it;

%Put the ss values in a vector consistent with Y and X vectors in model.m
% KC KI
xx  = [kc ki ]; 
yy  = [yc yi c ic it w rc ri p h h1 h2 kc1 kc2 ki1 ki2];
% YC YI C IC IT W RC RI P H H1 H2 KC1 KC2 KI1 KI2

ss  = [yy xx];
