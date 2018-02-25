% MODEL_REDUCED_SS - Return the steady state of the model (computed analytically)
%
% usage:
% 
% [ss, param] = model_ss(param)


function [ss,param] = model_very_reduced_ss(param)

%Parameters from param object
bet      = param.bet;  %Never call anything beta...it's a matlab function
alphe    = param.alphe; % share of equipment
alphs    = param.alphs; % share of structures
gamq     = param.gamq; % growth rate of equipment productivity
gamz     = param.gamz;    % growth rate of neutral technology
g        = param.g; % growth rate of the rest of the model
dele     = param.dele; % depreciation rate of capital equipment
dels     = param.dels; % depreciation rate of capital structures
thet     = param.thet; % preference parameter.
d        = param.d;    % the steady state ratio of z*q^alphe / (y^(1-alphe-alphs)). Open question, but seems to work ok so far.
chi      = param.chi;  % disutility of labor, taken from Ryan. 

%Use closed form expressions for the ss values. (My notes p. 170) - it works!
rs   = (bet/gamz)^(-1) -(1-dels);
ks_y = alphs/rs;
is_y = (gamz - (1-dels))*ks_y; %
c_y  = 1-is_y;
w_y  = chi*c_y;
l    = (1-alphs)/w_y;
ks   = ks_y^(1/(1-alphs))*gamz*l;
y    = ks/ks_y;
c    = (1-is_y)*y;
w    = w_y*y;
is   = is_y *y;

% %Ryan's st.st.
% rs = 1/(bet*gamz^-1)-1+dels;
% kl = (rs/alphs)^(1/(alphs-1));
% w = (1-alphs)*(kl)^alphs*gamz;
% c = w/chi;
% ck = kl^(alphs-1)-(gamz-1+dels);
% ks = ck^-1*c;
% is = (gamz-1+dels)*ks;
% l = kl^-1*ks/gamz;
% y = ks^alphs*(gamz*l)^(1-alphs); % I had to add this.

%Put the ss values in a vector consistent with Y and X vectors in model.m
% y c l w rs re ie is
yy  = [y c l w rs is];
xx  = [ks gamz]; 
% ke ks gamq gamz g q
ss  = [yy xx];
