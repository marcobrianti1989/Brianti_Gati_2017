% MODEL_SS - Return the steady state of the model (computed analytically)
%
% usage:
% 
% [ss, param] = model_ss(param)


function [ss,param] = model_ss(param)

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

%Use closed form expressions for the ss values. (My notes p. 147)
l  = 0.24; % labor in the steady state is given.
re = 0.1339;  % I just plugged some number in here... open question and seems to cause problems.

rs     = (bet/g)^(-1) -(1-dels); % 0.1339 with GHK parameterization %
req    = ((bet/g)^(-1))*gamq -(1-dele); 
q      = req/re;
ks_y   = alphs/rs;
ke_y   = alphe*gamq/req; 
is_y   = (g-(1-dels))*ks_y;
ie_y   = (g - (1-dele)/gamq)*ke_y;
c_y    = (thet/(1-thet))*(1/(1-alphe-alphs))*(1-l)/l; 
y      = ((ke_y^alphe + ks_y^alphs)/(c_y + ie_y + is_y))^(1/(1-alphe-alphs))*l; %
w      = (1-alphe-alphs)*y/l;
c      = c_y*y; %
ie     = ie_y*y;
is     = is_y*y;
ke     = ke_y*y;%
ks     = ks_y*y; %

%Put the ss values in a vector consistent with Y and X vectors in model.m
% y c l w rs re ie is
yy  = [y c l w rs re ie is q g ];
xx  = [ke ks gamq gamz]; 
% ke ks gamq gamz g q
ss  = [yy xx];
