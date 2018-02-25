% MODEL_REDUCED_SS - Return the steady state of the model (computed analytically)
%
% usage:
% 
% [ss, param] = model_ss(param)


function [ss,param] = model_reduced_ss(param)

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

rs     = (bet/gamz)^(-1) -(1-dels); % 0.1339 with GHK parameterization %
y_ks   = rs/alphs;
is_ks  = gamz - (1-dels);
c_ks   = y_ks-is_ks;
l_y    = y_ks^(alphs/(1-alphs)) /gamz;
w      = (1-alphs)*(l_y)^(-1);
c      = w*(thet/(1-thet))*(1-l);
y      = l/(l_y); %
ks     = (y_ks)^(-1)*y;
is     = is_ks*ks;


%Put the ss values in a vector consistent with Y and X vectors in model.m
% y c l w rs re ie is
yy  = [y c l w rs is];
xx  = [ks gamz]; 
% ke ks gamq gamz g q
ss  = [yy xx];
