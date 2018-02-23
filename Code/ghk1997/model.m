% MODEL_LINEAR - Linearized version of the simple RBC model example from 
% EC860.
%
% usage:
%
% [fy, fx, fyp, fxp] = model(param)
%
% where
%
% param = a parameters object, created in parmeters.m
%
% NOTES: This program, which generates the model matrices in cannonical form,
%        requires the MATLAB symbolic toolbox! The algorithm used to solve 
%        the model, however, is numerical and requires no m-files
%        beyond those that appear in this directory.
%
% Code by Ryan Chahrour, Boston College, 2014
 

function [fyn, fxn, fypn, fxpn] = model(param)

%Steady State
[ss, param] = model_ss(param);

%Declare parameters
bet      = param.bet;  %Never call anything beta...it's a matlab function
alphe    = param.alphe; % share of equipment
alphs    = param.alphs; % share of structures
gamq     = param.gamq; % growth rate of equipment productivity
gamz     = param.gamz;    % growth rate of neutral technology
g        = param.g; % growth rate of the rest of the model
dele     = param.dele; % depreciation rate of capital equipment
dels     = param.dels; % depreciation rate of capital structures
thet     = param.thet; % preference parameter.
d        = param.d;    % the steady state ratio of z*q^alphe / (y^(1-alphe-alphs)). Open question.

%Declare Needed Symbols
syms KE KS GAMQ G GAMZ
syms KE_p KS_p GAMQ_p G_p GAMZ_p
syms YY C L RE RS W IE IS Q
syms YY_p C_p L_p RE_p RS_p W_p IE_p IS_p Q_p
syms QRE QRE_p

%Declare X and Y vectors
% KE KS GAMQ GAMZ G Q
X  = [KE KS GAMQ GAMZ]; % vector of state variables  
XP = [KE_p KS_p GAMQ_p GAMZ_p]; % p signifies t+1 

% YY C L W RS RE IE IS
Y  = [YY C L W RS RE IE IS Q G ]; % vector of controls
YP = [YY_p C_p L_p W_p RS_p RE_p IE_p IS_p Q_p G_p ] ;


% Model Equations (my notes p. 160):
f(1)     = -W/C + ((1-thet)/thet)*(1/(1-L));                                % eq 1
f(end+1) = -1 + (bet/G)*(C/C_p)*(RS_p + 1-dels);                            % eq 2
f(end+1) = -1 + (bet/G)*(C/C_p)*(1/GAMQ_p)*(Q_p*RE_p + 1-dele);             % eq 3
f(end+1) = -W/YY + (1-alphe-alphs)/L;                                        % eq 4
f(end+1) = -YY/KS + RS/alphs;                                                % eq 5
f(end+1) = -alphe*YY/KE*GAMQ + Q*RE;                                         % eq 6
f(end+1) = -YY + d*KE^alphe*KS^alphs*L^(1-alphe-alphs);                      % eq 7
f(end+1) = -YY + C + IE + IS;                                                % eq 8
f(end+1) = - IS + KS*(G - (1-dels));                                        % eq 9
f(end+1) = - IE + KE*(G - (1-dele)/GAMQ);                                   % eq 10
f(end+1) = RS - Q*RE;                                                       % eq 11
f(end+1) = log(GAMQ_p/gamq) - .95*log(GAMQ/gamq); %eq. 12, taken directly from Ryan's example code.
f(end+1) = log(GAMZ_p/gamz) - .95*log(GAMZ/gamz); %eq. 14, taken directly from Ryan's example code.
f(end+1) = -G + GAMZ^(1/(1-alphe-alphs)) + GAMQ^(alphe/(1-alphe-alphs));    % eq. 13

%Check Computation of Steady-State Numerically
fnum = double(subs(f, [Y X YP XP], [ss, ss]));
disp('Checking steady-state equations:')
disp(fnum);

%Log-linear approx
log_var = [];
f = subs(f, log_var, exp(log_var)); 
   
%Differentiate
fx  = jacobian(f,X);
fy  = jacobian(f,Y);
fxp = jacobian(f,XP); 
fyp = jacobian(f,YP);

%Plug back into levels
fx =  subs(fx , log_var, log(log_var));
fy =  subs(fy , log_var, log(log_var));
fxp = subs(fxp, log_var, log(log_var));
fyp = subs(fyp, log_var, log(log_var));

%Compute numerical values
fxn =  double(subs(fx , [Y X YP XP], [ss, ss]));
fyn =  double(subs(fy , [Y X YP XP], [ss, ss]));
fxpn = double(subs(fxp, [Y X YP XP], [ss, ss]));
fypn = double(subs(fyp, [Y X YP XP], [ss, ss]));
