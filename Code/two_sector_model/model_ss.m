% MODEL_SS - Return the steady state of the model (computed analytically)
%
% usage:
% 
% [ss, param] = model_ss(param)


function [ss,param] = model_ss(param)

%Parameters from param object
bet      = param.bet;  %Never call anything beta...it's a matlab function
nu       = param.nu; % Frisch
alph1    = param.alph1;  %capital-share of final good PF 
alph2    = param.alph2 ; % labor-share of final good PF
gam      = param.gam; % spillover of IT capital
phi      = param.phi; % capital share of IT good PF
del1     = param.del1; % depreciation rate of capital
del2     = param.del2; % depreciation rate of IT capital 
%b        = param.b; % AR-coefficient of IT productivity
%c        = param.c; % AR-coefficient of common component of TFP
psi      = param.psi; % a final-good-specific technology component (should be a process, abstract from that for now)
lam      = param.lam; % IT productivity param
eta      = param.eta; % common tech param

%Use closed form expressions for the ss values. (My notes p. 147)
r     = 1/bet -(1-del1);
y_k1  = r/alph1;
wl1_y  = alph2;
sp_y  = (1-alph1-alph2)/r;
sp_k2 = r/(eta*lam*phi*del2);
k2_k1 = (eta*lam*phi*del2)/alph1*sp_y; % (eq. 10)
c_y   = wl1_y + (r-del1)*(y_k1)^(-1) + sp_y - del2*k2_k1*(y_k1)^(-1); % (eq. 11)   
wl2_y = eta*lam*(1-phi)*del2*sp_y; % (eq. 12)
y_w   = (wl1_y + wl2_y)^(1/nu -1)*(c_y)^(1/nu); %(eq. 13)
l1    = wl1_y*y_w;
l2    = wl2_y*y_w;
y_k2  = r/alph1*(k2_k1)^(-1); % (eq. 14)
s_y   = (y_k1)^(alph1/(1-alph1))*l1^(-alph2/(1-alph1)); % (eq. 14)
p     = sp_y*(s_y)^(-1); 
s_k2  = sp_k2/p; % (eq. 15)
k2    = (del2/(lam*eta)*s_k2*l2^(phi-1))^(1/(phi-1));
k1    = (k2_k1)^(-1)*k2;
y     = y_k2*k2;
c     = c_y*y;
s     = s_y*y;
w     = (y_w)^(-1)*y;
l     = l1+l2;
k     = k1+k2;
i     = del1*k;
it    = del2*s;

%Put the ss values in a vector consistent with Y and X vectors in model.m
yy  = [y c l l1 l2 r w i it p];
xx  = [k k1 k2 s];
ss  = [yy xx];
