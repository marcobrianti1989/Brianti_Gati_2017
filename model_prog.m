function [f fx fy fxp fyp eta R set dgam_dtheta deta_dtheta dR_dtheta xlag ylag] = model_prog(param, set)

param = param_trans(param);
%Assign parameter values to named variables.
rho = param(1);
phi = param(2);
rhoa = param(3);
alph = param(4);
sigz = param(5);
sige = param(6);

%Assign set values to named variables.
bet = set(1);
gam = set(2);
eta = set(3);
pis = set(4);
mu = set(5);
hbar = set(6);
rho_r = set(7);
rbar = set(8);
lam = set(9);
adiff = set(10);

%BEGIN_EXTRACT_HERE

%Use closed form expressions for the ss values.
w = (eta-1)/eta*gam;
h = ((eta-1)/(eta)*1/(1-rho/gam))^(mu/(mu+1));
c = gam*h;
x = c*(1-rho/gam);
R = pis*gam/bet;
d = c*(1-w/gam);
q = d*bet/(1-bet);
s = rho*c/gam;
rbar = R;


%Put the ss values in a vector consistent with X and Y vectors in model.m
Xss  = [gam s h c q d 1 R gam];
Yss  = [x   pis   w   h   q   R   d   c 1 gam gam gam];



%END_EXTRACT_HERE
%Compute Steady State
GAM= Xss(1);
S= Xss(2);
HL= Xss(3);
CL= Xss(4);
QL= Xss(5);
DL= Xss(6);
EPS= Xss(7);
RPL= Xss(8);
GAML= Xss(9);
XX= Yss(1);
PI= Yss(2);
W= Yss(3);
H= Yss(4);
Q= Yss(5);
RP= Yss(6);
D= Yss(7);
C= Yss(8);
GAMH= Yss(9);
GAMY= Yss(10);
GAMQ= Yss(11);
GAMD= Yss(12);
GAM_p= Xss(1);
S_p= Xss(2);
HL_p= Xss(3);
CL_p= Xss(4);
QL_p= Xss(5);
DL_p= Xss(6);
EPS_p= Xss(7);
RPL_p= Xss(8);
GAML_p= Xss(9);
XX_p= Yss(1);
PI_p= Yss(2);
W_p= Yss(3);
H_p= Yss(4);
Q_p= Yss(5);
RP_p= Yss(6);
D_p= Yss(7);
C_p= Yss(8);
GAMH_p= Yss(9);
GAMY_p= Yss(10);
GAMQ_p= Yss(11);
GAMD_p= Yss(12);

%Evaluate F.
f = [[XX - W/H^(1/mu), 1 - (RP*XX*bet)/(GAM*PI_p*XX_p), 1 - (XX*bet*(D_p + Q_p))/(Q*XX_p), S - C + XX, C - D - H*W, C*(eta*(W/GAM - 1) + 1) - PI*W*phi*(PI - pis) + (PI_p*W_p*XX*bet*phi*(PI_p - pis))/XX_p, C - GAM*H + (GAM*phi*(PI - pis)^2)/2, log(RP/rbar) - log(EPS) - rho_r*(log(RPL) - log(rbar)) + alph*(log(PI) - log(pis))*(rho_r - 1), log(GAM_p/gam) - rhoa*log(GAM/gam), S_p - (C*rho)/GAM, HL_p - H, CL_p - C, QL_p - Q, DL_p - D, GAMH - H/HL, GAMY - (C*GAML)/CL, GAMQ - (GAML*Q)/QL, GAMD - (D*GAML)/DL, EPS_p - 1, RPL_p - RP, GAML_p - GAM]];
%Evaluate derivative expressions.
fx = [[0, 0, 0, 0, 0, 0, 0, 0, 0]; [(RP*XX*bet)/(GAM*PI_p*XX_p), 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, S, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [-(C*W*eta)/GAM, 0, 0, 0, 0, 0, 0, 0, 0]; [(GAM*phi*(PI - pis)^2)/2 - GAM*H, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, -1, -rho_r, 0]; [-rhoa, 0, 0, 0, 0, 0, 0, 0, 0]; [(C*rho)/GAM, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, H/HL, 0, 0, 0, 0, 0, 0]; [0, 0, 0, (C*GAML)/CL, 0, 0, 0, 0, -(C*GAML)/CL]; [0, 0, 0, 0, (GAML*Q)/QL, 0, 0, 0, -(GAML*Q)/QL]; [0, 0, 0, 0, 0, (D*GAML)/DL, 0, 0, -(D*GAML)/DL]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [-GAM, 0, 0, 0, 0, 0, 0, 0, 0]];
fy = [[XX, 0, -W/H^(1/mu), (H*W)/(H^(1/mu + 1)*mu), 0, 0, 0, 0, 0, 0, 0, 0]; [-(RP*XX*bet)/(GAM*PI_p*XX_p), 0, 0, 0, 0, -(RP*XX*bet)/(GAM*PI_p*XX_p), 0, 0, 0, 0, 0, 0]; [-(XX*bet*(D_p + Q_p))/(Q*XX_p), 0, 0, 0, (XX*bet*(D_p + Q_p))/(Q*XX_p), 0, 0, 0, 0, 0, 0, 0]; [XX, 0, 0, 0, 0, 0, 0, -C, 0, 0, 0, 0]; [0, 0, -H*W, -H*W, 0, 0, -D, C, 0, 0, 0, 0]; [(PI_p*W_p*XX*bet*phi*(PI_p - pis))/XX_p, - PI^2*W*phi - PI*W*phi*(PI - pis), (C*W*eta)/GAM - PI*W*phi*(PI - pis), 0, 0, 0, 0, C*(eta*(W/GAM - 1) + 1), 0, 0, 0, 0]; [0, GAM*PI*phi*(PI - pis), 0, -GAM*H, 0, 0, 0, C, 0, 0, 0, 0]; [0, alph*(rho_r - 1), 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, -(C*rho)/GAM, 0, 0, 0, 0]; [0, 0, 0, -H, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, -C, 0, 0, 0, 0]; [0, 0, 0, 0, -Q, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, -D, 0, 0, 0, 0, 0]; [0, 0, 0, -H/HL, 0, 0, 0, 0, GAMH, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, -(C*GAML)/CL, 0, GAMY, 0, 0]; [0, 0, 0, 0, -(GAML*Q)/QL, 0, 0, 0, 0, 0, GAMQ, 0]; [0, 0, 0, 0, 0, 0, -(D*GAML)/DL, 0, 0, 0, 0, GAMD]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, -RP, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]];
fxp = [[0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [1, 0, 0, 0, 0, 0, 0, 0, 0]; [0, S_p, 0, 0, 0, 0, 0, 0, 0]; [0, 0, HL_p, 0, 0, 0, 0, 0, 0]; [0, 0, 0, CL_p, 0, 0, 0, 0, 0]; [0, 0, 0, 0, QL_p, 0, 0, 0, 0]; [0, 0, 0, 0, 0, DL_p, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, EPS_p, 0, 0]; [0, 0, 0, 0, 0, 0, 0, RPL_p, 0]; [0, 0, 0, 0, 0, 0, 0, 0, GAML_p]];
fyp = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [(RP*XX*bet)/(GAM*PI_p*XX_p), (RP*XX*bet)/(GAM*PI_p*XX_p), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [(XX*bet*(D_p + Q_p))/(Q*XX_p), 0, 0, 0, -(Q_p*XX*bet)/(Q*XX_p), 0, -(D_p*XX*bet)/(Q*XX_p), 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [-(PI_p*W_p*XX*bet*phi*(PI_p - pis))/XX_p, (PI_p^2*W_p*XX*bet*phi)/XX_p + (PI_p*W_p*XX*bet*phi*(PI_p - pis))/XX_p, (PI_p*W_p*XX*bet*phi*(PI_p - pis))/XX_p, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]];

eta = [[sigz, 0]; [0, 0]; [0, 0]; [0, 0]; [0, 0]; [0, 0]; [0, sige]; [0, 0]; [0, 0]];
R = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]];
