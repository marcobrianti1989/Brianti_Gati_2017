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
 
% This is OULTON 2010 WITH GROWTH

function [fyn, fxn, fypn, fxpn] = model_spillover_news(param)

%Steady State
[ss, param] = model_spillover_ss_news(param);

%Declare parameters
bet      = param.bet;  %Never call anything beta...it's a matlab function
a        = param.a; % share of standard capital
b        = param.b; % share of IT capital
biggami  = param.biggami; % growth rate of equipment productivity
biggamc  = param.biggamc;    % growth rate of neutral technology
di       = param.di; % depreciation rate of capital IT
dc       = param.dc; % depreciation rate of capital standard
chi      = param.chi; % preference parameter.
gam      = param.gam; % spillover elasticity

%Declare Needed Symbols
syms KC KI KC1 KC2 KI1 KI2 % 6 vars
syms KC_p KI_p KC1_p KC2_p KI1_p KI2_p 
syms BIGGAMC BIGGAMI % 2 vars
syms BIGGAMC_p BIGGAMI_p 
syms YC YI C IC IT W RC RI H H1 H2 % 11 vars
syms YC_p YI_p C_p IC_p IT_p W_p RC_p RI_p H_p H1_p H2_p  
syms P P_p % 1 var  % -> 20 vars
syms CL CL_p KIL KIL_p KCL KCL_p YCL YCL_p YIL YIL_p HL HL_p PL PL_p KC2L KC2L_p KI2L KI2L_p
syms GAMC GAMC_p GAMKI GAMKI_p GAMKC GAMKC_p GAMYC GAMYC_p GAMYI GAMYI_p GAMH GAMH_p GAMP GAMP_p GAMKC2 GAMKC2_p GAMKI2 GAMKI2_p
syms BIGGAMCL BIGGAMCL_p BIGGAMIL BIGGAMIL_p
syms EXPGC EXPGC_p EXPGI EXPGI_p
syms V0 V1 V2 V3 V4 V5 V6 V7 V8 V0_p V1_p V2_p V3_p V4_p V5_p V6_p V7_p V8_p
syms SIT SIT_p RIS RIS_p
syms N N_p

% %Declare X and Y vectors
X  = [KC KI BIGGAMC BIGGAMI ...
    CL BIGGAMCL BIGGAMIL YCL YIL HL PL KC2L KI2L KIL KCL ...
    V0 V1 V2 V3 V4 V5 V6 V7 V8 ...
    SIT ...
    N]; % vector of state variables  
XP = [KC_p KI_p BIGGAMC_p BIGGAMI_p ...
    CL_p BIGGAMCL_p BIGGAMIL_p YCL_p YIL_p HL_p PL_p KC2L_p KI2L_p KIL_p KCL_p ...
    V0_p V1_p V2_p V3_p V4_p V5_p V6_p V7_p V8_p ...
    SIT_p ...
    N_p]; % p signifies t+1 

Y  = [YC YI C IC IT W RC RI H H1 H2 KC1 KC2 KI1 KI2 P EXPGC EXPGI ...
    GAMC GAMKI GAMYC GAMYI GAMH GAMP GAMKC2 GAMKI2 RIS]; % vector of controls
YP = [YC_p YI_p C_p IC_p IT_p W_p RC_p RI_p H_p H1_p H2_p KC1_p KC2_p KI1_p KI2_p P_p EXPGC_p EXPGI_p ...
    GAMC_p GAMKI_p GAMYC_p GAMYI_p GAMH_p GAMP_p GAMKC2_p GAMKI2_p RIS_p] ;

%Make index variables for future use
make_index([Y,X])

% Model Equations 
f(1)     = -YC + KI^gam * N *(BIGGAMC)*H1^(1-a-b)*KC1^(a)*KI1^(b); 
f(end+1) = -YI + KI^gam * N * (BIGGAMI)*H2^(1-a-b)*KC2^(a)*KI2^(b); 
f(end+1) = -KC + KC1 + KC2;
f(end+1) = -KI + KI1 + KI2; %
f(end+1) = -H + H1 + H2;
f(end+1) = -EXPGC + BIGGAMC^((1-b-gam)/(1-a-b-gam))*BIGGAMI^((b+gam)/(1-a-b-gam));
f(end+1) = -EXPGI + BIGGAMC^(a/(1-a-b-gam))*BIGGAMI^((1-a)/(1-a-b-gam));
f(end+1) = -KC_p*(EXPGC) + (1-dc)*KC + IC;
f(end+1) = -KI_p*(EXPGI)+ (1-di)*KI + IT; 
f(end+1) = -YC + C + IC;
f(end+1) = -YI + IT;
% f(end+1) = -C + W/chi; % case of linear V(H) = -chi*H
% f(end+1) = -C -(1-H)*W/chi; Case of V(H) = chi*log(1-H)
f(end+1) = -C + W/(chi*H); % case of V(H) = chi/2 * H^2
f(end+1) = -1 + bet*C/C_p*1/EXPGC*(RC_p + 1-dc);
f(end+1) = -1 + bet*C/C_p*1/EXPGI*(SIT_p*RI_p/P_p + 1-di);
f(end+1) = -W + KI^gam * (1-a-b)*(BIGGAMC)*H1^(-a-b)*KC1^(a)*KI1^(b); %
f(end+1) = -RC + KI^gam * a*(BIGGAMC)*H1^(1-a-b)*KC1^(a-1)*KI1^(b); %
f(end+1) = -RI*SIT + KI^gam * b*(BIGGAMC)*H1^(1-a-b)*KC1^(a)*KI1^(b-1); %
f(end+1) = -W + KI^gam*(1-a-b)*(BIGGAMI)*H2^(-a-b)*KC2^(a)*KI2^(b)*P;%
f(end+1) = -RC + KI^gam*a*(BIGGAMI)*H2^(1-a-b)*KC2^(a-1)*KI2^(b)*P;%
f(end+1) = -RI*SIT + KI^gam*b*(BIGGAMI)*H2^(1-a-b)*KC2^(a)*KI2^(b-1)*P;%
f(end+1) = log(BIGGAMC_p/biggamc) - .8*log(BIGGAMC/biggamc); %taken directly from Ryan's example code.
f(end+1) = log(BIGGAMI_p/biggami) - .8*log(BIGGAMI/biggami); %taken directly from Ryan's example code.
f(end+1) = log(N_p) - 0.8*log(N) -V0; % common component of technology on which there is a news shock
% Approach for growth rates:
f(end+1) = BIGGAMCL_p - BIGGAMC;
f(end+1) = BIGGAMIL_p - BIGGAMI;
f(end+1) = CL_p - C;
f(end+1) = KIL_p - KI;
f(end+1) = KCL_p - KC;
f(end+1) = YCL_p - YC;
f(end+1) = YIL_p - YI;
f(end+1) = HL_p - H;
f(end+1) = PL_p -P;
f(end+1) = KC2L_p - KC2;
f(end+1) = KI2L_p - KI2;
f(end+1) = GAMC - C/CL*(BIGGAMCL^((1-b-gam)/(1-a-b-gam)) *BIGGAMIL^((b+gam)/(1-a-b-gam)));
f(end+1) = GAMKI - KI/KIL*(BIGGAMCL^((a)/(1-a-b-gam)) *BIGGAMIL^((1-a)/(1-a-b-gam)));
f(end+1) = GAMYC - YC/YCL*(BIGGAMCL^((1-b-gam)/(1-a-b-gam)) *BIGGAMIL^((b+gam)/(1-a-b-gam)));
f(end+1) = GAMYI - YI/YIL*(BIGGAMCL^((a)/(1-a-b-gam)) *BIGGAMIL^((1-a)/(1-a-b-gam)));
f(end+1) = GAMH - H/HL;
f(end+1) = GAMP - P/PL*BIGGAMCL/BIGGAMIL;
f(end+1) = GAMKC2 - KC2/KC2L*(BIGGAMCL^((1-b-gam)/(1-a-b-gam)) *BIGGAMIL^((b+gam)/(1-a-b-gam)));
f(end+1) = GAMKI2 - KI2/KI2L*(BIGGAMCL^((a)/(1-a-b-gam)) *BIGGAMIL^((1-a)/(1-a-b-gam)));
% Approach for news shocks:
f(end+1) = V8_p;
f(end+1) = V0_p - V1;
f(end+1) = V1_p - V2;
f(end+1) = V2_p - V3;
f(end+1) = V3_p - V4;
f(end+1) = V4_p - V5;
f(end+1) = V5_p - V6;
f(end+1) = V6_p - V7;
f(end+1) = V7_p - V8;
f(end+1) = log(SIT_p) - .8*log(SIT);
f(end+1) = SIT*RI - RIS;



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
