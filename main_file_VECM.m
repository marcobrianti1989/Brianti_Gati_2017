% This file corresponds to main_file_VECM, it is a tool to test if the
% function VECM is working

% Marco Brianti, Laura Gáti, Mar 13 2018
clear
close all

current_dir = pwd;
cd ../.. % go up 2 levels
base_path = pwd;
cd(current_dir)
addpath(base_path)

if exist([base_path '\Data'], 'dir')
      addpath([base_path '\Data']) %for Microsoft
else
      addpath([base_path '/Data']) %for Mac
end

%Data Reading and Transformation
filename = 'dataset_main';
sheet    = 'Data';
range    = 'B1:W287';
[data, varnames] = read_data2(filename, sheet, range);
shocknames = {'News Shock','IT Shock'};
varnames

%Technical Parameters
max_lags        = 10;
nvar            = size(data,2);

%%Checking the number of lags over BIC, AIC, and HQ (see 'Lecture2M' in our folder)
[AIC,BIC,HQ] = aic_bic_hq(data,max_lags);
if AIC >= 4
    nlags = 2; % 2;
    warning(['AIC > 4, setting nlags to ', num2str(nlags) ])
else
    nlags = AIC;
end

%Run reduced form VECM
r = 1; %number of cointegrating vectors
%[alph_hat2, bet_hat2,Gam_hat2,sigma2] = RRR(data, nlags, r);
[alph_hat,bet_hat,Pi,Gam_hat,res,sigma] = redu_VECM(data, nlags, r); %It seems to be correct...
%compared to a different method.

%Run structural VECM
diff_BBp_sigma = 1;
j = 0;
while diff_BBp_sigma >= 10^(-6)
    j = j + 1
    [B, Xi, A, diff_BBp_sigma] = structural_VECM(alph_hat,bet_hat,Gam_hat,res,sigma,nlags,r);
    diff_BBp_sigma
end

% Implement the "just IT" ID strategy in a VAR
which_variable = 3;
which_shock = 1;
[impact, impact_IT_opt, gam_opt]  = just_IT_ID(which_variable,which_shock,B);

return
close all
H = 40;
sig1 = 0;
sig2 = 0;
[IRFs, ub1, lb1, ub2, lb2] = genIRFs_VECM(B,0,A',0,H,sig1,sig2);
print_figs = 'no';
which_shock = [1 2 3 4];
simple_plotIRFs(IRFs,ub1,lb1,H,which_shock, print_figs)

