% This file corresponds to main_file_VECM, it is a tool to test if the
% function VECM is working

% Marco Brianti, Laura Gáti, Mar 13 2018

clear
close all

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
    nlags = 1; % 2;
    warning(['AIC > 4, setting nlags to ', num2str(nlags) ])
else
    nlags = AIC;
end

%Run reduced form VECM
r = nvar-1; %number of cointegrating vectors
%[alph_hat, bet_hat,Gam_hat,sigma] = RRR(data, nlags, r);
[alph_hat,bet_hat,Pi,Gam_hat,res,sigma] = redu_VECM(data, nlags, r);

%Run structural VECM
[B, Xi, A] = structural_VECM(alph_hat,bet_hat,Gam_hat,res,sigma,nlags,r);

return

H = 40;
sig1 = 0;
sig2 = 0;
[IRFs, ub1, lb1, ub2, lb2] = genIRFs_VECM(B,0,A',0,H,sig1,sig2);
print_figs = 'no';
which_shock = [1,2,3];
simple_plotIRFs(IRFs,ub1,lb1,H,which_shock, print_figs)

