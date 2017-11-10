% This file corresponds to main_file_SVAR, except it performs our
% identification strategy to identify news and IT shocks. Namely, both max
% FEV of TFP but we set short-run impact restrictions in order to
% distinguish the two.

% Marco Brianti, Laura Gáti, Oct 28 2017

clear all
close all

tic
%Data Reading and Transformation
filename = 'dataset_23_sept_2017';
sheet    = 'Sheet1';
range    = 'B126:K283';
[data, shocknames,varnames, which_shocks, pos_rel_prices] = ...
    read_data(filename, sheet, range); %q (pos_rel_prices) position for relative prices

return

%Technical Parameters
max_lags        = 10;
nburn           = 0; %with the Kilian correction better not burning!!!
nsimul          = 500; %5000
nvar            = size(data,2);
sig             = 0.95; % significance level
H               = 40; %40; % horizon for generation of IRFs
h               = H; %40; % horizon for IRF plots
which_variable  = 1; % select TFP as the variable whose FEV we wanna max

%%Checking the number of lags over BIC, AIC, and HQ (see 'Lecture2M' in our folder)
[AIC,BIC,HQ] = aic_bic_hq(data,max_lags);
if AIC >= 4
    nlags = 2;
    warning(['AIC > 4, setting nlags to ', num2str(nlags) ])
else
    nlags = AIC;
end

%Run VAR imposing Cholesky
[A,B,res,sigma] = sr_var(data, nlags);

%Checking if the VAR is stationary
test_stationarity(B');

% TO DO: check info sufficiency

% Implement our strategy
q = 3; %Position of Michigan Index
[impact, FEV_opt, ~, gamma_opt, FEV_news, FEV_IT] = ...
    FEVmax_sr_ID(which_variable,which_shocks,H,B,A,q);
%Nesting gamma inside a fake matrix of zeros to have a square matrix.
%Similar to Barsky&uncleSims this is a partial ID strategy
fake_impact                    = zeros(nvar,nvar);
fake_impact(:,which_shocks)    = impact;

% %Calculate IRFs, bootstrapped CI and plot them
which_correction    = 'blocks'; %[none, blocks] --> Choose whether to draws residuals in blocks or not.
which_ID            = 'FEVmax_sr'; %choose between different identification strategy
blocksize           = 5; % size of block for drawing in blocks
[A_boot, B_boot] = bootstrappedCIs(B, nburn, res, nsimul, which_correction, blocksize, nvar, ...
    nlags,pos_rel_prices, which_shocks,which_variable,H,which_ID);

[IRFs, ub, lb] = genIRFs(fake_impact,A_boot,B,B_boot,H,sig);
plotIRFs(IRFs,ub,lb,h,which_shocks,shocknames,varnames)

% Show the variance decomposition share
disp(['FEV share of News Shock on TFP is ' num2str(FEV_news)])
disp(['FEV share of IT Shock on TFP is ' num2str(FEV_IT)])


toc

