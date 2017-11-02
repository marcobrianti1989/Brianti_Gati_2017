% This file corresponds to main_file_SVAR, except it performs Ryan's
% identification strategy to identify news and IT shocks

% Marco Brianti, Laura Gáti, Oct 28 2017

% TOOK 46.8612 MIN TO RUN!!

clear all
close all

tic
%Data Reading and Transformation
filename = 'dataset_23_sept_2017';
sheet    = 'Sheet1';
range    = 'B126:K283';
[data, shocknames,varnames, which_shocks, pos_rel_prices] = ...
    read_data(filename, sheet, range); %q position for relative prices

%Technical Parameters
max_lags        = 10;
nburn           = 0; %with the Kilian correction better not burning!!!
nsimul          = 500; %5000
nvar            = size(data,2);
sig             = 0.90; % significance level
H               = 100; %40; % horizon for generation of IRFs
h               = H; %40; % horizon for IRF plots
which_variable  = 1; % select TFP as the variable whose FEV we wanna max

%%Checking the number of lags over BIC, AIC, and HQ (see 'Lecture2M' in our folder)
[AIC,BIC,HQ] = aic_bic_hq(data,max_lags);
if AIC >= 4
    nlags = 1;
    warning(['AIC > 4, setting nlags to ', num2str(nlags) ])
else
    nlags = AIC;
end

%Run VAR imposing Cholesky
[A,B,res,sigma] = sr_var(data, nlags);

%Checking if the VAR is stationary
test_stationarity(B');

% TO DO: check info sufficiency

% dbstop in Ryan_two_stepsID at 51
% Implement Ryan's ID strategy
% [impact, FEV_opt, IRFs, gamma_opt, FEV_news, FEV_IT] = ryansID(which_variable,which_shocks,H,B,A,q);
[impact, FEV_opt, ~, gam3_opt, FEV_news, FEV_IT] ...
        = Ryan_two_stepsID(which_variable,which_shocks,H,B,A,pos_rel_prices);
% impact is the nvar x 2 impact matrix of news and IT.

% %Calculate IRFs, bootstrapped CI and plot them

% % Do bootstrap: TO DO: correct bootstrappedCIs function to...
% % 1.) recover impact matrix Q = A*gamma
% % 2.) Be general.
which_correction = 'blocks'; % [none, blocks] --> Choose whether to draws residuals in blocks or not.
blocksize = 5; % size of block for drawing in blocks
[A_boot, B_boot, fake_impact_boot] = ...
    bootstrappedCIs(B, nburn, res, nsimul, which_correction, blocksize, nvar, ...
    nlags,pos_rel_prices, which_shocks,which_variable,H); % automatically does Kilian correction

fake_impact = zeros(nvar,nvar);
fake_impact(:,which_shocks) = impact;

[IRFs, ub, lb] = genIRFs(fake_impact,fake_impact_boot,B,B_boot,H,sig);
plotIRFs(IRFs,ub,lb,h,which_shocks,shocknames,varnames)

toc

