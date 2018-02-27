% This file corresponds to main_file_VECM, it is a tool to test if the
% function VECM is working

% Marco Brianti, Laura Gáti, Feb 22 2018

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

%Run VECM

constant = 0; %no contant in the VECM
[Pi,B,res,sigma] = VECM(data, nlags, constant);
