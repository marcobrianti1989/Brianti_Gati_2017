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
    nlags = 2; % 2;
    warning(['AIC > 4, setting nlags to ', num2str(nlags) ])
else
    nlags = AIC;
end

[LR_MET LR_TT] = Johansen_Test(data, nlags);