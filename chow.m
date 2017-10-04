function [reject, pval] = chow(series, tau, max_lags, time_vector)
% Chow test for structural breaks at a known point
% series = single series on which to test for structural break 
% tau = candidate time period for break (specify as an time index in data, i.e. not year, but e.g. 50th observation)
% max_lags = maximum number of lags for lag selection
% time_vector = the datenumbers of time periods, of which tau is the period
% we're interested in

candidate = time_vector(tau);

% Split sample into (before and incl. tau) and (after tau)
series1 = series(1:tau);
N1 = length(series1);
series2 = series(tau+1:end);
N2 = length(series2);

% Estimate an AR(k)
% Whole sample
% Select a lag order for AR:
[AIC,BIC,HQ] = aic_bic_hq(series,max_lags);
nlags = AIC;
[~,~,res,~] = sr_var(series, nlags);
RSS_p = sum(res.^2);

% Sample up to and including tau
% Select a lag order for AR:
[~,~,res,~] = sr_var(series1, nlags);
RSS_1 = sum(res.^2);

% Sample up after tau
% Select a lag order for AR:
[~,~,res,~] = sr_var(series2, nlags);
RSS_2 = sum(res.^2);

% Create Chow F-stat
k = AIC;
F = ((RSS_p -(RSS_1 + RSS_2))/k)/((RSS_1 + RSS_2)/(N1 + N2 - 2*k));
v1 = k;
v2 = (N1 + N2 - 2*k);
pval = fcdf(F,v1,v2);
if pval >= 0.9
      reject = 1;
      disp(['H0 of no structural break at ', datestr(candidate, 'yyyy'), ' is rejected at ' ...
            num2str(100*pval), '%.']);
else
      reject = 0;
      disp(['H0 of no structural break at ', datestr(candidate, 'yyyy'), ' is NOT rejected.']);
end


