function [AIC,BIC,HQ] = aic_bic_hq(dataset,max_lags) 

 TO BE FINISHED

%These criteria minimize the following objective functions:
%AIC(p) = T ln|Omega_hat| + 2(n^2*p)
%BIC(p) = T ln|Omega_hat| + (n^2*p)lnT
%HQ(p)  = T ln|Omega_hat| + 2(n^2*p)lnlnT
%where n is the number of variables, T is the length of the dataset, and p
%is the number of lags

T = size(dataset,1);
nvar = size(dataset,2);
res_test = zeros(T,nvar,max_lags);
for nlags = 1:max_lags
      [A,B,res] = sr_var(dataset,nlags);
      AIC(nlags) = sum(sum(res));
end
BIC = 0;
HQ = 0;
end