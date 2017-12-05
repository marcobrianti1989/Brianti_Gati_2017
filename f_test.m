function pvalue = f_test(res_r, res_ur, q, n, k)
% Inputs
% res_r  = residuals from restricted regression
% res_ur = residuals from unrestricted regression
% q      = no. of restrictions
% k      = no. of regressors
% n      = length of dataset (e.g time dimension)
% Outputs
% sig_level = the significance at which the null would be rejected

SSR_R  = sum(res_r.^2);
SSR_UR = sum(res_ur.^2);

numerator   = (SSR_R - SSR_UR)/q;
denominator = SSR_UR / (n-k-1);

F = numerator / denominator;
sig_level = fcdf(F, q, n-k-1);
pvalue = 1 - sig_level;

end