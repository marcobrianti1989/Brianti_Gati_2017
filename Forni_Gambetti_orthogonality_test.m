function [pvalue_news, pvalue_IT] = ...
      Forni_Gambetti_orthogonality_test(filename,...
      sheet,range,first_n_PCs,A,gamma,resid,which_shocks,nlags)
% We test if the structural shock of interest is orthogonal to the 
% the available information in the economy. 
% filename, sheet, and range are the name, the sheet, and 
% the coordinates of the xls file with the huge dataset.
% first_n_PCs: is the number which points out how many 
% principal components we want to use as a regressor to 
% test the orthogonality of the structural shock time series.
% A: is the standar Cholesky. (nvar,nvar)
% gamma: is the column vectors such that A*gamma is the impact 
% responses of the system to the shock of interst. (nvar,nshocks_recovered)
% resid: is the time series of reduced form residuals from the VAR
% (T,nvar).
T = size(resid,1);

% Read the data from a large dataset (with nvar variables)
[data, ~] = read_data2(filename, sheet, range);
if size(data,1) > T
      gap = size(data,1) - T;
      data = data(1+gap:end,:);
elseif size(data,1) < T
      error('Orthogonality Test cannot be performed. Some Time Series of the large dataset are too short.')
end

%Obtain nvar principal components
PCs = get_principal_components(data);

%Select the first first_n_PCs principal components
first_PCs = PCs(:,1:first_n_PCs);

%Obtain the structural shock time series
[structurals, ~] = get_structral_shocks_alternative(A,gamma,resid,which_shocks);

% Regress the structural shock time series on the first ...
% first_n_PCs principal components
structural_news          = structurals(which_shocks(1),:)';
structural_IT            = structurals(which_shocks(2),:)';
[~, residual_news]       = lagged_ols(structural_news,first_PCs,nlags);
[~, residual_IT]         = lagged_ols(structural_IT,first_PCs,nlags);

% Test if the structural shock is orthogonal to the ...
% first_n_PCs principal components
k                        = first_n_PCs;
q                        = first_n_PCs;
res_restricted_news      = structural_news;
res_unrestricted_news    = residual_news;
res_restricted_IT        = structural_IT;
res_unrestricted_IT      = residual_IT;
pvalue_news              = f_test(res_restricted_news, ...
      res_unrestricted_news, q, T, k);
pvalue_IT                = f_test(res_restricted_IT, ...
      res_unrestricted_IT, q, T, k);

end