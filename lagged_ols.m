function [B, res] = lagged_ols(y,x,nlags)
% This function run OLS: y = Bx + res where both y and x are timeseries and
% x can be multidimensional vectors, y is a single variable. This function regress the
% past of x on y using nlags. x and y should be inputed at time t and the
% code is going to properly lag x.

T = size(x,1); 
T1 = size(y,1);
if abs(T - T1) > 0
      error('x and y have two different time dimensions.')
end
nvar = size(x,2); %number of regressors

%Taking the lag on the dataset
for ilag = 1:nlags+1
    eval(['reg_first_', num2str(ilag), ' = x(2+nlags-ilag:end+1-ilag,:);'])
end

%Building the matrix of regressor (which are the lags of the left hand side)
reg = zeros(T-nlags,nlags*nvar);
for ilag = 2:nlags+1
    for ivar = 1:nvar
        eval(['reg(:,ivar+size(dataset,2)*(ilag-2)) = reg_first_',num2str(ilag),'(:,ivar);'])
    end
end

% Add constant
reg = [ones(T-nlags,1) reg]; % (T-nlags, 1+nvar*nlags)

% Readjusting y series given the number of lags
dependent_variable = y(nlags+1:end);

%OLS in the Reduced Form Vector Autoregression
B = (reg'*reg)\(reg'*dependent_variable);

%Evaluate the residuals (res = yhat - y)
res = dependent_variable - reg*B;

%Compute the variance-covariance matrix
% sigma = res'*res/(T-nlags);












end