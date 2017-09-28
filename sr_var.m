function [A,B,res] = sr_var(dataset, nlags)

% Rearrange data to run a VAR, outputting A, the impact matrix identified using Cholesky SR restriction.
% res = residuals from VAR estimation

T = size(dataset,1); % time periods
nvar = size(dataset,2);

%Taking the lag on the dataset
for ilag = 1:nlags+1
      eval(['data_', num2str(ilag), ' = dataset(2+nlags-ilag:end+1-ilag,:);'])
end

%Building the matrix of regressor (which are the lags of the left hand side)
reg = zeros(T-nlags,nlags*nvar);
for ilag = 2:nlags+1
      for ivar = 1:nvar
            eval(['reg(:,ivar+size(dataset,2)*(ilag-2)) = data_',num2str(ilag),'(:,ivar);'])
      end
end

% Add constant
reg = [ones(T-nlags,1) reg];

%OLS in the Reduced Form Vector Autoregression
B = (reg'*reg)^(-1)*(reg'*data_1);

%Evaluate the residuals (res = yhat - y)
res = data_1 - reg*B;

%Compute the variance-covariance matrix
sigma = res'*res/(T-nlags);

%Static rotation matrix
A = chol(sigma)';

% Proving that A is what we are looking for A'u'uA = res*res where u'u = I
zero1 = A*A' - sigma;
zero2 = sum(sum(zero1.^2));
if zero2 > 10^(-16)
      error('The rotation matrix is not correct. Check the code')
end