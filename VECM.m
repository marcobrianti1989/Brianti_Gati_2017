function [Pi,reg,res,sigma] = VECM(dataset, nlags)

% Rearrange data to run a VECM, outputting Pi, the correction matrix.
% reg the betas of the lagged terms. Similar to Reduced form VAR
% res = residuals from VAR estimation
% sigma the variance covariance matrix.

T = size(dataset,1); % time periods
nvar = size(dataset,2);
diff_dataset = diff(dataset);

%Taking the lag on the dataset
for ilag = 1:nlags+1
      eval(['data_', num2str(ilag), ' = diff_dataset(2+nlags-ilag:end+1-ilag,:);'])
end

%Building the matrix of regressor (which are the lags of the left hand side)
reg = zeros(T-nlags-1,nlags*nvar);
for ilag = 2:nlags+1
      for ivar = 1:nvar
            eval(['reg(:,ivar+size(diff_dataset,2)*(ilag-2)) = data_',num2str(ilag),'(:,ivar);'])
      end
end

% Add constant and correction term
reg = [ones(T-nlags-1,1) dataset(3:end,:) reg];

%OLS in the Reduced Form Vector Autoregression
B = (reg'*reg)\(reg'*data_1);

%Evaluate the residuals (res = yhat - y)
res = data_1 - reg*B;

%Compute the variance-covariance matrix
sigma = res'*res/(T-nlags);

%Compute the correction matrix
Pi = B(2:nvar+1,:);
B = B(nvar+2:end,:);

end