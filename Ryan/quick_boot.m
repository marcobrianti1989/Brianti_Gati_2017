% QUICK_BOOT - Generate a single bootstrapped timeseries
% using the procedure described in class.
%
% usage
%
% ysample = quick_boot(beta,mu,T,nburn,y0)
%
% where
%
% beta  = (ny*nlag)-by-ny coefficient of VAR (i.e. VAR model written as
%         y(t) = [y(t-1) y(t-2)...y(t-nlag)]*beta + mu(t))
% mu    = T-by-ny matrix of reduced form residuals
% T     = Number of periods to simulate (excluding burn in)
% nburn = Burn-in period
% y0    = intial values of y (nlag-by-ny matrix)


function ysamp = quick_boot(beta,mu,T,nburn,y0)

%Beta should be (ny*nlag)-by-ny;

%How many variables
ny = size(beta,2);

%How many lags
nlag = size(beta,1)/ny;

%Initialize output matrix
ysamp = zeros(T+nburn,ny);

%If fixing intial values
if nargin > 4
    ysamp(1:nlag,:) = y0;
end

%Draw mu's with replacement
musamp = replace(mu,(T+nburn-1));

%Generate new data
for j = nlag:(T+nburn-1)
    %My trick for stacking the data appropriately
    tmp = ysamp(j:-1:j+1-nlag,:)';tmp = tmp(:)';
    
    %Simulate new data
    ysamp(j+1,:) = tmp*beta + musamp(j-nlag+1,:);
    
end

ysamp = ysamp(nburn+1:end,:);


% REPLACE - This function takes samples with replacement
% Usage:
% [s, index] = replace (data, n)
%
% Inputs:
% data = data set, with variables in cols, obs in rows
% n = the number of observations to return
%
% Ouput:
% s = the sampled data
% index = the indexes of the sampled data

function [s, index] = replace(data, n)

%If simulated sample is of same length as data
if nargin < 2
    n = size(data,1);
end

%Select indexes
index = 1+floor((size(data,1))*rand(n,1));
s = data(index,:);
