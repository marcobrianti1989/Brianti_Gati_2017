% QUICK_VAR - Estimate a VAR model using OLS. The VAR model is written as
%
%         y(t) = [y(t-1) y(t-2) ... y(t-nlag)]*beta + mu(t)
%
% usage
%
% [beta, c, mu] = quick_var(Y,nlag)
%
% where
%
% Y = a T-by-ny column vector of observables
% nlag = the number of lags to include in the VAR
%
% and
%
% beta  = (ny*nlag)-by-ny coefficient of VAR (i.e. VAR model written as
%         y(t) = [y(t-1) y(t-2)...y(t-nlag)]*beta + mu(t))\
% c     = the choleski decomposition of the covariance matrix of residuals
% mu    = T-by-ny matrix of reduced form residuals
% omega = the covariance matrix of the residuals
%
% NOTES: By default, quick VAR includes a constant term in the VAR. The
% constant is eliminated in the output argument beta.

function [beta, c, mu, omega] = quick_var(Y,nlag)

%Number of Variables
nvar = size(Y,2);
T = length(Y)-nlag;
const = 1;

%*****************
% Arrange Data
%*****************

%Stack lagged X
X = zeros(T, nlag*nvar);
for j =1:nlag
    X(:, (j-1)*nvar+1:j*nvar) = Y(nlag+1-j:end-j,:);
end

%Include a constant in VAR
X = [X, ones(length(X),const)];

%Exclude lags from DY
YY = Y(nlag+1:end,:);

%*********************************************
% Regress, get covariance matrix and cholesky
%*********************************************

%VAR coefficients, stderrs, etc
beta = (X'*X)\(X'*YY);
mu = YY - X*beta;
omega = mu'*mu/(T);
c = chol(omega)';
beta = beta(1:end-const,:);


