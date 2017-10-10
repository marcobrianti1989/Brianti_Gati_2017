function [B, yhat, res] = quick_ols(Y,X)

B    = (X'*X)\(X'*Y);
yhat = X*B;
res = Y - yhat;
