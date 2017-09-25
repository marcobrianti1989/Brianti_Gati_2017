% QUICK_IR - Compute impulse response of a VAR given impact response x0. The 
% VAR model is written 
%  
%       y(t) = [y(t-1) y(t-2) ... y(t-nlag)]*beta + mu(t)
%
% usage
%
% out = quick_IR(beta, nt, x0)
%
% where
% 
% beta = the VAR coefficient
% nt   = the number of period to compute the impulse response
% x0   = the impact response of the shock (a.k.a. initial state of the
%        system)


function [out,F] = quick_IR(beta, nt, x0)

nvar = size(beta,2);
nlag = size(beta,1)/nvar;

%***********************
% Get Companion Matrix
%***********************
F = zeros(nvar*nlag);
F(1:nvar,:) = beta';
for j =1:nlag-1
    F(j*nvar+1:(j+1)*nvar,(j-1)*nvar+1:(j)*nvar) = eye(nvar);
end

%***************************
%Compute IR Quickly
%***************************
out = zeros(nvar*nlag, nt);
out(1:nvar,1) = x0';
for j = 1:nt-1
    out(:,j+1) = F*out(:,j);
end
out = out(1:nvar,:)';
