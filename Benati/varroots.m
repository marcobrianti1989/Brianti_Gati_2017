function R=varroots(p,k,B)
%
%		Luca Benati
%		Bank of England
%		Monetary Assessment and Strategy Division
%		January 2003
%
% function R=varroots(p,k,B)
% This program computes the roots for the companion form of a VAR(p). 
%
%                                           Input of the program is
% p  = the order of the VAR
% k  = the number of variables in the VAR
% B = the VAR coefficients. They need to be stored as [MU; PHI(1); PHI(2); ...; PHI(p)], where MU is the vector of the
%     intercepts, and the PHI(k)'s are the autoregressive matrices.
% 
%                                         Output of the program is
%
% R = the roots of the VAR
%
[r,c]=size(B);
F=[B(:,2:c); eye((p-1)*k) zeros((p-1)*k,k)];
[V,R]=eig(F);
R=diag(R);