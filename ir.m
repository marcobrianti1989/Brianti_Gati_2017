function [IR, iry, irx]=ir(gx,hx,x0,T)
%IR=ir(gx,hx,x0,T) computes T-period 
%impulse responses of the 
%vector [y x], whose law of 
%motion is:
% x(t+1) = hx x(t)
% y(t) = gx x2(t)
% with initial condition x0
% Inputs: gx, hx, x0, and 
%T (optional, default value 10)
%(c) Stephanie Schmitt-Grohe and Martin Uribe, August 18, 1993. 

if nargin < 4
T =10;
end 

x0=x0(:);
pd=length(x0);
MX=[gx;eye(pd)]; % MX is a vector of [gx(:) hx(:)] which in each iteration gets multiplied by the impulse of that period.
IR=[]; % (T, nvar) with the jumps first and then the states.
x=x0;
for t=1:T
IR(t,:)=(MX*x)';
x = hx * x;
end
ny = size(gx,1);
iry = IR(:,1:ny);
irx = IR(:,ny+1:end); 
