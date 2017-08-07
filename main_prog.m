% ====================SIMPLIFIED VERSION OF COMIN AND GERTLER (2006)=======================
%                       WITH NOISY R&D EFFICIENCY
% =================== 28 of June, 2017.                 ==================

%Code provided by Marco Brianti & Laura Gati on the basis of
%the Stephanie Schmitt-Grohé and Martín Uribe algorithm

clear all
close all

%Parameterization
param

%% 1.) Model with full information

% Compute [fyn, fxn, fypn, fxpn] 
check_steadystate = 1; % 1 if we want to verify that our steady state is analytically correct
[fyn,fxn,fypn,fxpn] =model_1(check_steadystate); % model_fullinfo;

[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);

save('gxhx.mat', 'gx', 'hx')

%Eigenvalues of hx
disp('Computing eigenvalues of hx');
disp(eig(hx))

return
%% 2.) Model with incomplete information

tol = 1.0e-06;
R   = 0; % measurement error in state equation
eta = [0, 0, sigma]'; % vector of shocks
Q   = eta*eta'; % covariance matrix of shocks
X0  = [1 2 3]'; % initial E(X). For X0 and P0 I'm putting bullshit right now, to change later!
P0  = X0*X0'; % initial Var(X)

[P,K,Om, numiter] = kalmanfilter_generic(hx,gx, X0, P0, R, Q, tol); %it runs but is not 100% yet.

% Compute [fyn, fxn, fypn, fxpn] again
% model_incinfo; % to be written yet.

% % Compute model solution for incomplete information (ii) case
% [gx_ii,hx_ii]=gx_hx_alt(fyn,fxn,fypn,fxpn);
% 
% save('gxhx_ii.mat', 'gx_ii', 'hx_ii')
% 
% %Eigenvalues of hx
% disp('Computing eigenvalues of hx_ii');
% disp(eig(hx_ii))