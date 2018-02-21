%**************************************************************
% MAIN_PROG - Solves the neoclassical model with a random walk 
% TFP solved in class in EC860.
%
% Code by Ryan Chahrour, Boston College, 2014
%**************************************************************


clear all

%Load Parameters
param = parameters;

%Compute the first-order coefficiencients of the model
[fyn, fxn, fypn, fxpn] = model(param);

%Compute the transition and policy functions, using code by
%Stephanie Schmitt-Grohé and Martín Uribe (and available on their wedsite.)
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);

save('gxhx.mat', 'gx', 'hx')

%Eigenvalues of hx
disp('Computing eigenvalues of hx');
disp(eig(hx))

