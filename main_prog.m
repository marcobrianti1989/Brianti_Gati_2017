% ====================SIMPLIFIED VERSION OF COMIN AND GERTLER (2006)=======================
%                       WITH NOISY R&D EFFICIENCY
% =================== 28 of June, 2017.                 ==================

%Code provided by Marco Brianti & Laura Gati on the basis of
%the Stephanie Schmitt-Grohé and Martín Uribe algorithm

clear all
close all

%Parameterization
beta = 0.98; %discount factor
delta = 0.025; %Depreciation rate of capital - set equal to 1 to check if the model is correct: I = K!
alpha = 1/3; %capital share - Cobb-Douglass production function
sigma = 0.01; %standard deviation of the TFP shock
phi = 0.9; %rate of survival of the old technology - obsolescense parameter
rhog = 0.8; %persistence of the exogenous TFP shock
theta = 0.5; %share between standard production and technology 
fit = 0.7; %Efficiency of the R&D sector <--- LG: this guy should become a variable for incomplete info case.

ss = 0.1*ones(1,11); %vector of all the steady state variables <--- LG: Need to calculate steady state properly
xx = ones(3);
yy = ones(8);

%% 1.) Model with full information

% Compute [fyn, fxn, fypn, fxpn] 
model_fullinfo;

[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);

save('gxhx.mat', 'gx', 'hx')

%Eigenvalues of hx
disp('Computing eigenvalues of hx');
disp(eig(hx))

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