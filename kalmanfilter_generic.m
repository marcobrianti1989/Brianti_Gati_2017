function [P,K,Om,numiter] = kalmanfilter_generic(hx,gx, X0, P0, R, Q, tol)
% Code relies on "Kalman Filter Once and For all"
% General variable names:
% X  = vector of states
% Y  = vector of jumps
% Sx = selector matrix: which of the states is being formed beliefs upon (input here as Sx*hx already)
% Sy = selector matrix: which of the jumps is being seen and used for
% inference on noisy state (input here as Sy*gx already)
% Inputs:
% hx = state transition matrix
% gx = jump dependence matrix on states
% X0 = initial belief on state = X(t|t-1). Usually initialized at E(X).
% P0 = initial FEV of state = Var(X(t)). Usually initialized at Var(X).
% R  = Cov(Wt), where W = measurement error of observation equation.
% Q  = Cov(Epsilont), where Epsilon = error in state equation
% tol = convergence error tolerance
% Outputs:
% P  = forecast error variance of belief on state
% K  = Kalman gain
% Om = forecast error variance of jump (Omega from class)
% numiter = number of iterations until convergence

Yt = 4; %data arriving, can be anything. Check on this!
%% Initialize
Xtt_1 = X0;
Ptt_1 = P0;

err = 10;
numiter = 0;
while err > tol
%% Step 1. Forecast Y(t|t-1) using X(t|t-1)
Ytt_1  = gx*Xtt_1;
Omtt_1 = gx*Ptt_1*gx' + R;

%% Step 2. Infer state at t using info arriving at time t. (i.e. get X(t|t))
K   = Ptt_1*gx'*Omtt_1^(-1);
Xtt = Xtt_1 + K*(Yt - Ytt_1);
Ptt = Ptt_1 - Ptt_1*gx'*Omtt_1^(-1)*gx*Ptt_1;

%% Forecast states at t+1 using t info 
Xt1t = hx*Xtt;
Pt1t = hx*Ptt*hx' + Q;

%% Create error and update X(t|t-1) and P(t|t-1)
err = Pt1t - Ptt;
Xtt_1 = Xt1t;
Ptt_1 = Pt1t;
numiter = numiter+1;
end

% Once P has converged, set its steady state value, as well as that of Om.
% K is simply given as the last value in the iteration
P = Ptt_1;
Om = Omtt_1;

