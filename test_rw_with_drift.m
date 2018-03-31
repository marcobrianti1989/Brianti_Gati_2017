function [a, b, a_counter, b_counter] = test_rw_with_drift(drift, rho, T)

% Inputs
% drift = a value for the drift
% rho   = AR-coefficient
% T     = length of process
% The RW with drift process is log(A'/drift) = rho*log(A/drift) + e'  (' stands for t+1)
% where A = b'/b, so A is the by-construction stationarized level of tech,
% and b is the actual level of the tech process (TFP level).

% Outputs:
% a = series of stationarized TFP levels with shock (= gross growth rate of b)
% b = series of actual TFP levels with shock
% a_counter = counterfactual a w/o a shock
% b_counter = counterfactual b w/o a shock

a0  = drift; % we start in steady state, so initial value = drift
bc0 = 1;
a   = zeros(T,1);
b   = zeros(T,1);
b(1) = a0*bc0;
a(1)= drift* exp(log((a0/drift)^rho) + 1); % a shock of 1
a_counter   = zeros(T,1);
b_counter   = zeros(T,1);
a_counter(1)= drift* exp(log((a0/drift)^rho)); % no shock
b_counter(1)= b(1);

for i=2:T
    a(i) = drift* (a(i-1)/drift)^rho;
    b(i) = a(i-1)*b(i-1);
    a_counter(i) = drift* (a_counter(i-1)/drift)^rho;
    b_counter(i) = a_counter(i-1)*b_counter(i-1);
end

b= log(b);
b_counter = log(b_counter);