function hd = historical_decomposition(s,A_IRF,B, shock_number, variable_number)
% Based on Kilian "SVAR Chapter 4" p. 115
% Inputs:
% s     = structural shock (T,1)
% A_IRF = some impact matrix you need to input to genIRF
% B     = reduced form beta
% shock_number = this function does the decomposition for one shock at a
% time, and this variable tells it which shock to consider
% variable_number = tells the function the position of the variable for
% which the decomposition is being performed.

T = size(s,1);
H = T;
[IRFs, ~, ~, ~, ~] = genIRFs(A_IRF,0,B,0,H, 0.95, 0.9);

% Here I select from the matrix of IRFs the IRF of the variable we want to
% the shock we want
IRF = squeeze(IRFs(variable_number,:,shock_number));

hd   = zeros(T,1);
hd_i = zeros(T,1);
for t = 1:T
    for i=1:t
    hd_i(i) = IRF(t-i+1)*s(i);
    end
    hd(t) = sum(hd_i);
end