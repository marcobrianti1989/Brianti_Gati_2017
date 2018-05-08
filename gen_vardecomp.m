function [vardec] = gen_vardecomp(IRFs,m,H)
% Perform variance decomposition after estimating a VAR
% m = the horizon of the variance decomposition (can equal, but not exceed H)
% H = the number of IRs generated in plotIRFs.m

% Marco Brianti, Laura Gati    Oct 3, 2017

if m > H
    error('number of periods for variance decomposition cannot exceed the number of IRs generated')
end

nvar = size(IRFs,1);
nshocks = nvar;
% cut out those IRs we don't need
IRFs = IRFs(:,1:m,:);

identification = 'complete';
switch identification
    case 'complete'
        % Generate denominator
        den = sum(sum(IRFs.^2,3),2); % one value for each variable, so (nvar x 1)
        
    case 'partial'
        % Alternative denominator for partial identification case
        
end

% Generate numerator
num = sum(IRFs.^2,2);
num = squeeze(num); % get rid of unnecessary dimension
% one value for each variable AND each shock, so (nvar x nshocks)

% Generate variance decomposition table
vardec = zeros(nvar, nshocks);
for i_var = 1:nvar
    vardec(i_var,:) = num(i_var,:)/den(i_var);
end

% check if they sum up to 1
if sum(abs(sum(vardec,2) - ones(nvar,1))) > 10^(-14)
    error('Variance decomposition doesn''t sum to 1')
end
