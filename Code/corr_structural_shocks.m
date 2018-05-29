% Investigate correlations between structural shock time series obtained
% from different identification strategies

clear
load 'correct_favourite_macrolunch_structural_shocks.mat'

s_old_news = structural_shock_RyanID(:,2);

load 's_shock_just_IT.mat'

s_justIT = s_shock_just_IT(1:end-1);

load 'structural_shocks_justIT_controllingNEWS.mat'

s_justIT_controlling_news = structural_shock_RyanID(1:end-1,2);

size(s_old_news)
size(s_justIT)
size(s_justIT_controlling_news)

structural_shocks = [s_old_news s_justIT s_justIT_controlling_news];

corr(structural_shocks)