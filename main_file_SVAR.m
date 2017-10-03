%Run a SVAR and identify news shocks and R&D shocks

% Marco Brianti, Laura Gáti, Sep 28 2017

clear all
close all

data = xlsread('dataset_23_sept_2017','Sheet1','B126:H283');
% Cumulate growth variables to levels (log levels to be precise, b/c growth
% rates are calculated as log diffs)
data_levels(:,1) = cumsum(data(:,1)); % TFP
% data_levels(:,1) = data(:,1); % TFP % I tried TFP grt.

data_levels(:,3) = log(data(:,5)); % this series was levels to start out with so don't cumsum <-- taking logs here induces stationarity of VAR - DISCUSS! If VAR nonstat and not cointegrated, estimation not possible.
% 5 means we take IT investment instead of R&D
data_levels(:,2) = data(:,4); % the Mich index 
data_levels(:,4) = data(:,6); % real GDP % whether this guy's in logs or not doesn't seem to make a diff
data_levels(:,5) = data(:,7); % real cons % whether this guy's in logs or not doesn't seem to make a diff

% Have an initial look at data
figure
hold on
plot(data_levels(:,1),'k')
plot(data_levels(:,2), 'b')
plot(data_levels(:,3), 'r')
plot(data_levels(:,4), 'g')
grid on 
hold off

%Inputs
nlags = 4;
nburn = 200;
nsimul = 5000;
nvar = size(data_levels,2);

% Run VAR imposing Cholesky
[A,B,res] = sr_var(data_levels, nlags);

%Checking if the VAR is stationary
test_stationarity(B');

asdfgh

% Generate bootstrapped data samples
dataset_boot = data_boot(B, nburn, res, nsimul); % <--- TO DO: draw shocks in blocks

% Redo VAR nsimul times on the bootstrapped datasets
A_boot = zeros(nvar,nvar,nsimul);
B_boot = zeros(nvar*nlags+1,nvar,nsimul);
for i_simul = 1:nsimul
    [A_boot(:,:,i_simul), B_boot(:,:,i_simul), ~] = sr_var(dataset_boot(:,:,i_simul), nlags);
end

%Kilian correction - IT IS NOT WORKING VERY NICELY. DONT KNOW WHY!
% for i_simul = 1:nsimul
%       average_B_boot = mean(B_boot,3);
%       average_A_boot = mean(A_boot,3);
%       B_boot(:,:,i_simul) = 2*B_boot(:,:,i_simul) - average_B_boot;
%       A_boot(:,:,i_simul) = 2*A_boot(:,:,i_simul) - average_A_boot;
% end

% Calculate IRFs, bootstrapped CI and plot them
h=80;
which_shock = [2 3];
names = {'News shock','R&D shock'}; % shock names in order of appearance
% varnames = {'TFP','Mich index','R&D'} ; % variable names in order of appearance
% varnames = {'TFP','Mich index','R&D', 'GDP'} ; % alternative specs 1
varnames = {'TFP','Mich index','R&D', 'GDP', 'C'} ; % alternative specs 2

sig = 0.90; % significance level
plotIRFs(A,A_boot,B,B_boot,h,which_shock, names, varnames,sig);

%------------------------------------------------------------------------------------------


