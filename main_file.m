%% This file is the main file of running the first small-scale VAR.

clear all

%% Load data and arrange
data = xlsread('dataset_23_sept_2017','Sheet1','B182:D283');
% 1st column: TFP (annual growth rates, quarterly)
% 2nd column: R&D (levels, quarterly)
% 3rd column: SPF long run productivity annual growth expectation (growth rates, annual)

% Extrapolate annual variables (to make them quarterly)
for i=1:4:length(data)-4
    extrap = linspace(data(i,3),data(i+4,3),5);
    data(i+1:i+3,3) = extrap(2:4); 
end

% Cumulate growth variables to levels (log levels to be precise, b/c growth
% rates are calculated as log diffs)
data_levels(:,1) = cumsum(data(:,1));
data_levels(:,2) = log(data(:,2)); % this series was levels to start out with so we take logs
data_levels(:,3) = cumsum(data(:,3));

% Cut out the last line where we have no SPF value
data_levels = data_levels(1:end-1,:);

%Wrong order in the dataset! RD is the last! TFP first and SPF second!
y      = zeros(size(data_levels,1),size(data_levels,2));
y(:,1) = data_levels(:,1);
y(:,2) = data_levels(:,3);
y(:,3) = data_levels(:,2);

data_levels = y;

figure(1)
hold on
plot(data_levels(:,1),'LineWidth',2)
plot(data_levels(:,2),'LineWidth',2)
% plot(data_levels(:,3),'LineWidth',2)
grid on
hold off

%% Run a simple VAR
nlag = 1;
nt = 80;
[beta, c, mu] = quick_var(data_levels,nlag);
news_shock = c*[0 1 0]';
rd_shock   = c*[0 0 1]';
IR_news = quick_IR(beta, nt, news_shock);
IR_rd   = quick_IR(beta, nt, rd_shock);

%% Bootstrap and generate bootstrapped CIs
nsimul = 5000;
burnin = 200;
T = length(data_levels);
simsample = zeros(size(data_levels,1), size(data_levels,2)); % a single simulated sample
IR_news_boot = zeros(nt,size(data_levels,2),nsimul);
IR_rd_boot = zeros(nt,size(data_levels,2),nsimul);

for j=1:nsimul
    % generate a simulated sample
    simsample = quick_boot(beta,mu,T,burnin,0);
    % estimate a bootstrapped beta
    [beta_boot, c_boot, mu_boot] = ...
        quick_var(simsample,nlag);
    news_shock_boot = c_boot*[0 1 0]';
    rd_shock_boot   = c_boot*[0 0 1]';
    IR_news_boot(:,:,j) = ...
        quick_IR(beta_boot, nt, news_shock_boot);
    IR_rd_boot(:,:,j)   = ...
        quick_IR(beta_boot, nt, rd_shock_boot);
end
IR_news_sorted = sort(IR_news_boot,3);
IR_rd_sorted   = sort(IR_rd_boot,3);

ub_news = IR_news_sorted(:,:,ceil(nsimul*0.99));
lb_news = IR_news_sorted(:,:,floor(nsimul*0.01));

ub_rd = IR_rd_sorted(:,:,ceil(nsimul*0.99));
lb_rd = IR_rd_sorted(:,:,floor(nsimul*0.01));

%% Plot IRFs

periods = 1:nt;

figure
%News shock on TFP
subplot(3,2,1)
hold on
plot(periods, IR_news(:,1),'Color','r') 
plot(periods, ub_news(:,1),'Color','k')
plot(periods, lb_news(:,1),'Color','k')
title('News shock on TFP')
grid on
hold off

%News shock on SPF
subplot(3,2,3)
hold on
plot(periods, IR_news(:,2),'Color','r') 
plot(periods, ub_news(:,2),'Color','k')
plot(periods, lb_news(:,2),'Color','k')
title('News shock on SPF')
grid on
hold off

%News shock on RD
subplot(3,2,5)
hold on
plot(periods, IR_news(:,3),'Color','r') 
plot(periods, ub_news(:,3),'Color','k')
plot(periods, lb_news(:,3),'Color','k')
title('News shock on RD')
grid on
hold off

%RD shock on TFP
subplot(3,2,2)
hold on
plot(periods, IR_rd(:,1),'Color','r') 
plot(periods, ub_rd(:,1),'Color','k')
plot(periods, lb_rd(:,1),'Color','k')
title('RD shock on TFP')
grid on
hold off

%RD shock on SPF
subplot(3,2,4)
hold on
plot(periods, IR_rd(:,2),'Color','r') 
plot(periods, ub_rd(:,2),'Color','k')
plot(periods, lb_rd(:,2),'Color','k')
title('RD shock on SPF')
grid on
hold off

%RD shock on RD
subplot(3,2,6)
hold on
plot(periods, IR_rd(:,3),'Color','r') 
plot(periods, ub_rd(:,3),'Color','k')
plot(periods, lb_rd(:,3),'Color','k')
title('RD shock on RD')
grid on
hold off

disp('Done.')