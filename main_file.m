%% This file is the main file of running the first small-scale VAR.

clear all

%% Load data and arrange
data = xlsread('dataset_23_sept_2017','Sheet1','B126:F283');

% % try doing with Forni's dataset: 
% data = xlsread('DataTestRevb.xlsx','Quarterly','b11:bk214');
% code = xlsread('DataTestRevb.xlsx','Quarterly','b5:bk5');
% 
% datalevel = data(2:end,:);
% 
% for o = 1:length(code);
%  if code(o) == 2
%      datalevel(:,o) = log(data(2:end,o));
% 
% elseif code(o) == 3
%      datalevel(:,o) = diff(data(:,o));
%  end
% end
% 
% % In this case, the variables of interest are: 
% % datalevel(:,49) % TFP (log-levels)
% % datalevel(:,55) % Mich index (log-levels)
% data_levels = [datalevel(:,49) datalevel(:,55)];

% Our stuff:
% 1st column: TFP (annual growth rates, quarterly)
% 2nd column: R&D (levels, quarterly)
% 3rd column: SPF long run productivity annual growth expectation (growth
% rates, annual) (DON'T USE RIGHT NOW!)
% 4th column: Michigan index of business conditions 5-years ahead (levels, quarterly)
% 5th column: nonresidential fixed investment in various IT stuff (levels, quarterly)

% Extrapolate annual variables (to make them quarterly)
do_extrap = 0;
if do_extrap == 1
    for i=1:4:length(data)-4
        extrap = linspace(data(i,3),data(i+4,3),5);
        data(i+1:i+3,3) = extrap(2:4);
    end
else % Turn quarterly variables into annual
    % NOT DOING THAT EITHER
end


% Cumulate growth variables to levels (log levels to be precise, b/c growth
% rates are calculated as log diffs)
data_levels(:,1) = cumsum(data(:,1));
data_levels(:,2) = data(:,5); % this series was levels to start out with so we take logs
% 5 means we take IT investment instead of R&D
% data_levels(:,3) = cumsum(data(:,3)); %ignore this since it's SPF
data_levels(:,3) = data(:,4); % the Mich index take the third column in data_levels.


% % Cut out the last line where we have no SPF value
% data_levels = data_levels(1:end-1,:);

%Wrong order in the dataset! Rearrange to have: RD is the last! TFP first and SPF second!
y      = zeros(size(data_levels,1),size(data_levels,2));
y(:,1) = data_levels(:,1);
y(:,2) = data_levels(:,3);
y(:,3) = data_levels(:,2);

data_levels = y;

% figure(1)
% hold on
% plot(data_levels(:,1),'LineWidth',2)
% plot(data_levels(:,2),'LineWidth',2)
% plot(data_levels(:,3),'LineWidth',2)
% grid on
% hold off


%% Run a simple VAR
nlag = 2;
nt = 80;
[beta, c, mu] = quick_var(data_levels,nlag);
news_shock = c*[0 1 0]';
rd_shock   = c*[0 0 1]';
% news_shock = c*[0 1]'; % Forni data check
IR_news = quick_IR(beta, nt, news_shock);
IR_rd   = quick_IR(beta, nt, rd_shock);

%% Bootstrap and generate bootstrapped CIs
nsimul = 5000;
burnin = 200;
T = length(data_levels);
simsample = zeros(size(data_levels,1), size(data_levels,2)); % a single simulated sample
IR_news_boot = zeros(nt,size(data_levels,2),nsimul);
IR_rd_boot = zeros(nt,size(data_levels,2),nsimul);
beta_boot = zeros(nlag*size(data_levels,2),size(data_levels,2),nsimul);

for j=1:nsimul
    % generate a simulated sample
    simsample = quick_boot(beta,mu,T,burnin,0);
    % estimate a bootstrapped beta
    [beta_boot(:,:,j), c_boot, mu_boot] = ...
        quick_var(simsample,nlag);
    news_shock_boot = c_boot*[0 1 0]';
%     news_shock_boot = c_boot*[0 1]'; % Forni check 

    rd_shock_boot   = c_boot*[0 0 1]';
    IR_news_boot(:,:,j) = ...
        quick_IR(beta_boot(:,:,j), nt, news_shock_boot);
    IR_rd_boot(:,:,j)   = ...
        quick_IR(beta_boot(:,:,j), nt, rd_shock_boot);
end


% Do Kilian correction 
beta_i = mean(beta_boot,3); % average bootstrapped beta, see lect. 9 Ryan
beta_boot_k = 2.*beta_boot-beta_i.*ones(size(beta_boot,1), size(beta_boot,2),nsimul);

for j = 1:nsimul
    beta_boot_k_check(:,:,j) = 2*beta_boot(:,:,j) - beta_i;
end
if sum(sum(sum(abs(beta_boot_k - beta_boot_k_check)))) > 10^-6
    error('Kilian-corrected beta_boot is wrong')
end

beta_boot = beta_boot_k;

IR_news_sorted = sort(IR_news_boot,3);
IR_rd_sorted   = sort(IR_rd_boot,3);

ub_news = IR_news_sorted(:,:,ceil(nsimul*0.9));
lb_news = IR_news_sorted(:,:,floor(nsimul*0.1));

ub_rd = IR_rd_sorted(:,:,ceil(nsimul*0.9));
lb_rd = IR_rd_sorted(:,:,floor(nsimul*0.1));

% For checking purposes: an average beta_boot
avg_beta_boot = mean(beta_boot,3);

return
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