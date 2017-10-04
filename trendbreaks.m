%% Identifying breaks in trend

% This file aims to construct our "Empirical Fact" figure: a figure that
% shows that some measure of "investment in learning" (--> that increases
% the endogenous component of TFP), which we simply refer to as "IT",
% exhibits a break in the trend around 2005-6, so that what is relevant for
% the evolution of TFP is NOT the cyclical variation in this variable, but
% the trend.

%% Read in a series
clear all

%Technical name
file = 'data_trendbreaks.xlsx';
sheet ='Sheet1';
range = 'B2:E333';
time_range = 'A2:A333';

%Downloading Data
dataset = xlsread(file,sheet,range);
nvar = size(dataset,2);
varnames = {'Lessors of Nonfinancial Intangible Assets','Value of Manufacturers New Orders for IT Industries','Value of Manufacturers Total Inventories for IT Industries', 'San Francisco Tech Pulse'};

for i_var = 1:nvar
data = dataset(:,1);
varname = varnames{i_var};

%Defining the variable time to label the x-axis
T = size(data,1);
time = datetime(1990,1,1) + calmonths(0:T-1); %ALWAYS CHECK IF THE INITIAL PERIOD IS THE SAME

% %Alternative way to define time to label x-axis
% date_numbers_excel = xlsread(file,sheet,time_range);
% datenumbers = x2mdate(date_numbers_excel,0);
% periods = datestr(datenumbers, 'yyyy-mm');

% Create a moving average smoothed trend
weigths = [1/24;repmat(1/12,11,1);1/24];% choosing 1/24 weights for end terms, 1/12 weight for interior ones b/c data is monthly
weights2 = [1/60; repmat(1/30,28,1); 1/60 ];
trend = conv(data,weigths,'valid'); % 'valid'  leaves out the end observations as those cannot be smoothed and would thus distort things
trend = vertcat([nan nan nan nan nan nan]', trend, [nan nan nan nan nan nan]');

%Plot moving average
figure(i_var*2 - 1)
hold on
plot(time, data,'b', 'linewidth', 2)
plot(time, trend, 'r', 'linewidth', 1.5)
grid on
legend('Original data', 'Trend')
x1 = time(133);
y1=get(gca,'ylim');
plot([x1 x1],y1)
title([varname])
hold off


% Fit a quadratic trend from OLS reg
t = (1:T)';
X = [ones(T,1) t t.^2];
b = X\data;
tH = X*b; % quadratic trend estimate

%Plot Fit quadratic trend
figure(i_var*2)
hold on
plot(time, data,'b', 'linewidth', 2)
% plot(datenumbers, tH, 'r', 'linewidth', 2)
plot(time, tH, 'r', 'linewidth', 1.5)
grid on
legend('Original data', 'Trend')
x1 = time(133);
y1=get(gca,'ylim');
plot([x1 x1],y1)
title([varname])
hold off

end

%% TO DO: tests for structural breaks
