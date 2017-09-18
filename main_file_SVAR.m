clear all 
close all

dataset000 = xlsread('Data_CEE','quarterly','b5:j283'); % reading data
dataset001 = log(dataset000(:,1:6)); % taking log
dataset002 = log(dataset000(:,8:9)); % taking log
m2 = diff(dataset002(:,2)); % construct m2 growth
m2 = [zeros(1,1);  m2];
dataset02 = [dataset001 dataset000(:,7) dataset002(:,1) m2]; % whole dataset
dataset = dataset02(260:end,1:end);
lag_number = 4;
total_extractions = 100;
which_shock = 7;

[A, mshock, A99, A1, A95, A5, A86, A16] = ...
      cholboot(dataset,lag_number,which_shock,total_extractions);