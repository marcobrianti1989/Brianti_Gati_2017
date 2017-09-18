clear all
close all

dataset000 = xlsread('Data_CEE','quarterly','b5:j283'); % reading data
dataset001 = log(dataset000(:,1:6)); % taking log
dataset002 = log(dataset000(:,8:9)); % taking log
m2 = diff(dataset002(:,2)); % construct m2 growth
m2 = [zeros(1,1);  m2];
dataset02 = [dataset001 dataset000(:,7) dataset002(:,1) m2]; % whole dataset
dataset = dataset02(71:end,1:end); %260
lag_number = 4;
total_extractions = 100;
which_shock = 7;

tic
[A, B, B_boot, B_boot_Kilian, mean_B_boot_check, mshock, A99, A1, A95, A5, A86, A16] = ...
    cholboot(dataset,lag_number,which_shock,total_extractions);
toc