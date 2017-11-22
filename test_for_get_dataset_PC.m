% Load in all the series on FRED (81 series) and save them to individual
% excel files

% It takes about 50 sec.

clear
close all

tic
get_dataset_test_pc;
toc