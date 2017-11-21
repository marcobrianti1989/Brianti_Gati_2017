clear all
close all

filename = 'dataset_main';
sheet    = 'Data';
range    = 'B1:M286';
[data2, variables] = read_da ta2(filename, sheet, range);