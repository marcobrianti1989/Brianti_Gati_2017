clear
close all

%Data Reading and Transformation
filename = 'dataset_main';
sheet    = 'Data';
range    = 'B1:B286';
[data, varnames] = read_data2(filename, sheet, range);

TFP = data(:);
T = size(data,1);
time = datetime(1947,4,1) + calquarters(0:T);
ghjk
hold on
plot(time,TFP,'linewidth',1.5)
title('Total Factor Productivity')
grid on
hold off