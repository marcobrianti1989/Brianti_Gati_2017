clear all
close all


%Data Reading and Transformation
filename = 'dataset_IT_prices';
sheet    = 'Sheet1';
range    = 'B1:F324';

data = xlsread(filename,sheet,range);
data(:,1) = data(:,1)./10;

[T,n] = size(data);

time = datetime(1990,1,1) + calmonths(0:T-1);

figure(1)
for i = 1:n
    if i == 2
        % don't plot this one, it's too short
    else
plot(time, data(:,i), 'linewidth', 2); hold on
grid on
    end
end
title('4 different IT price indeces (PPI subsectors) - blue is normalized by 1/10')

h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperPosition', [0 0 11 8]); %[1 1 28 19]
filename = 'IT_prices trends 14 Nov 2017';
print(filename,'-dpdf')


