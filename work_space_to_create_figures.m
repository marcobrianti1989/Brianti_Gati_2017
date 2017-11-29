clear
close all

%Data Reading and Transformation
filename = 'dataset_main';
sheet    = 'Data';
range    = 'B1:P286';
[data, varnames] = read_data2(filename, sheet, range);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING TFP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_TFP = datetime(1989,4,1) + calquarters(0:size(data,1)-1);
gcf = figure(1);
set(gcf,'color','w');
set(gcf,'position',[1 41 1920 963])
hold on
plot(time_TFP,data(:,1),'linewidth',1.5)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 24)
brk = time_TFP(68);
y1 = get(gca,'ylim');
plot([brk brk],y1)
txt1 = '2006:q1';
xtx = brk + 70;
ytx = 345;
text(xtx,ytx,txt1,'fontsize',24)
leg = legend('Utilization-adjusted TFP');
set(leg,'fontsize',28)
grid on
hold off
%invoke_export_fig('TFP','macrolunch')

clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% PLOTTING R&D GROWTH %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = 'dataset_main';
sheet    = 'Data';
range    = 'B1:C286';
[data, varnames] = read_data2(filename, sheet, range);
data2 = [NaN NaN ; diff(log(data))];
time_RD = datetime(1947,4,1) + calquarters(0:size(data,1)-1);
time_RD = time_RD';


figure(1)
gcf = figure(1);
set(gcf,'color','w');
set(gcf,'position',[1 41 1920 963])
hold on
plot(time_RD(80:end),data2(80:end,:),'linewidth',1.5)
leg = legend('Utilization-adjusted TFP','Real R&D expenditure');
set(leg,'fontsize',28)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 24)
grid on
hold off
%invoke_export_fig('RD','macrolunch')
%corr(data(80:end,1),data(80:end,2))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% PLOTTING R&D LEVEL %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

filename = 'dataset_main';
sheet    = 'Data';
range    = 'B1:C286';
[data, varnames] = read_data2(filename, sheet, range);
data2 = log(data);
time_RD = datetime(1947,4,1) + calquarters(0:size(data,1)-1);
time_RD = time_RD';


figure(1)
gcf = figure(1);
set(gcf,'color','w');
set(gcf,'position',[1 41 1920 963])
hold on
yyaxis right
plot(time_RD(192:end),data2(192:end,2),'linewidth',1.5)
yyaxis left
plot(time_RD(192:end),data2(192:end,1),'linewidth',1.5);
leg = legend('Utilization-adjusted TFP','Real R&D expenditure');
brk1 = time_RD(236);
brk2 = time_RD(245);
y1 = get(gca,'ylim');
plot([brk1 brk1],y1,'--','linewidth',1.5,'color','k')
plot([brk2 brk2],y1,'--','linewidth',1.5,'color','k')
set(leg,'fontsize',28)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 24)
grid on
hold off
invoke_export_fig('RD_level','macrolunch')
%corr(data(80:end,1),data(80:end,2))














