clear
close all

% Create the correct path
base_path = pwd;
if exist([base_path '\Data'], 'dir')
      addpath([base_path '\Data']) %for Microsoft
else
      addpath([base_path '/Data']) %for Mac
end

% Reading Data
filename                    = 'Compustat_Data_Nov_2018';
sheet                       = 'WRDS';
range                       = 'A1:R629747';
[dataset, var_names]        = xlsread(filename, sheet, range);
tf                          = isreal(dataset);
if tf == 0
      warning('Dataset has complex variables in it.')
end

% Assess names to each variable as an array
CompanyKey    = dataset(:,1);
Quarters      = var_names(2:end,11);
Assets        = dataset(:,13);
CapSoftware   = dataset(:,14);
Equity        = dataset(:,15);
Investment    = dataset(:,16);
RDexpense     = dataset(:,17);

time = [1961:0.25:2018.25]';
year = [1961:1:2018]';
qrt  = [1:1:4]';
j    = 1;
for iy = 1:length(year)
      for iq = 1:length(qrt)
            index = find(strcmp(Quarters, [num2str(year(iy)) 'Q' num2str(qrt(iq))]));
            for ii = 1:length(index)
                  subInv(ii)                               = Investment(index(ii));
                  subAssets(ii)                            = Assets(index(ii));
                  subSoft(ii)                              = CapSoftware(index(ii));
                  subEqui(ii)                              = Equity(index(ii));
                  subRD(ii)                                = RDexpense(index(ii));
                  if isnan(subInv(ii)) == 1 || isnan(subAssets(ii)) == 1 || isnan(subSoft(ii)) == 1 || isnan(subEqui(ii)) == 1 || isnan(subRD(ii)) == 1
                        subInv(ii)                 = [];
                        subAssets(ii)              = [];
                        subSoft(ii)                = [];
                        subEqui(ii)                = [];
                        subRD(ii)                  = [];
                  end
            end
            % Cash
            sumInv                               = sum(subInv);
            AggInv(j)                            = sumInv;
            % Total Assets
            sumAssets                            = sum(subAssets);
            AggAssets(j)                         = sumAssets;
            % Total Equity
            sumSoft                              = sum(subSoft);
            AggSoft(j)                           = sumSoft;
            % Cash + Short-Term Investment
            sumEqui                              = sum(subEqui);
            AggEqui(j)                           = sumEqui;
            % RD Investment
            sumRD                                = sum(subRD);
            AggRD(j)                             = sumRD;
            % Counter
            j                                    = j + 1;
      end
end

% Remove last two quarters since there are no data on that.
AggInv              = AggInv(2:end-2);
AggEqui             = AggEqui(2:end-2);
AggSoft             = AggSoft(2:end-2);
AggAssets           = AggAssets(2:end-2);
AggRD               = AggRD(2:end-2);
time                = time(2:end); %this is just to show that now time is aligned with other variables

% Remove trend using the seven term henderson filter
[AggEquityDt, ~]    = seven_term_henderson_filter(AggEqui);
[AggInvestDt, ~]    = seven_term_henderson_filter(AggInv);
[AggSoftDt, ~]      = seven_term_henderson_filter(AggSoft);
[AggAssetsDt, ~]    = seven_term_henderson_filter(AggAssets);
[AggRDDt, ~]        = seven_term_henderson_filter(AggRD);

RD2Assets           = (AggRDDt./AggAssetsDt)';
RD2Inv              = (AggRDDt./AggInvestDt)';
RD2Equity           = (AggRDDt./AggEquityDt)';
Soft2Assets         = (AggSoftDt./AggAssetsDt)';
Soft2Inv            = (AggSoftDt./AggInvestDt)';
Soft2Equity         = (AggSoftDt./AggEquityDt)';

corr([RD2Assets(171:end) RD2Equity(171:end) RD2Inv(171:end)]);
corr([Soft2Assets(171:end) Soft2Equity(171:end) Soft2Inv(171:end)]);

figure(1)
plot(time,RD2Assets,'linewidth',1.5)
hold on
plot(time,Soft2Assets,'linewidth',1.5)
grid on
legend('RD2Assets','Soft2Assets')

figure(2)
plot(time,RD2Inv,'linewidth',1.5)
hold on
plot(time,Soft2Inv,'linewidth',1.5)
grid on
legend('RD2Inv','SoftInv')

figure(3)
plot(time,RD2Equity,'linewidth',1.5)
hold on
plot(time,Soft2Equity,'linewidth',1.5)
grid on
legend('RD2Inv','SoftInv')



