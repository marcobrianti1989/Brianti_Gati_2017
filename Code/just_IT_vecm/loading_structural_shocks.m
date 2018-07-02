clear
close all

load main_just_IT_noHOURS_1LAG.mat
sICT1 = s_shock_just_IT;

load Workspace_Just_IT_SVAR_1LAG.mat
sICT2 = s_shock_just_IT;

%just for Microsoft (for now)
addpath('C:\Users\th3\Documents\BG_2017\Code\just_IT_controlling_NEWS') 

load Workspace_Just_IT_controllingNEWS_shortrunimpactsRP_maxNEWS_maxIT_1LAG.mat
sICT3      = structural_shock_RyanID(:,2);
sNEWS1     = structural_shock_RyanID(:,1);

load Workspace_Just_IT_controllingNEWS_zeroimpactIT_maxNEWS_maxIT_1LAG.mat
sICT4      = s_shock_just_IT(:,2);
sNEWS2     = s_shock_just_IT(:,1);

load Workspace_Just_IT_SVAR_2Variables_1LAG.mat
sICT5      = s_shock_just_IT;

load Workspace_Just_IT_SVAR_3Variables_1LAG.mat
sICT6      = s_shock_just_IT;

quarter = 1989.25:0.25:2017;
quarter = quarter';
corr(sICT1,sICT6)
plot(quarter,sICT1,'linewidth',1.8)
hold on
plot(quarter,sICT6,'--','linewidth',2.3)
set(gcf,'color','w');
grid on
LEG = legend('5-dimensional system','3-dimesional system');
LEG.FontSize = 20;
legend boxoff
title('Series of ICT shocks for different specifications','fontsize',40)
axis tight


