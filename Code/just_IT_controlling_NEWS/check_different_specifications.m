clear all

%load twotousandbootstrap_controllingNews_RPzero.mat
load Workspace_Just_IT_controllingNEWS_shortrunimpactsRP_maxNEWS_maxIT_1LAG.mat

data1 = data;
B1 = B;
A1 = A;
res1 = res;
imp1 = impact;
IRFs1 = IRFs;
gam1  = gam_opt;


%load Workspace_Just_IT_controllingNEWS_shortrunimpactsRP_maxNEWS_maxIT_1LAG.mat
load checkcheck3
data2 = data;
B2 = B;
A2 = A;
res2 = res;
imp2 = impact;
IRFs2 = IRFs;
gam2 = gam;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



check_data = sum(sum((data1 - data2).^2));
check_B = sum(sum((B1 - B2).^2))
check_A = sum(sum((A1 - A2).^2))
check_res = sum(sum((res1 - res2).^2))
check_imp = sum(sum((imp1 - imp2).^2))
check_IRFs = sum(sum(sum((IRFs1 - IRFs2).^2)))























