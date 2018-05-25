 % MODEL_SS - Return the steady state of the model (computed analytically)
%
% usage:
% 
% [ss, param] = model_ss(param)
% STEADY STATE OF OULTON 2010 WITH GROWTH

function [ss,param] = model_spillover_ss_news(param)

%Parameters from param object
bet      = param.bet;  %Never call anything beta...it's a matlab function
a        = param.a; % share of standard capital
b        = param.b; % share of IT capital
biggami  = param.biggami; % growth rate of equipment productivity
biggamc  = param.biggamc;    % growth rate of neutral technology
di       = param.di; % depreciation rate of capital IT
dc       = param.dc; % depreciation rate of capital standard
chi      = param.chi; % preference parameter.
gam      = param.gam; % spillover elasticity
siggami  = param.siggami; % variance of BIGGAMI
sige     = param.sige; % variance of signal on BIGGAMI
%Use closed form expressions for the ss values. 

%Step 0 - always run it!
p        = (biggamc)/(biggami); 
expgc    = biggamc^((1-b-gam)/(1-a-b-gam))*biggami^((b+gam)/(1-a-b-gam));
expgi    = biggamc^(a/(1-a-b-gam))*biggami^((1-a)/(1-a-b-gam));
gc       = expgc - (1-dc); 
gi       = expgi - (1-di);
rc       = 1/bet * expgc - (1-dc); 
ri       = (1/bet* expgi -(1-di))*p; 

%Step 1
kc_bar          = @(wx) wx/rc*a/(1-a-b); %kc = kc1 = kc2 = Kc/h = Kc1/h1 = Kc2/h2
ki_bar          = @(wx) wx/ri*b/(1-a-b); %ki = ki1 = ki2 = Ki/h = Ki1/h1 = Ki2/h2
ki              = @(wx) (rc/(a*biggamc)*kc_bar(wx)^(1-a)*ki_bar(wx)^(-b))^(1/gam);
Ki_check        = @(wx) (ri/(b*biggamc)*kc_bar(wx)^(-a)*ki_bar(wx)^(1-b))^(1/gam); %check 
Ki_check2       = @(wx) (wx/((1-a-b)*biggamc)*kc_bar(wx)^(-a)*ki_bar(wx)^(-b))^(1/gam); %check
check           = ((ki(15) - Ki_check(15)) + (ki(15) - Ki_check2(15)))^2;
if check > 10^(-16)
      error('Ki is wrong')
end
Ki_Kc = b/a*rc/ri; %Ki_Kc = Ki/Kc = Ki1/Kc1 = Ki2/Kc2 = ki/kc

%Step 2
h2 = @(wx) gi/biggami*ki(wx)^(1-gam)*kc_bar(wx)^(-a)*ki_bar(wx)^(-b);
ki2 = @(wx) ki_bar(wx)*h2(wx);
kc2 = @(wx) kc_bar(wx)*h2(wx);

%Step 3
% Case of V(h) = -chi*H
h1 = @(wx) (1/chi*wx/ki(wx) + gc/Ki_Kc)/biggamc*ki(wx)^(1-gam)*kc_bar(wx)^(-a)*ki_bar(wx)^(-b);
h1_check = @(wx) (wx/chi + gc/Ki_Kc*ki(wx))/biggamc*ki(wx)^(-gam)*kc_bar(wx)^(-a)*ki_bar(wx)^(-b);
% % Case of V(h) = chi*log(1-h)
% h1 = @(wx) (wx*h2(wx)/chi - wx/chi + gc/Ki_Kc*ki(wx)) / (biggamc*ki(wx)^(gam)*kc_bar(wx)^(a)*ki_bar(wx)^(b) -wx/chi);
% h1_check = @(wx) (biggamc*ki(wx)^(gam)*kc_bar(wx)^(a)*ki_bar(wx)^(b) -wx/chi)^(-1) ...
%     *(gc/Ki_Kc*ki(wx) + biggamc*ki(wx)^(gam)*kc_bar(wx)^(a)*ki_bar(wx)^(b)*h2(wx) -wx/chi) - h2(wx);
% Case of V(h) = chi/2*H^2
h1 = @(wx) (gc/Ki_Kc*ki(wx) + biggamc*ki(wx)^(gam)*kc_bar(wx)^(a)*ki_bar(wx)^(b)*h2(wx)  ...
    + sqrt(  (gc/Ki_Kc*ki(wx) + biggamc*ki(wx)^(gam)*kc_bar(wx)^(a)*ki_bar(wx)^(b)*h2(wx))^2  ...
    + 4*biggamc*ki(wx)^(gam)*kc_bar(wx)^(a)*ki_bar(wx)^(b) *wx/chi)) ...
    / (2*biggamc*ki(wx)^(gam)*kc_bar(wx)^(a)*ki_bar(wx)^(b)) -h2(wx);

% check_h = (h1(150) - h1_check(150))^2;
% if check_h > 10^(-16)
%       error('h1 is wrong')
% end
ki1 = @(wx) ki_bar(wx)*h1(wx);
kc1 = @(wx) kc_bar(wx)*h1(wx);

%Step 4
h = @(wx) h1(wx) + h2(wx);
kc = @(wx) kc1(wx) + kc2(wx);
yc = @(wx) biggamc*ki(wx)^(gam)*h1(wx)^(1-a-b)*kc1(wx)^(a)*ki1(wx)^(b);
yi = @(wx) biggami*ki(wx)^(gam)*h2(wx)^(1-a-b)*kc2(wx)^(a)*ki2(wx)^(b);
ii = @(wx) gi*ki(wx);
ic = @(wx) gc*kc(wx);
c = @(wx) yc(wx) - ic(wx);

%Step 5 - Minimization
%Objective: w = @(wx) chi*c(wx)
options = optimoptions('fmincon'); 
% Set OptimalityTolerance to 1e-15
%options = optimoptions(options, 'OptimalityTolerance', 1e-18); 
% Set the Display option to 'iter' and StepTolerance to 1e-4
options.Display = 'iter';
%options.StepTolerance = 1e-18;
% objw = @(wx) (wx/chi - c(wx))^2; % case of linear V(H)
objw = @(wx) (chi*h(wx) - wx/c(wx))^2; % case of  V(H) quadratic
wx0 = 100;
wstar = fmincon(objw,wx0,[],[],[],[],[],[],[],options);
wstar
% A check on how h(w) behaves
wgrid = linspace(0.5,0.9,100);
for j=1:length(wgrid)
w_figure(j) = objw(wgrid(j));
h_figure(j) = h(wgrid(j));
c_figure(j) = c(wgrid(j));
ic_figure(j) = ic(wgrid(j));
yc_figure(j) = yc(wgrid(j));
ki_figure(j) = ki(wgrid(j));
kc_figure(j) = kc(wgrid(j));
ki1_figure(j) = ki1(wgrid(j));
kc1_figure(j) = kc1(wgrid(j));
ki2_figure(j) = ki2(wgrid(j));
kc2_figure(j) = kc2(wgrid(j));
h1_figure(j) = h1(wgrid(j));
end

% plot(wgrid, h_figure); hold on
 %plot(wgrid, w_figure)
 %plot(wgrid, c_figure)
 %legend('objective')
% Try to use also fzero
% objw2 = @(wx)  wx - chi*c(wx);
% wstar2 = fzero(objw2,wx0);
res_obj_wstar = objw(wstar)

%Step 5
kc = kc(wstar); %small kc back to be original Kc - drop previous notation
ki = ki(wstar); %small ki back to be original Ki - drop previous notation
yc = yc(wstar);
yi = yi(wstar);
c = c(wstar);
ic = ic(wstar);
it = ii(wstar);
w = wstar;
h = h(wstar);
h1 = h1(wstar);
h2 = h2(wstar);
kc1 = kc1(wstar); %small kc1 back to be original Kc1 - drop previous notation
kc2 = kc2(wstar); %small kc2 back to be original Kc2 - drop previous notation
ki1 = ki1(wstar); %small ki1 back to be original Ki1 - drop previous notation
ki2 = ki2(wstar); %small ki2 back to be original Ki2 - drop previous notation

% Computation of NIPA-consistent GDP and TFP
gamgdp = expgc + p*expgi;
% gdp = yc + p*yi;
tfp = gamgdp - 2*(1-a-b)*w -2*a*rc - 2*b*ri;

%Put the ss values in a vector consistent with Y and X vectors in model.m
xx  = [kc ki biggamc biggami ...
    c biggamc biggami yc yi h p kc2 ki2 ki kc ri w 0 0 0 0 0 0 0 0 0 1 ...
    biggami biggami ...
    1]; 
yy  = [yc yi c ic it w rc ri h h1 h2 kc1 kc2 ki1 ki2 p expgc expgi ...
    expgc expgi expgc expgi 1 expgc/expgi expgc expgi expgc/expgi expgc gamgdp tfp];

ss  = [yy xx];



