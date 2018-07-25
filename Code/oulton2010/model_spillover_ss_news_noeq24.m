 % MODEL_SS - Return the steady state of the model (computed analytically)
%
% usage:
% 
% [ss, param] = model_ss(param)
% STEADY STATE OF OULTON 2010 WITH GROWTH
% SAME AS model_spillover_ss_news.M, ONLY AVOIDS EVALUATING EQ. 24 FROM
% STAIONARIZED_SYSTEM_SPILLOVER.PDF, WHICH INCORPORATES 1/GAM, GIVING RISE
% TO THE UNDEFINED OBJECTIVE FUNCTION WHEN GAM = 0.
% 25 July 2018.

function [ss,param] = model_spillover_ss_news_noeq24(param)

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

% NOTE: % x = [wx kix]

%Step 1
kc_bar          = @(x) x(1)/rc*a/(1-a-b); %kc = kc1 = kc2 = Kc/h = Kc1/h1 = Kc2/h2
ki_bar          = @(x) x(1)/ri*b/(1-a-b); %ki = ki1 = ki2 = Ki/h = Ki1/h1 = Ki2/h2 
% ki              = @(wx) (rc/(a*biggamc)*kc_bar(wx)^(1-a)*ki_bar(wx)^(-b))^(1/gam); % THIS IS WHERE GAM=0 MAKES OBJECTIVE UNDEFINED.
% Ki_check        = @(wx)(ri/(b*biggamc)*kc_bar(wx)^(-a)*ki_bar(wx)^(1-b))^(1/gam); 
% Ki_check2       = @(wx) ( wx/((1-a-b)*biggamc)*kc_bar(wx)^(-a)*ki_bar(wx)^(-b) )^(1/gam); %check
% check           = ((ki(15) - Ki_check(15)) + (ki(15) - Ki_check2(15)))^2;
% if check > 10^(-16)
%       warning('Ki is wrong')
% end
Ki_Kc = b/a*rc/ri; %Ki_Kc = Ki/Kc = Ki1/Kc1 = Ki2/Kc2 = ki/kc

%Step 2
h2 = @(x) gi/biggami*x(2)^(1-gam)*kc_bar(x)^(-a)*ki_bar(x)^(-b);
ki2 = @(x) ki_bar(x)*h2(x);
kc2 = @(x) kc_bar(x)*h2(x);

%Step 3
% Case of V(h) = -chi*H
% h1 = @(wx) (1/chi*wx/ki(wx) + gc/Ki_Kc)/biggamc*ki(wx)^(1-gam)*kc_bar(wx)^(-a)*ki_bar(wx)^(-b);
% h1_check = @(wx) (wx/chi + gc/Ki_Kc*ki(wx))/biggamc*ki(wx)^(-gam)*kc_bar(wx)^(-a)*ki_bar(wx)^(-b);
% % Case of V(h) = chi*log(1-h)
% h1 = @(wx) (wx*h2(wx)/chi - wx/chi + gc/Ki_Kc*ki(wx)) / (biggamc*ki(wx)^(gam)*kc_bar(wx)^(a)*ki_bar(wx)^(b) -wx/chi);
% h1_check = @(wx) (biggamc*ki(wx)^(gam)*kc_bar(wx)^(a)*ki_bar(wx)^(b) -wx/chi)^(-1) ...
%     *(gc/Ki_Kc*ki(wx) + biggamc*ki(wx)^(gam)*kc_bar(wx)^(a)*ki_bar(wx)^(b)*h2(wx) -wx/chi) - h2(wx);
% Case of V(h) = chi/2*H^2
h1 = @(x) (gc/Ki_Kc*x(2) + biggamc*x(2)^(gam)*kc_bar(x)^(a)*ki_bar(x)^(b)*h2(x)  ...
    + sqrt(  (gc/Ki_Kc*x(2) + biggamc*x(2)^(gam)*kc_bar(x)^(a)*ki_bar(x)^(b)*h2(x))^2  ...
    + 4*biggamc*x(2)^(gam)*kc_bar(x)^(a)*ki_bar(x)^(b) *x(1)/chi)) ...
    / (2*biggamc*x(2))^(gam)*kc_bar(x)^(a)*ki_bar(x)^(b) -h2(x);
aa = @(x) biggamc * x(2)^gam * kc_bar(x)^a * ki_bar(x)^b; % my notes4, p. 54
bb = @(x) biggamc * x(2)^gam * kc_bar(x)^a * ki_bar(x)^b * h2(x) - Ki_Kc^(-1)*x(2)*gc;
cc = @(x) -(x(1)/chi + Ki_Kc^(-1)*x(2)*gc * h2(x));
h1_check = @(x) (-bb(x) + sqrt(bb(x)^2 -4*aa(x)*cc(x)))/(2*aa(x)); % taking the positive root
check_h = (h1([150,150]) - h1_check([150,150]))^2;
if check_h > 10^(-14)
      warning('h1 is wrong')
end
ki1 = @(x) ki_bar(x)*h1(x);
kc1 = @(x) kc_bar(x)*h1(x);

%Step 4
h  = @(x) h1(x) + h2(x);
kc = @(x) kc1(x) + kc2(x);
yc = @(x) biggamc*x(2)^(gam)*h1(x)^(1-a-b)*kc1(x)^(a)*ki1(x)^(b);
yi = @(x) biggami*x(2)^(gam)*h2(x)^(1-a-b)*kc2(x)^(a)*ki2(x)^(b);
ii = @(x) gi*x(2);
ic = @(x) gc*kc(x);
c  = @(x) yc(x) - ic(x);

%Step 5 - Minimization
%Objective: w = @(wx) chi*c(wx)
options = optimoptions('fmincon'); 
% Set OptimalityTolerance to 1e-15
% options = optimoptions(options, 'OptimalityTolerance', 1e-18); 
% Set the Display option to 'iter' and StepTolerance to 1e-4
options.Display = 'iter';
%options.StepTolerance = 1e-18;
% objw = @(wx) (wx/chi - c(wx))^2; % case of linear V(H)
% x = [wx kix]
objw = @(x) (chi*h(x) - x(1)/c(x))^2; % case of  V(H) quadratic
wx0 = 0.1;
kix0 = 0.1;
x0 = [0.3, 0.3]; % [0.3, 0.3]
opt = fmincon(objw,x0,[],[],[],[],[],[],[],options);
wstar = opt(1)
kistar= opt(2)

res_obj = objw(opt)

%Step 5
kc = kc(opt); %small kc back to be original Kc - drop previous notation
ki = kistar; %small ki back to be original Ki - drop previous notation
yc = yc(opt);
yi = yi(opt);
c = c(opt);
ic = ic(opt);
it = ii(opt);
w = wstar;
h = h(opt);
h1 = h1(opt);
h2 = h2(opt);
kc1 = kc1(opt); %small kc1 back to be original Kc1 - drop previous notation
kc2 = kc2(opt); %small kc2 back to be original Kc2 - drop previous notation
ki1 = ki1(opt); %small ki1 back to be original Ki1 - drop previous notation
ki2 = ki2(opt); %small ki2 back to be original Ki2 - drop previous notation

% Computation of NIPA-consistent GDP and TFP
wi = p*yc/(yc+p*yi);
gamgdp = expgc^(1-wi) * expgi^wi;
gamtfp = expgi^gam * biggamc^(1-wi) * biggami^wi;

%Put the ss values in a vector consistent with Y and X vectors in model.m
xx  = [kc ki biggamc biggami ...
    c biggamc biggami yc yi h p kc2 ki2 ki kc ri w 0 0 0 0 0 0 0 0 0 ...
    1 ...
    biggami biggami ...
    1 ...
    h1 h2 kc1 ki1 it]; 
yy  = [yc yi c ic it w rc ri h h1 h2 kc1 kc2 ki1 ki2 p expgc expgi ...
    expgc expgi expgc expgi 1 expgc/expgi expgc expgi expgc/expgi expgc ...
    gamgdp gamtfp expgc ...
    1 1 expgc expgi expgi]; % need to check if grt rate of h h1 h2 is really 1 in st.st... I think it should be 0. 

ss  = [yy xx];

save ss_old.mat xx yy

