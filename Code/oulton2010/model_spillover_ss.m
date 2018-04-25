% MODEL_SS - Return the steady state of the model (computed analytically)
%
% usage:
% 
% [ss, param] = model_ss(param)
% STEADY STATE OF OULTON 2010 WITH GROWTH

function [ss,param] = model_spillover_ss(param)

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
h1 = @(wx) (1/chi*wx/ki(wx) + gc/Ki_Kc)/biggamc*ki(wx)^(1-gam)*kc_bar(wx)^(-a)*ki_bar(wx)^(-b);
h1_check = @(wx) (wx/chi + gc/Ki_Kc*ki(wx))/biggamc*ki(wx)^(-gam)*kc_bar(wx)^(-a)*ki_bar(wx)^(-b);
check_h = (h1(150) - h1_check(150))^2;
if check > 10^(-16)
      error('h1 is wrong')
end
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
options = optimoptions(options, 'OptimalityTolerance', 1e-18); 
% Set the Display option to 'iter' and StepTolerance to 1e-4
% options.Display = 'iter';
options.StepTolerance = 1e-18;
objw = @(wx) (wx - chi*c(wx))^2;
wx0 = 3.5;
wstar = fmincon(objw,wx0,[],[],[],[],[],[],[],options);
wstar
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


%Put the ss values in a vector consistent with Y and X vectors in model.m
xx  = [kc ki biggamc biggami ...
    c biggamc biggami yc yi h p kc2 ki2 ki kc]; 
yy  = [yc yi c ic it w rc ri h h1 h2 kc1 kc2 ki1 ki2 p expgc expgi ...
    expgc expgi expgc expgi 1 p expgc expgi];

ss  = [yy xx];



