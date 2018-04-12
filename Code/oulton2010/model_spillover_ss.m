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
kc          = @(wx) wx/rc*a/(1-a-b); %kc = kc1 = kc2 = Kc/h = Kc1/h1 = Kc2/h2
ki          = @(wx) wx/ri*b/(1-a-b); %ki = ki1 = ki2 = Ki/h = Ki1/h1 = Ki2/h2
Ki          = @(wx) (rc/(a*biggamc)*kc(wx)^(1-a)*ki(wx)^(-b))^(1/gam);
Ki_check    = @(wx) (ri/(b*biggamc)*kc(wx)^(-a)*ki(wx)^(1-b))^(1/gam); %check 
Ki_check2   = @(wx) (wx/((1-a-b)*biggamc)*kc(wx)^(-a)*ki(wx)^(-b))^(1/gam); %check
check       = ((Ki(15) - Ki_check(15)) + (Ki(15) - Ki_check2(15)))^2;
if check > 10^(-20)
      error('Ki is wrong')
end
Ki_Kc = b/a*rc/ri; %Ki_Kc = Ki/Kc = Ki1/Kc1 = Ki2/Kc2 = ki/kc

%Step 2
h2 = @(wx) gi/biggami*Ki(wx)^(1-gam)*kc(wx)^(-a)*ki(wx)^(-b);
Ki2 = @(wx) ki(wx)*h2(wx);
Kc2 = @(wx) kc(wx)*h2(wx);

%Step 3
h1 = @(wx) (1/chi*wx/Ki(wx) + gc/Ki_Kc)/biggamc*Ki(wx)^(1-gam)*kc(wx)^(-a)*ki(wx)^(-b);
h1_check = @(wx) (wx/chi + gc/Ki_Kc*Ki(wx))/biggamc*Ki(wx)^(-gam)*kc(wx)^(-a)*ki(wx)^(-b);
check_h = (h1(150) - h1_check(150))^2;
if check > 10^(-20)
      error('h1 is wrong')
end
Ki1 = @(wx) ki(wx)*h1(wx);
Kc1 = @(wx) kc(wx)*h1(wx);

%Error is above this point

%Step 4
h = @(wx) h1(wx) + h2(wx);
Kc = @(wx) Kc1(wx) + Kc2(wx);
yc = @(wx) biggamc*Ki(wx)^(gam)*h1(wx)^(1-a-b)*Kc1(wx)^(a)*Ki1(wx)^(b);
yi = @(wx) biggami*Ki(wx)^(gam)*h2(wx)^(1-a-b)*Kc2(wx)^(a)*Ki2(wx)^(b);
Ii = @(wx) gi*Ki(wx);
Ic = @(wx) gc*Kc(wx);
c = @(wx) yc(wx) - Ic(wx);

%Step 5 - Minimization
%Objective: w = @(wx) chi*c(wx)
options = optimoptions('fmincon'); 
% Set OptimalityTolerance to 1e-15
options = optimoptions(options, 'OptimalityTolerance', 1e-23); 
% Set the Display option to 'iter' and StepTolerance to 1e-4
% options.Display = 'iter';
options.StepTolerance = 1e-30;
objw = @(wx) (wx - chi*c(wx))^2;
wx0 = 1.5;
wstar = fmincon(objw,wx0,[],[],[],[],[],[],[],options);
% Try to use also fzero
% objw2 = @(wx)  wx - chi*c(wx);
% wstar2 = fzero(objw2,wx0);
% res_obj_wstar = objw2(wstar2);

%Step 5
kc = Kc(wstar); %small kc back to be original Kc - drop previous notation
ki = Ki(wstar); %small ki back to be original Ki - drop previous notation
yc = yc(wstar);
yi = yi(wstar);
c = c(wstar);
ic = Ic(wstar);
it = Ii(wstar);
w = wstar;
h = h(wstar);
h1 = h1(wstar);
h2 = h2(wstar);
kc1 = Kc1(wstar); %small kc1 back to be original Kc1 - drop previous notation
kc2 = Kc2(wstar); %small kc2 back to be original Kc2 - drop previous notation
ki1 = Ki1(wstar); %small ki1 back to be original Ki1 - drop previous notation
ki2 = Ki2(wstar); %small ki2 back to be original Ki2 - drop previous notation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECOND STRATEGY

% %Step 1
% w = @(Kix) (  (1-a-b)^(1-a-b) * biggamc * Kix^gam * (a/rc)^a * (b/ri)^b  )^(1/(1-a-b));
% kc = @(Kix) a/(1-a-b) * w(Kix)/rc; %kc = kc1 = kc2 = Kc/h = Kc1/h1 = Kc2/h2
% ki = @(Kix) b/(1-a-b) * w(Kix)/ri; %ki = ki1 = ki2 = Ki/h = Ki1/h1 = Ki2/h2
% Ki_Kc = b/a*rc/ri; %Ki_Kc = Ki/Kc = Ki1/Kc1 = Ki2/Kc2 = ki/kc
% 
% %Step 2
% h2 = @(Kix) gi/biggami * Kix^(1-gam) * kc(Kix)^(-a) * ki(Kix)^(-b);
% Ki2 = @(Kix) ki(Kix)*h2(Kix);
% Kc2 = @(Kix) kc(Kix)*h2(Kix);
% 
% %Step 3
% h1 = @(Kix) (w(Kix)/chi + gc/Ki_Kc*Kix)/biggamc*Kix^(-gam)*kc(Kix)^(-a)*kc(wx)^(-b);


%Put the ss values in a vector consistent with Y and X vectors in model.m
% KC KI BIGGAMC BIGGAMI
xx  = [kc ki biggamc biggami]; 
yy  = [yc yi c ic it w rc ri h h1 h2 kc1 kc2 ki1 ki2 p expgc expgi];
% YC YI C IC IT W RC RI H H1 H2 KC1 KC2 KI1 KI2 P


ss  = [yy xx];



