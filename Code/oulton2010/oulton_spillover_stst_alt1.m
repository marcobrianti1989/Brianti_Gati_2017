function ss = oulton_spillover_stst_alt1(param)
% MIRRORS model_spillover_ss_news.m BUT IMPLEMENTS LSQNONLIN AS A SOLVER
% set global params
global bet a b biggami biggamc di dc chi gam siggami sige sigitlev rhoitlev expgc expgi gc gi
% set global st.st. values
global yc yi h1 h2 h kc1 kc2 kc ki1 ki2 ki ic it c w p
global rc ri % these two are given in closed form.
% set global functions
global eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10 
global eq13 eq14 eq15 eq16 eq17 eq18

% This converts the parameters from the struct object to one big cell.
struct_params = structvars(param);
for i = 1:size(struct_params,1)
    eval(struct_params(i,:)); % here they are each saved as a double
end

%Step 0 - define what we already know
p        = (biggamc)/(biggami); % maybe we need to put p amongst the stuff to solve for
expgc    = biggamc^((1-b-gam)/(1-a-b-gam))*biggami^((b+gam)/(1-a-b-gam)); % Roman I, p.2 of my notes4, in exp().
expgi    = biggamc^(a/(1-a-b-gam))*biggami^((1-a)/(1-a-b-gam)); % Roman II, p.3 of my notes4, in exp().
gc       = expgc - (1-dc); 
gi       = expgi - (1-di);
rc       = 1/bet * expgc - (1-dc); % eq11
ri       = (1/bet* expgi -(1-di))*p; % eq12

%Step 1
kc_bar          = @(wx) wx/rc*a/(1-a-b); %kc = kc1 = kc2 = Kc/h = Kc1/h1 = Kc2/h2
ki_bar          = @(wx) wx/ri*b/(1-a-b); %ki = ki1 = ki2 = Ki/h = Ki1/h1 = Ki2/h2 
ki              = @(wx) (rc/(a*biggamc)*kc_bar(wx)^(1-a)*ki_bar(wx)^(-b))^(1/gam);
Ki_check        = @(wx)(ri/(b*biggamc)*kc_bar(wx)^(-a)*ki_bar(wx)^(1-b))^(1/gam); 
Ki_check2       = @(wx) ( wx/((1-a-b)*biggamc)*kc_bar(wx)^(-a)*ki_bar(wx)^(-b) )^(1/gam); %check
check           = ((ki(15) - Ki_check(15)) + (ki(15) - Ki_check2(15)))^2;
if check > 10^(-16)
      warning('Ki is wrong')
end
Ki_Kc = b/a*rc/ri; %Ki_Kc = Ki/Kc = Ki1/Kc1 = Ki2/Kc2 = ki/kc

%Step 2
h2 = @(wx) gi/biggami*ki(wx)^(1-gam)*kc_bar(wx)^(-a)*ki_bar(wx)^(-b);
ki2 = @(wx) ki_bar(wx)*h2(wx);
kc2 = @(wx) kc_bar(wx)*h2(wx);

%Step 3
% Case of V(h) = -chi*H
% h1 = @(wx) (1/chi*wx/ki(wx) + gc/Ki_Kc)/biggamc*ki(wx)^(1-gam)*kc_bar(wx)^(-a)*ki_bar(wx)^(-b);
% h1_check = @(wx) (wx/chi + gc/Ki_Kc*ki(wx))/biggamc*ki(wx)^(-gam)*kc_bar(wx)^(-a)*ki_bar(wx)^(-b);
% % Case of V(h) = chi*log(1-h)
% h1 = @(wx) (wx*h2(wx)/chi - wx/chi + gc/Ki_Kc*ki(wx)) / (biggamc*ki(wx)^(gam)*kc_bar(wx)^(a)*ki_bar(wx)^(b) -wx/chi);
% h1_check = @(wx) (biggamc*ki(wx)^(gam)*kc_bar(wx)^(a)*ki_bar(wx)^(b) -wx/chi)^(-1) ...
%     *(gc/Ki_Kc*ki(wx) + biggamc*ki(wx)^(gam)*kc_bar(wx)^(a)*ki_bar(wx)^(b)*h2(wx) -wx/chi) - h2(wx);
% Case of V(h) = chi/2*H^2
h1 = @(wx) (gc/Ki_Kc*ki(wx) + biggamc*ki(wx)^(gam)*kc_bar(wx)^(a)*ki_bar(wx)^(b)*h2(wx)  ...
    + sqrt(  (gc/Ki_Kc*ki(wx) + biggamc*ki(wx)^(gam)*kc_bar(wx)^(a)*ki_bar(wx)^(b)*h2(wx))^2  ...
    + 4*biggamc*ki(wx)^(gam)*kc_bar(wx)^(a)*ki_bar(wx)^(b) *wx/chi)) ...
    / (2*biggamc*ki(wx)^(gam)*kc_bar(wx)^(a)*ki_bar(wx)^(b)) -h2(wx);
aa = @(wx) biggamc * ki(wx)^gam * kc_bar(wx)^a * ki_bar(wx)^b; % my notes4, p. 54
bb = @(wx) biggamc * ki(wx)^gam * kc_bar(wx)^a * ki_bar(wx)^b * h2(wx) - Ki_Kc^(-1)*ki(wx)*gc;
cc = @(wx) -(wx/chi + Ki_Kc^(-1)*ki(wx)*gc * h2(wx));
h1_check = @(wx) (-bb(wx) + sqrt(bb(wx)^2 -4*aa(wx)*cc(wx)))/(2*aa(wx)); % taking the positive root
check_h = (h1(150) - h1_check(150))^2;
if check_h > 10^(-14)
      warning('h1 is wrong')
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

OPTIONS = optimoptions('lsqnonlin'); % Note: for all 3 options, 1e-10 is the default value
OPTIONS.OptimalityTolerance = 1e-15;
OPTIONS.StepTolerance = 1e-25;
OPTIONS.FunctionTolerance = 1e-40;

x0 = 100;
lb = [];
ub = [];

[x,res] = lsqnonlin(@problem,x0,lb,ub,OPTIONS);

wstar  = x(1);

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


check1 = yc - biggamc * h1^(1-a-b) * kc1^a * ki1^b;
check2 = yi - biggami * h2^(1-a-b) * kc2^a * ki2^b;
check3 = kc -kc1 -kc2;
check4 = ki -ki1 -ki2;
check5 =  h -h1 -h2;
check6 = kc * gc - ic;
check7 = ki * gi - it;
check8 = yc - c -ic;
check9 = yi -it;
check10 = w/c -chi*h;
check13 =  w -(1-a-b)*biggamc * h1^(-a-b) * kc1^(a) * ki1^(b);
check14  =  rc -(a)*biggamc * h1^(1-a-b) * kc1^(a-1) * ki1^(b);
check15  = ri -(b)*biggamc * h1^(1-a-b) * kc1^(a) * ki1^(b-1);
check16  = w -(1-a-b)*biggami * h2^(-a-b) * kc2^(a) * ki2^(b);
check17  = rc -(a)*biggami * h2^(1-a-b) * kc2^(a-1) * ki2^(b);
check18  = ri -(b)*biggami * h2^(1-a-b) * kc2^(a) * ki2^(b-1);

check = [check1 check2 check3 check4 check5 check6 check7 check8 check9 check10 check13 check14 check15 check16 check17 check18];

if sum(check.^2) > 10^(-12)
    warning('Steady state solution not reached.')
%     error('Steady state solution not reached.')
end
sum(check.^2)

%Put the ss values in a vector consistent with Y and X vectors in model.m
xx  = [kc ki biggamc biggami ...
    c biggamc biggami yc yi h p kc2 ki2 ki kc]; 
yy  = [yc yi c ic it w rc ri h h1 h2 kc1 kc2 ki1 ki2 p expgc expgi ...
    expgc expgi expgc expgi 1 p expgc expgi];

ss  = [yy xx];
end

function out = problem(X)

global chi h c

wx  = X(1);

out(1,1)     = chi*h(wx) - wx/c(wx);


end
