function ss = oulton_spillover_stst_alt2(param)
% DOES A SMALLER EQ. SYSTEM THAT'S BEEN SIMPLIFIED SOMEWHAT - DOESN'T SOLVE
% B/C MORE UNKNOWNS THAN EQUATIONS
% set global params
global bet a b biggami biggamc di dc chi gam siggami sige sigitlev rhoitlev expgc expgi gc gi p
% set global st.st. values
global yc yi h1 h2 h kc1 kc2 kc ki1 ki2 ki ic it c w 
global kc_h ki_h
global rc ri % these two are given in closed form.
% set global functions
global eq1 eq2 
global eq13 eq14 eq15 

% This converts the parameters from the struct object to one big cell.
struct_params = structvars(param);
for i = 1:size(struct_params,1)
    eval(struct_params(i,:)); % here they are each saved as a double
end

%Step 0 - define what we already know (eq numbers correspond to my notes4, p. 6.)
p        = (biggamc)/(biggami); % eq (13) and (18)
expgc    = biggamc^((1-b-gam)/(1-a-b-gam))*biggami^((b+gam)/(1-a-b-gam)); % See stationarized_system_spillover.pdf
expgi    = biggamc^(a/(1-a-b-gam))*biggami^((1-a)/(1-a-b-gam)); % See stationarized_system_spillover.pdf
gc       = expgc - (1-dc); % See stationarized_system_spillover.pdf
gi       = expgi - (1-di); % See stationarized_system_spillover.pdf
rc       = 1/bet * expgc - (1-dc); % eq11
ri       = (1/bet* expgi -(1-di))*p; % eq12

% Step 2 - define all the remaining st.st. equations as global functions 
% (eq numbers correspond to my notes4, p. 6. and stationarized_system_spillover.pdf)
eq1  = @(w, kc_h, ki_h, ki, h1, h2) h1-(w/chi + gc*(ki_h/kc_h)^(-1) *ki)*1/biggamc * ki^(-gam)*kc_h^(-a)*ki_h^(-b); % eq. 29 in the pdf.
eq2  = @(w, kc_h, ki_h, ki, h1, h2) h2-gi/biggami*ki^(1-gam)*kc_h^(-a)*ki_h^(-b); % eq. 27 in the pdf.
eq13 = @(w, kc_h, ki_h, ki, h1, h2) w -(1-a-b)*biggamc * kc_h^(a) * ki_h^(b) * ki^gam;
eq14 = @(w, kc_h, ki_h, ki, h1, h2) rc -(a)*biggamc * kc_h^(a-1) * ki_h^(b) * ki^gam;
eq15 = @(w, kc_h, ki_h, ki, h1, h2) ri -(b)*biggamc * kc_h^(a) * ki_h^(b-1) * ki^gam;
% 5 eqs

OPTIONS = optimoptions('lsqnonlin'); % Note: for all 3 options, 1e-10 is the default value
OPTIONS.OptimalityTolerance = 1e-15;
OPTIONS.StepTolerance = 1e-25;
OPTIONS.FunctionTolerance = 1e-40;

%w, kc_h, ki_h, ki, h1, h2
%     w    kc_h  ki_h   ki    h1   h2   
x0 = [0.6, 1.1,  1.05,  0.6, 0.8, 0.1];
% lb = [0,0,0,0,0.00001];
% ub = [1,1,1,1,1000000];
% x0 = 1*ones(1,6);
lb = zeros(1,6);
ub = 1000*ones(1,6);

[x,res] = lsqnonlin(@problem,x0,lb,ub,OPTIONS);

w     = x(1);
kc_h  = x(2);
ki_h  = x(3);
ki    = x(4);
h1    = x(5);
h2    = x(6);

check1  = h1-(w/chi + gc*(ki_h/kc_h)^(-1) *ki)*1/biggamc * ki^(-gam)*kc_h^(-a)*ki_h^(-b); 
check2  = h2-gi/biggami*ki^(1-gam)*kc_h^(-a)*ki_h^(-b);
check13 =  w -(1-a-b)*biggamc * h1^(-a-b) * kc1^(a) * ki1^(b);
check14  =  rc -(a)*biggamc * h1^(1-a-b) * kc1^(a-1) * ki1^(b);
check15  = ri -(b)*biggamc * h1^(1-a-b) * kc1^(a) * ki1^(b-1);


check = [check1 check2 check13 check14 check15];

if sum(check.^2) > 10^(-12)
    warning('Steady state solution not reached.')
%     error('Steady state solution not reached.')
end
sum(check.^2)

% %Put the ss values in a vector consistent with Y and X vectors in model.m
% xx  = [kc ki biggamc biggami ...
%     c biggamc biggami yc yi h p kc2 ki2 ki kc]; 
% yy  = [yc yi c ic it w rc ri h h1 h2 kc1 kc2 ki1 ki2 p expgc expgi ...
%     expgc expgi expgc expgi 1 p expgc expgi];
% 
% ss  = [yy xx];
ss = [w, kc_h, ki_h, ki, h1, h2];
end

function out = problem(X)

global eq1 eq2 
global eq13 eq14 eq15

w     = X(1);
kc_h  = X(2);
ki_h  = X(3);
ki    = X(4);
h1    = X(5);
h2    = X(6); 

out(1,1)     = eq1(w, kc_h, ki_h, ki, h1, h2);
out(end+1,1) = eq2(w, kc_h, ki_h, ki, h1, h2);
out(end+1,1) = eq13(w, kc_h, ki_h, ki, h1, h2);
out(end+1,1) = eq14(w, kc_h, ki_h, ki, h1, h2);
out(end+1,1) = eq15(w, kc_h, ki_h, ki, h1, h2);

end
