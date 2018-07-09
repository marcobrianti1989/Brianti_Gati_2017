function ss = oulton_spillover_stst(param)
% SOLVES THE ENTIRE SS SYSTEM NUMERICALLY. IS RELATIVELY IMPRECISE; IT
% MAKES SENSE TO REDUCE THE SYSTEM.
% set global params
global bet a b biggami biggamc di dc chi gam siggami sige sigitlev rhoitlev expgc expgi gc gi p
% set global st.st. values
global yc yi h1 h2 h kc1 kc2 kc ki1 ki2 ki ic it c w 
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

% Step 2 - define all the remaining st.st. equations as global functions (eq numbers correspond to my notes3, p. 183.)
eq1  = @(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w) yc - biggamc * ki^gam * h1^(1-a-b) * kc1^a * ki1^b;
eq2  = @(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w) yi - biggami * ki^gam * h2^(1-a-b) * kc2^a * ki2^b;
eq3  = @(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w) kc -kc1 -kc2;
eq4  = @(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w) ki -ki1 -ki2;
eq5  = @(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w) h -h1 -h2;
eq6  = @(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w) kc * gc - ic;
eq7  = @(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w) ki * gi - it;
eq8  = @(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w) yc - c -ic;
eq9  = @(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w) yi -it;
eq10 = @(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w) w/c -chi;
eq13  = @(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w) ...
    w -(1-a-b)*biggamc * ki^gam * h1^(-a-b) * kc1^(a) * ki1^(b);
eq14  = @(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w) ...
    rc -(a)*biggamc * ki^gam * h1^(1-a-b) * kc1^(a-1) * ki1^(b);
eq15  = @(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w) ...
    ri -(b)*biggamc * ki^gam * h1^(1-a-b) * kc1^(a) * ki1^(b-1);
eq16  = @(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w) ...
    w/p -(1-a-b)*biggami * ki^gam * h2^(-a-b) * kc2^(a) * ki2^(b);
eq17  = @(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w) ...
    rc/p -(a)*biggami * ki^gam * h2^(1-a-b) * kc2^(a-1) * ki2^(b);
eq18  = @(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w) ...
    ri/p -(b)*biggami * ki^gam * h2^(1-a-b) * kc2^(a) * ki2^(b-1);
% 16 eqs

OPTIONS = optimoptions('lsqnonlin'); % Note: for all 3 options, 1e-10 is the default value
OPTIONS.OptimalityTolerance = 1e-15;
OPTIONS.StepTolerance = 1e-25;
OPTIONS.FunctionTolerance = 1e-40;

%        yc   yi   h1   h2   h    kc1  kc2  kc    ki1  ki2   ki   ic    it   c    w 
x0 = [0.69, 0.1, 0.8, 0.1, 0.9, 0.7, 0.1, 0.79, 0.5, 0.1, 0.6, 0.06, 0.1, 0.6, 0.6];
% lb = [0,0,0,0,0.00001];
% ub = [1,1,1,1,1000000];
% x0 = 1*ones(1,15);
lb = zeros(1,15);
ub = 1000*ones(1,15);

[x,res] = lsqnonlin(@problem,x0,lb,ub,OPTIONS);

yc  = x(1);
yi  = x(2);
h1  = x(3);
h2  = x(4);
h   = x(5);
kc1 = x(6);
kc2 = x(7);
kc  = x(8);
ki1 = x(9);
ki2 = x(10);
ki  = x(11);
ic  = x(12);
it  = x(13);
c   = x(14);
w   = x(15); 


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

global eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10 
global eq13 eq14 eq15 eq16 eq17 eq18

yc  = X(1);
yi  = X(2);
h1  = X(3);
h2  = X(4);
h   = X(5);
kc1 = X(6);
kc2 = X(7);
kc  = X(8);
ki1 = X(9);
ki2 = X(10);
ki  = X(11);
ic  = X(12);
it  = X(13);
c   = X(14);
w   = X(15); 

out(1,1)     = eq1(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w);
out(end+1,1) = eq2(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w);
out(end+1,1) = eq3(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w);
out(end+1,1) = eq4(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w);
out(end+1,1) = eq5(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w);
out(end+1,1) = eq6(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w);
out(end+1,1) = eq7(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w);
out(end+1,1) = eq8(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w);
out(end+1,1) = eq9(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w);
out(end+1,1) = eq10(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w);
out(end+1,1) = eq13(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w);
out(end+1,1) = eq14(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w);
out(end+1,1) = eq15(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w);
out(end+1,1) = eq16(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w);
out(end+1,1) = eq17(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w);
out(end+1,1) = eq18(yc, yi, h1, h2, h, kc1, kc2, kc, ki1, ki2, ki, ic, it, c, w);


end
