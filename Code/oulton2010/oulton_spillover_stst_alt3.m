function ss = oulton_spillover_stst_alt3(known_ss_values)
% A COPY OF OULTON_SPILLOVER_STST.M, WHICH SOLVES THE ENTIRE SS SYSTEM
% NUMERICALLY. HERE THE DIFFERENCE IS THAT THE CODE TAKES SOME STEADY STATE
% VALUES AS GIVEN AND SOLVES FOR PARAMETERS INSTEAD, NOT THE OPPOSITE. 
% set global params
global bet a b biggami biggamc di dc chi gam siggami sige sigitlev rhoitlev expgc expgi gc gi p
% set global st.st. values
global yc yi h1 h2 h kc1 kc2 kc ki1 ki2 ki ic it c w rc ri %kcbar kibar
% set global functions
global eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10 eq11 eq12
global eq13 eq14 eq15 eq16 eq17 eq18

% See my notes4, p. 59.
ic = known_ss_values(1);
it = known_ss_values(2);
c  = known_ss_values(3);
w  = known_ss_values(4);
rc = known_ss_values(5);
ri = known_ss_values(6);
p  = known_ss_values(7);
h  = known_ss_values(8);
kc = known_ss_values(9);

%Step 0 - define what we already know. Blue stuff in my notes.
yc = c + ic;
yi = it;
chi = w/c;
bet = 0.97; 
gam = 0; % for now.


% Step 2 - define all the remaining st.st. equations as global functions (eq numbers correspond to my notes3, p. 183.)
eq1  = @(biggamc,biggami,ki,kcbar,kibar,a,b, h1, h2,expgc,expgi,dc,di) yc - biggamc * ki^gam * h1^(1) * kcbar^a * kibar^b;
eq2  = @(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di) yi - biggami * ki^gam * h2^(1) * kcbar^a * kibar^b;
eq3  = @(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di) kc * (expgc -1+dc) - ic;
eq4  = @(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di) ki * (expgi -1+di) - it;
eq10  = @(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di) h -h1 -h2;
eq11  = @(biggamc,biggami,ki,kcbar,kibar,a,b,h1,h2,expgc,expgi,dc,di) rc - 1/bet * expgc - (1-dc);
eq12  = @(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di) ri - (1/bet* expgi -(1-di))*p;
eq13  = @(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di) w -(1-a-b)*biggamc * ki^gam  * kcbar^(a) * kibar^(b);
eq14  = @(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di) rc -(a)*biggamc * ki^gam * kcbar^(a-1) * kibar^(b);
eq15  = @(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di) ri -(b)*biggamc * ki^gam * kcbar^(a) * kibar^(b-1);
eq16  = @(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di) p-biggamc/biggami;
eq17  = @(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di) expgc - biggamc^((1-b-gam)/(1-a-b-gam))*biggami^((b+gam)/(1-a-b-gam));
eq18  = @(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di) expgi - biggamc^(a/(1-a-b-gam))*biggami^((1-a)/(1-a-b-gam));
% 13 eqs, 13 unknowns

OPTIONS = optimoptions('lsqnonlin'); % Note: for all 3 options, 1e-10 is the default value
OPTIONS.OptimalityTolerance = 1e-15;
OPTIONS.StepTolerance = 1e-25;
OPTIONS.FunctionTolerance = 1e-40;

%        yc   yi   h1   h2   h    kc1  kc2  kc    ki1  ki2   ki   ic    it   c    w 
% x0 = [0.69, 0.1, 0.8, 0.1, 0.9, 0.7, 0.1, 0.79, 0.5, 0.1, 0.6, 0.06, 0.1, 0.6, 0.6];
% lb = [0,0,0,0,0.00001];
% ub = [1,1,1,1,1000000];
x0 = 1*ones(1,13);
lb = zeros(1,13);
ub = 1000*ones(1,13);

[x,res] = lsqnonlin(@problem,x0,lb,ub,OPTIONS);

biggamc  = x(1);
biggami  = x(2);
ki       = x(3);
kcbar    = x(4);
kibar    = x(5);
a        = x(6);
b        = x(7);
h1       = x(8);
h2       = x(9);
expgc    = x(10);
expgi    = x(11);
dc       = x(12);
di       = x(13);

check1  = yc - biggamc * ki^gam * h1^(1) * kcbar^a * kibar^b;
check2  = yi - biggami * ki^gam * h2^(1) * kcbar^a * kibar^b;
check3  = kc * (expgc -1+dc) - ic;
check4  = ki * (expgi -1+di) - it;
check10  = h -h1 -h2;
check11  = rc - 1/bet * expgc - (1-dc);
check12  = ri - (1/bet* expgi -(1-di))*p;
check13  = w -(1-a-b)*biggamc * ki^gam  * kcbar^(a) * kibar^(b);
check14  = rc -(a)*biggamc * ki^gam * kcbar^(a-1) * kibar^(b);
check15  = ri -(b)*biggamc * ki^gam * kcbar^(a) * kibar^(b-1);
check16  = p-biggamc/biggami;
check17  = expgc - biggamc^((1-b-gam)/(1-a-b-gam))*biggami^((b+gam)/(1-a-b-gam));
check18  = expgi - biggamc^(a/(1-a-b-gam))*biggami^((1-a)/(1-a-b-gam));

check = [check1 check2 check3 check4 check10 check11 check12 check13 check14 check15 check16 check17 check18];

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

global eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10 eq11 eq12
global eq13 eq14 eq15 eq16 eq17 eq18

biggamc  = X(1);
biggami  = X(2);
ki       = X(3);
kcbar    = X(4);
kibar    = X(5);
a        = X(6);
b        = X(7);
h1       = X(8);
h2       = X(9);
expgc    = X(10);
expgi    = X(11);
dc       = X(12);
di       = X(13);

out(1,1)     = eq1(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di);
out(end+1,1) = eq2(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di);
out(end+1,1) = eq3(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di);
out(end+1,1) = eq4(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di);
out(end+1,1) = eq10(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di);
out(end+1,1) = eq11(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di);
out(end+1,1) = eq12(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di);
out(end+1,1) = eq13(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di);
out(end+1,1) = eq14(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di);
out(end+1,1) = eq15(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di);
out(end+1,1) = eq16(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di);
out(end+1,1) = eq17(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di);
out(end+1,1) = eq18(biggamc,biggami,ki,kcbar,kibar,a,b, h1,h2,expgc,expgi,dc,di);


end
