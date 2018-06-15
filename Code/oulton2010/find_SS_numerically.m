function output = find_SS_numerically(params)

global h chi xi del k bet A rho sigma_z psi_c tau b
global uu p q qp pp c cp
global fun1 fun2 fun3 fun4 fun5

bla = structvars(params);
for i1 = 1:size(bla,1)
    eval(bla(i1,:));
end

uu           = log(h+b);
p            = @(thet) thet/(1+thet^chi)^(1/chi);
q            = @(thet) 1/((1+thet^chi)^(1/chi));
qp           = @(thet) -thet^(chi-1)/((1+thet^chi)^(1/chi+1));
pp           = @(thet) q(thet) + thet*qp(thet);
c            = @(s) A*((1/(1-s))^(1+psi_c)-1)/(1+psi_c) - A*s;
cp           = @(s) A*(1/(1-s))^(2+psi_c)- A;

fun1 = @(s,phimax,phimin,L,thet) log(h+b) - log(phimin) - ...
    bet*( c(s) + (1 - s*p(thet) - del) *  cp(s)/p(thet) ); %same of latex

fun2 = @(s,phimax,phimin,L,thet) - k/q(thet) + 1/2*(1-phimax)^2 + ...
    (1-phimin)*bet*(1-del)*k/q(thet) ; %same of latex

fun3 = @(s,phimax,phimin,L,thet) - cp(s)/p(thet) + log(phimax) + phimin*(1-log(phimin)) ...
    - phimax + (1-phimin)*(- uu + bet*(c(s) +  (1 - s*p(thet) - del) *  cp(s)/p(thet)));

fun4 = @(s,phimax,phimin,L,thet) - L + ( (1-del)*L + p(thet)*s*(1-L) )*(1-phimin);

fun5 = @(s,phimax,phimin,L,thet) - xi*k*thet/phimax + (1-xi)*cp(s);

%Evaluate the steady state for a simplified version of Mitman and
%Rabinovich
OPTIONS = optimoptions('lsqnonlin');
OPTIONS.OptimalityTolerance = 1e-10;
OPTIONS.StepTolerance = 1e-10;
OPTIONS.FunctionTolerance = 1e-10;

x0 = [0.1,0.9,0.1,0.9,10];
lb = [0,0,0,0,0.00001];
ub = [1,1,1,1,1000000];
[x,res] = lsqnonlin(@problem,x0,lb,ub,OPTIONS);

s            = x(1);
phimax       = x(2);
phimin       = x(3);
L            = x(4);
thet         = x(5);

%Variable names in dynare are different
output.sss        = s;
output.wss   = phimax;
output.wminss   = phimin;
output.Lss        = L;
output.thetss     = thet;

max(abs([fun1(s,phimax,phimin,L,thet);...
    fun2(s,phimax,phimin,L,thet);...
    fun3(s,phimax,phimin,L,thet);...
    fun4(s,phimax,phimin,L,thet);...
    fun5(s,phimax,phimin,L,thet)]));

check1 = log(h+b) - log(phimin) - bet*( c(s) + (1 - s*p(thet) - del) *  cp(s)/p(thet) );
check2 = - cp(s)/p(thet) + log(phimax) - phimax + phimin*(1 - log(phimin)) ...
    + (1 - phimin)*( - log(h+b) + bet*( c(s) + (1 - s*p(thet) - del) *  cp(s)/p(thet) )  );
check3 = L - ((1-del)*L + p(thet)*s*(1-L))*(1-phimin);
check4 = xi*k*thet/phimax - (1-xi)*cp(s);
check5 = + k/q(thet) - 1/2*(1-phimax)^2 - (1-phimin)*bet*(1-del)*k/q(thet);

check = [check1 check2 check3 check4 check5];

if sum(check.^2) > 10^(-12)
    error('Steady state solution not reached.')
end
sum(check.^2)


end





