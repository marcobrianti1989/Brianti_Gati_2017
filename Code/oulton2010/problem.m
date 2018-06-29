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

out(1,1) = fun1(s,phimax,phimin,L,thet);
out(2,1) = fun2(s,phimax,phimin,L,thet);
out(3,1) = fun3(s,phimax,phimin,L,thet);
out(4,1) = fun4(s,phimax,phimin,L,thet);
out(5,1) = fun5(s,phimax,phimin,L,thet);

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