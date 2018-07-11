clear 
close all

b = 0.05;
a = 1/3 - b;
gam = 0; 
gc_mean = 0.002942;
gp_mean = - 0.0126;
w1 = (1-a-gam)/(1-a-b-gam);
w2 = (b+gam)/(1-a-b-gam);


gamC = (gc_mean - w2*gp_mean)/(w1+w2);
gamI = gamC - gp_mean;

deltyear = 0.15;
deltquarter= 1 - (1 - deltyear)^(1/4)
