function F = oulton_spillover_stst(param)

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