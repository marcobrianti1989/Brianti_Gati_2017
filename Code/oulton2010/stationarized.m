%stationarized System. My notes p. 228
f(1)    = -YC + (1+GAMC)*H1^(1-a-b)*KC1^(a)*KI1^(b);
f(end+1)= -YI + (1+GAMI)*H2^(1-a-b)*KC2^(a)*KI2^(b);
f(end+1)= -KC + KC1 + KC2;
f(end+1)= -KI + KI1 + KI2;
f(end+1)= -H + H1 + H2;
f(end+1)= -KC_p*(1+G) + (1-dc)*KC + IC;
f(end+1)= -KI_p*(1+GI)+ (1-di)*KI + IT; 
f(end+1)= -YC + C + IC;
f(end+1)= -YI + IT;
f(end+1)= -W + C*chi;
f(end+1)= -1 + bet*C/C_p*1/(1+G)*(RC_p + 1-dc);
f(end+1)= -1 + bet*C/C_p*1/(1+G)*(RI_p/(1+GI_p) + 1-di);
f(end+1)= -W + (1-a-b)*(1+GAMC)*H1^(-a-b)*KC1^(a)*KI1^(b);
f(end+1)= -RC + a*(1+GAMC)*H1^(1-a-b)*KC1^(a-1)*KI1^(b);
f(end+1)= -RI + b*(1+GAMC)*H1^(1-a-b)*KC1^(a)*KI1^(b-1);
f(end+1)= -W + (1-a-b)*(1+GAMI)*(1+GP)*H2^(-a-b)*KC2^(a)*KI2^(b);
f(end+1)= -RC + a*(1+GAMI)*(1+GP)*H2^(1-a-b)*KC2^(a-1)*KI2^(b);
f(end+1)= -RI + b*(1+GAMI)*(1+GP)*H2^(1-a-b)*KC2^(a)*KI2^(b-1);
f(end+1)= % GAMC
f(end+1)= % GAMI
f(end+1)= -GP + GAMC -GAMI; % eq. I p. 222
f(end+1)= -G + ((1-b)*GAMC + b*GAMI)/(1-a-b);  % eq. II
f(end+1)= -GI + ((a)*GAMC + (1-a)*GAMI)/(1-a-b);  % eq. III

% steady state, starting at my notes p. 231, green equations.
% need gamc and gami
gp = gamc - gami; % eq. I p. 222
g = ((1-b)*gamc + b*gami)/(1-a-b);  % eq. II
gi = ((a)*gamc + (1-a)*gami)/(1-a-b);  % eq. III
rc = (1+g)/bet - (1-dc); % green 1
ri = ((1+g)/bet - (1-di))*(1+gp); % green 2
ki1_ki2 = (1+gami)/(1+gi -(1-di))*(ri)/(b*(1+gp)) -1; % green 9
kc1_kc2 = ki1_ki2; % green 8
h1_h2   = ki1_ki2:  % green 8
w_kc1   = chi*rc/a - chi*(1+g -(1-dc))*(1+kc1_kc2^(-1)); % green 10
h1      = w_kc1 * a/(1-a-b)*1/rc; % green 11
h2      = h1/h1_h2;
kc1     = ((a/rc)^(1-b) * (b/ri)^b *(1+gamc))^(1/(1-a-b))*h1; % green 14
kc2     = kc1/kc1_kc2;
ki1_h1  = (b*(1+gamc)/ri)^(1/(1-b))*(kc1/h1)^(a/(1-b)); % green 13
ki1     = ki1_h1*h1;
ki2     = ki1/ki1_ki2;
w       = w_kc1 *kc1; 