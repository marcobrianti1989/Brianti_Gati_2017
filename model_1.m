% ====================SIMPLIFIED VERSION OF COMIN AND GERTLER (2006)=======================

% =================== 3 of August, 2017. Brianti and Gati - Attempt Number 1 ==================

%Code provided by Laura Gati and Marco Brianti on the basis of
%and Stephanie Schmitt-Grohé and Martín Uribe algorithm

%Clearing and closing previous projects
clear all
close all

%Parameterization
param

%Steady State
ss = model_1_ss;

%Declare Needed Symbols
syms K K_p Z Z_p 
syms C C_p R R_p PI PI_p Q Q_p M M_p S S_p V V_p LAMB LAMB_p W W_p L L_p

%Declare X and Y vectors
X  = [K   Z   ]; %states at time t
XP = [K_p Z_p ]; %states at time t+1
Y  = [Q   C   L   S   W   R   M   LAMB   PI   V  ]; %jumps at time t
YP = [Q_p C_p L_p S_p W_p R_p M_p LAMB_p PI_p V_p]; %jumps at time t+1

%Model Equations 
f(1)     = -1/C + chi*L^(thet-1)/W;
f(end+1) = -1 + bet*C/C_p*(1-delt + R_p);
f(end+1) = -C - K_p -S + W*L + (1-delt + R)*K + PI;
f(end+1) = -W + (1-gamm)*(1-alph)*Q/L;
f(end+1) = -R + (1-gamm)*alph*Q/K;
f(end+1) = -nu + gamm*Q/M;
f(end+1) = -PI + (nu-1)*M;
f(end+1) = -V + PI + phi*bet*C/C_p*V_p;
f(end+1) = -Z_p + LAMB*S + phi*Z;
f(end+1) = -1/LAMB + phi*bet*C/C_p*V_p;
f(end+1) = -Q + K^(alph*(1-gamm)) *L^((1-alph)*(1-gamm)) *Z^(nu*gamm) * M^gamm;
f(end+1) = -LAMB + lambx * K^(-rho) *S^(rho-1);


%Verify if the computation of Steady-State is correct
fnum = double(subs(f, [Y X YP XP], [ss, ss])); %It is substituting the ss values into the f system.
%Morevore, it is also calculating f(i) for all i in order to check if the
%ss values are correct. If they are correct than f() should be a vector of
%zeros at working precision.
disp('Checking steady-state equations:'); %It reminds me that it is checking if the ss values are...
%correct or not
disp(fnum); %it displays me the vector f() to show me if the steady states we evaluated above are...
%correct or nor
if sum(fnum) > 10^(-10)
    warning('The steady state solution derived analytically is wrong')
end

return
%Log-linear approx
log_var = [X Y XP YP];
f = subs(f, log_var, exp(log_var)); 
   
%Differentiate
fx  = jacobian(f,X);
fy  = jacobian(f,Y);
fxp = jacobian(f,XP);
fyp = jacobian(f,YP);

%Plug back into levels
f =   subs(f ,  log_var, log(log_var));
fx =  subs(fx , log_var, log(log_var));
fy =  subs(fy , log_var, log(log_var));
fxp = subs(fxp, log_var, log(log_var));
fyp = subs(fyp, log_var, log(log_var));

%Compute numerical values
fxn =  double(subs(fx , [Y X YP XP], [ss, ss]));
fyn =  double(subs(fy , [Y X YP XP], [ss, ss]));
fxpn = double(subs(fxp, [Y X YP XP], [ss, ss]));
fypn = double(subs(fyp, [Y X YP XP], [ss, ss]));

%Compute the transition and policy functions, using code by
%Stephanie Schmitt-Grohé and Martín Uribe (and available on their wedsite.)
[gx,hx]=gx_hx_alt(fyn,fxn,fypn,fxpn);
%Remember that x_t+1 = hx*x_t and y_t+1 = gx*x_t

%Compute the IRF to a positive technology shock
j=20; %Defining the prefered time horizon
%Defining matrices for the IRFs to a technology shock
irx = zeros(length(xx),j); %Notice that we are doing it for all states (gamma and k)
iry = zeros(length(yy),j); %Notice thet we are doing it for all variables. 
irf_y = zeros(1,j); %We add one vector to add also output which we want to see its IRF
%Evaluating the IRFs. We use hx to evaluate first the effect on the states
for i=1:j
    irx(:,i) = hx^(i-1)*[0; 0; -sigma]; %[sigma; 0] is the vector of shocks. We have just one shock...
    %the first element which is gamma. No shock to capital, it does not
    %make sense in this setting
    iry(:,i) = gx*irx(:,i); %Use gx to transform the IRF on the states into the jumps
end

% %Plot the IRF - complete IRF
% subplot(5,2,1)
% plot(iry(5,:))
% title('Output')
% grid on
% subplot(5,2,2)
% plot(iry(1,:))
% title('Consumption')
% grid on
% subplot(5,2,3)
% plot(iry(3,:))
% title('Investment')
% grid on
% subplot(5,2,4)
% plot(iry(7,:))
% title('R&D Effort')
% grid on
% subplot(5,2,5)
% plot(irx(1,:))
% title('Capital')
% grid on
% subplot(5,2,6)
% plot(irx(2,:))
% title('Endogenous TFP')
% grid on
% subplot(5,2,7)
% plot(iry(4,:))
% title('Intermediate Profits')
% grid on
% subplot(5,2,8)
% plot(iry(8,:))
% title('PDV Invention')
% grid on
% subplot(5,2,9)
% plot(iry(6,:))
% title('Intermediate Good')
% grid on
% subplot(5,2,10)
% plot(irx(3,:))
% title('Wealth Shock')
% grid on




%Plot the IRF
subplot(2,2,1)
plot(iry(3,:))
title('Investment')
grid on
subplot(2,2,2)
plot(iry(7,:))
title('R&D Effort')
grid on
subplot(2,2,3)
plot(irx(1,:))
title('Capital')
grid on
subplot(2,2,4)
plot(irx(2,:))
title('Endogenous TFP')
grid on
