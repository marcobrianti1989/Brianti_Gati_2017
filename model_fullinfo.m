%Declare Needed Symbols
syms K K_p Z Z_p G G_p
syms C C_p R R_p I I_p PI PI_p O O_p M M_p S S_p V V_p

%Declare X and Y vectors
X  = [K   Z   G  ]; %states at time t
XP = [K_p Z_p G_p]; %states at time t+1
Y  = [C   R   I   PI   O   M   S   V]; %jumps at time t
YP = [C_p R_p I_p PI_p O_p M_p S_p V_p] ; %jumps at time t+1

%Model Equations - notice that we stationarized the economy dividing for
%its own growth rate gamma
f(1) = 1 - (beta*(C/C_p)*(R_p + 1 - delta)); %Intertemporal EE
f(2) = O - K^(1-theta)*Z*M^theta; %technology
f(3) = R - (1-theta)*Z*K^(-theta)*M^theta; %Capital demand 
f(4) = (Z_p)/(theta*S) - theta*Z*K^(1-theta)*M^(theta-1); %Intermediate input demand
f(5) = PI - ((Z_p)/(theta*S) - 1)*M; %monopolistic profits of the intermediate producer
%notice in f(5) that I am assuming that the cost of producing M is 1 and
%the price I can sell it is (Z_p)/(theta*S). Unfortunately, I need this
%strange price rule to ensure the existence of a  stationary steady state.
%This is a limit of the model since also C&G make a similar assumption in
%their more complicated version. In any case, the intuition is the
%following, as S (R&D expenditure) increases then the price decreases for
%decreasing return to scale, more competition in the sector,...However, as
%the stock of technology increases the the price will be higher, the
%intuition is that the tech level of the economy is overall higher an thus
%we can sell the intermediate good at a higher price. theta is a parameter.
f(6) = (1-delta)*K + I - K_p; %Law of motion of capital
f(7) = C + K_p + S - (1 - delta)*K - R*K - Z*PI - G; %Budget constraint
f(8) = V - PI - phi*beta*(C/C_p)*V_p; %present discount value of the flow profits of monopolist
f(9) = 1/fit - phi*beta*(C/C_p)*V_p; %free entry condition. 1/fit is the price I have to pay...
%in order to be sure to invent a new intermediary good and become a
%monopolist. fit is the efficiency of the R&D sector, it is this guy that
%cannot be observed and it will be noisy...
f(10) = Z_p - fit*S - phi*Z; %law of motion of endogenous TFP component
f(11) = log(G_p) - rhog*log(G); %process of exogenous component of TFP

%Verify if the computation of Steady-State is correct
fnum = double(subs(f, [Y X YP XP], [ss, ss])); %It is substituting the ss values into the f system.
%Morevore, it is also calculating f(i) for all i in order to check if the
%ss values are correct. If they are correct than f() should be a vector of
%zeros at working precision.
disp('Checking steady-state equations:'); %It reminds me that it is checking if the ss values are...
%correct or not
disp(fnum); %it displays me the vector f() to show me if the steady states we evaluated above are...
%correct or nor

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