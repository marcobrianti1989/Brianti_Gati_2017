%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
path(path,'C:\EmpiricalMacro')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,B]=xlsread('C:\EmpiricalMacro\ZLB.xls','UnitedStates','A95:J793');
[T,N]=size(X);
Time=X(2:T,1);
LogCPI=log(X(:,2));
UnemploymentRate=X(2:T,4);
FEDFUNDS=X(2:T,5);
GS10=X(2:T,10);
CPIInflation=((diff(LogCPI)+1).^12-1)*100;
Spread=GS10-FEDFUNDS;
% The matrix of series with which we want to estimate the Bayesian VAR
X=[FEDFUNDS Spread CPIInflation UnemploymentRate]; 
[T,N]=size(X);
%
Dates=[1957+4/12 1971+8/12 1979+8/12 2007+8/12 2008+9/12]';
%
T=length(Time);
for hh=1:length(Dates)
    Regimes(:,hh)=(Time>=Dates(hh))*(-1)^hh;
end
if size(Dates,1)==1
    Regimes=((Regimes')'*2+1)*1000000;
else
    Regimes=(sum(Regimes')'*2+1)*1000000;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumberOfDraws=2000; % The number of draws from the posterior distribution
Intercept='Y'; % 'Y' if you want an intercept, 'N' otherwise
Stationary='Y'; % 'Y' if you want to impose stationarity, 'N' otherwise
Trend=0; % This is the order of the deterministic time trend to be included in the VAR.
LagOrder=12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[B,Sigma,U,TimeTrend]=BayesianVARp(X,LagOrder,Intercept,Trend,NumberOfDraws,Stationary);
if Trend>0
    CoefficientsOnTimeTrend=B(:,size(B,2),:);
    B=B(:,1:size(B,2)-1,:);
end
if Intercept=='Y'
    MU=B(:,1,:);
    B=B(:,2:size(B,2),:);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        The matrices where we store the results we are interested in:
%
% These are the signs we want to impose: 
% 1  means 'non-negative'
% -1 means 'non-positive'
% -9999 means 'I don't know, and therefore I am not imposing any sign here ...'
%
SIGNS=[1 -9999 1 -9999; -1 -1 -9999 -9999; -1 1 1 -1; 1 -1 -1 -1];
% Here we vectorize it, so that it is easier to wortk with:
SIGNS=vec(SIGNS);
% The we only keep the entries coresponding to the signs we want to
% impose:
Indices=find(SIGNS~=-9999);
SIGNS=SIGNS(Indices);
%
% These are the structural impact matrices:
AA00=zeros(N,N,NumberOfDraws); 
% This is an indicator which is equal to 1 if we draw a matrix which
% satisfies the signs we want to impose, and it is equal to zero otherwise: 
A0IN=zeros(NumberOfDraws,1); 
% This is the long-run impact of structural shocks:
LRIM=zeros(N,N,NumberOfDraws); 
%
% What we are interested in is the signs of the colums of the structural impact matrix A0.
% Here we create indices of all the possible permutations of the columns of
% a NxN matrix (in our case, 4x4), so that for each stochastically drawn
% structural impact matrix A0, we look at all possible ways of re-ordering the 
% columns in such a way that it satisfies the signs we want to impose ... 
V=(1:1:N)';
PermutationsIndices=perms(V);
%
% Order of the variables: X=[FEDFUNDS Spread CPIInflation UnemploymentRate]; 
% Order of the shocks: Monetary, Spread, Demand, Supply
%
% Here we start the loop in which we randomly draw structural impact
% matrices A0, and we then check for their signs, in order to see whether
% the signs satisfy a the pattern we want to impose ...
%
MaxNumberTry=2000; % You'll see later on what this is ...
% Variables to be saved:
DFILE=['C:\EmpiricalMacro\USSpreadShockMonthlyData'];
varname(1,:)=['LRIM'];
varname(2,:)=['BBBB'];
varname(3,:)=['SSSS'];
varname(4,:)=['UUUU'];
varname(5,:)=['MUMU'];
varname(6,:)=['AA00'];
varname(7,:)=['A0IN'];
%
ss=1;
while ss<=NumberOfDraws
    ss
    % Here we get, for each single draw from 
    % the posterior, the elements of the VAR:
    mu=MU(:,:,ss);
    bb=B(:,:,ss);
    var=Sigma(:,:,ss);
    %
    % Here the matrix C we need to get the long-run impacts of shocks:
    if LagOrder>1
        InvC=(eye(N)-sum(reshape(bb,N,N,LagOrder),3));
    else
        InvC=(eye(N)-bb);
    end
    %
    C=inv(InvC);
    % Here we randomly draw (NxN) matrices from the Normal(0,1)
    % distribution, which we are going to use to get the randomly drawn
    % structural impact matruces A0:
    KK=normrnd(0,1,N,N,MaxNumberTry);
    %
    % Here we take the eigenvalue-eigenvector decomposition of the VAR's
    % covariance matrix, which we are going to combine with 
    [P,D]=eig(var);
    %
    % This is an index which keeps track of the fact that we have not found the
    % matrix we are looking for yet: when we find it, we set the index to
    % 1, and the lopop automaically stops ...
    IND=0; 
    %   
    xx=1;
    while IND==0 & xx<MaxNumberTry % Notice that we set a maximum number of try:
        %  if one particular draw from the posterior is 'unlucky' -- in the
        %  sense that, conditional on that draw, it takes 10 years to get
        %  an impact matrix satisfying  the signs we want to impose, we do
        %  not go on drawing for 10 years: we try up to MaxNumberTry, and i
        %  we don't get anything, we move on to the nect draw ...
        %
        K=KK(:,:,xx); % The randomly draw (NxN) matrix from the Normal(0,1) distribution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This is the Rubio-Waggoner-Zha procedure: we take the QR
        % decomposition of the random matrix K:
        [Q,R]=qr(K);
        for i=1:N;
            if R(i,i)<0
                Q(:,i)=-Q(:,i);
            end
        end
        % This is our candidate impact matrix A0:
        aa0=P*D.^0.5*Q';
        % Check out that:
        %                 aa0*aa0'-var
        % is negligible. This is indeed the case ...
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Here we consider all of the possible permutations
        % of the columns of the candidate impact matrix A0:
        for zz=1:size(PermutationsIndices,1)
            A0=aa0(:,PermutationsIndices(zz,:)');
            % Here we impose the zero restriction on impact
            Omega=atan(-A0(1,2)/A0(1,1));
            RR=[cos(Omega) sin(Omega); -sin(Omega) cos(Omega)];
            RR=[RR zeros(2,2); zeros(2,2) eye(2)];
            A0=A0*RR;
            % Here we check the signs:
            SignImpactAt0=vec(sign(A0));
            SignImpactAt0=SignImpactAt0(Indices);
            % If the signs are fine, we terminate the search:
            if any(SignImpactAt0-SIGNS)==0
                IND=1;
                break
            elseif any(-SignImpactAt0-SIGNS)==0
                A0=-A0;
                IND=1;
                break
            else
            end
        end
        xx=xx+1;
    end
    if IND==1 % If we are successful, we store the results:
        AA00(:,:,ss)=A0;
        A0IN(ss)=1;
        LRIM(:,:,ss)=C*A0;
    end
    ss=ss+1;
end
BBBB=B;
SSSS=Sigma;
UUUU=U;
if Intercept=='Y'
    MUMU=MU;
else
    MUMU=-9999;
end
save(DFILE,varname(1,:),varname(2,:),varname(3,:),varname(4,:),varname(5,:),varname(6,:),varname(7,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   Now we compute IRFs, variance decompositons, etc.
%
% load C:\EmpiricalMacro\USSpreadShockMonthlyData
%
% Here we extract the indices corresponding to the successful draws:
Indices=ExtractIndices(A0IN);
NN=length(Indices);
NumberOfSuccessfulDraws=NN;
Percentiles=fix(NN*[0.5 0.16 0.84 0.05 0.95]');
N=size(LRIM,1);
%
disp('------------------------------------------------------------------')
disp('This application is very simple, so all draws are successful:')
NumberOfSuccessfulDraws
disp('In more complex applications this is typically not the case, and')
disp('you just have to work with thise draws which have been successful ')
disp('------------------------------------------------------------------')
%
% Here we extract the successful draws. In our case it is irrelevant because
% all have been successful, but in general you have to do it:
LRIM=LRIM(:,:,Indices);
BBBB=BBBB(:,:,Indices);
SSSS=SSSS(:,:,Indices);
UUUU=UUUU(:,:,Indices);
MUMU=MUMU(:,:,Indices);
AA00=AA00(:,:,Indices);
%
% Here we compute the IRFs:
% Let's just focus on the IRFs to a monetary shock and to a pure
% compression in the spread:
Horizon=10*12; % The horizon for the IRFs, in months (10 years)
%
IRFsMonetaryShock=zeros(Horizon+1,4,NumberOfDraws);
IRFsSpreadShock=zeros(Horizon+1,4,NumberOfDraws);
%
ss=1;
while ss<=NN
    ss
    bb=B(:,:,ss);
    a0=AA00(:,:,ss);
    %
    [IRFsOnlyFirstShock,IRFsOnlySecondShock]=GetIRFsMonetaryAndSpreadShocks(bb,a0,N,LagOrder,Horizon);
    % X=[FEDFUNDS Spread CPIInflation UnemploymentRate];
    %
    % Here we convert inflation from month-on-month to the annual rate, by
    % taking 12-month rolling averages:
    IRFsOnlyFirstShock(:,3)=RollingMean([zeros(11,1); IRFsOnlyFirstShock(:,3)],12);
    IRFsOnlySecondShock(:,3)=RollingMean([zeros(11,1); IRFsOnlySecondShock(:,3)],12);
    %
    % Here we normalize the IRFs. In particular, 
    % (1) we normalize the IRFs to a monetary shock so that the median of the
    %     distribution of the impact at t=0 on the FED Funds rate is 25 basis
    %     points, and
    % (2) we normalize the IRFs to a spread shock so that the median of the
    %     distribution of the impact at t=0 on the spread is 25 basis points.
    IRFsOnlyFirstShock=(IRFsOnlyFirstShock/IRFsOnlyFirstShock(1,1))*0.25;
    IRFsOnlySecondShock=(IRFsOnlySecondShock/IRFsOnlySecondShock(1,2))*0.25;
    %
    IRFsMonetaryShock(:,:,ss)=IRFsOnlyFirstShock;
    IRFsSpreadShock(:,:,ss)=IRFsOnlySecondShock;
    %
    ss=ss+1;
end
% %
% Order of the shocks: Monetary, Spread, Demand, Supply
% Variables: X=[FEDFUNDS Spread Inflation UnemploymentRate];
%
HOR=(0:1:Horizon)';
%
xx=1;
while xx<=N
    %
    IRF=sort(squeeze(IRFsMonetaryShock(:,xx,:)),2);
    IRF=IRF(:,Percentiles);
    figure(1)
    subplot(2,N,xx)
    plot(HOR,IRF(:,2:3),'r',HOR,IRF(:,1),'k',HOR,zeros(size(HOR)),'b:','LineWidth',2)
    xlim([0 Horizon])
    if xx==1
        title('Federal funds rate')
        ylabel('Monetary policy shock')
    elseif xx==2
        title('Spread')
    elseif xx==3
        title('Inflation')
    elseif xx==4
        title('Unemployment rate')
    else
    end
    %
    IRF=sort(squeeze(IRFsSpreadShock(:,xx,:)),2);
    IRF=IRF(:,Percentiles);
    figure(1)
    subplot(2,N,N+xx)
    plot(HOR,IRF(:,2:3),'r',HOR,IRF(:,1),'k',HOR,zeros(size(HOR)),'b:','LineWidth',2)
    xlim([0 Horizon])
    if xx==1
        ylabel('Spread shock')
    end
    xx=xx+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   Exploring the role played by spread shocks
%                                   within the context of the financial crisis
%
Inflation=RollingMean([zeros(11,1); CPIInflation],12);
figure(2)
subplot(1,2,1)
plot(Time,Spread,'r',Time,Inflation,'k',Time,UnemploymentRate,'b',Time,Regimes,'b',Time,zeros(size(Time)),'k:','LineWidth',2)
axis([2000 Time(T) -2 10.5])
xlabel('Red = spread; Black = inflation; Blue = unemployment rate')
%
disp('------------------------------------------------------------------')
disp('Notice the counter-cyclicality of the spread (it closely co-moves with')
disp('the unemployment rate ...')
disp('This is an incredibly robust stylized fact, holding for different countries')
disp('time periods, and alternative monetary regimes ...')
disp('------------------------------------------------------------------')
%
disp('------------------------------------------------------------------')
disp('The spread starts increasing immediately after the beginning of the crisis,')
disp('in August 2007, and this suggests that such creeping up may have played a role')
disp('in fostering the recession and the increase in the unemployment rate ...')
disp('------------------------------------------------------------------')
%
disp('------------------------------------------------------------------')
disp('So lets see the time series of the shocks to the spread')
disp('------------------------------------------------------------------')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     Getting the identified spread shocks
StructuralShocks=zeros(size(UUUU));
% 
ss=1;
while ss<=NN
    A0=AA00(:,:,ss);
    StructuralShocks(:,:,ss)=A0\UUUU(:,:,ss);
    ss=ss+1;
end
% This is the posterior distribution of the identified spread shocks:
PosteriorDistributionOfSpreadShocks=-squeeze(StructuralShocks(2,:,:));
% This is the fraction of identified spread shocks which had been positive in each month:
FractionOfPositiveSpreadShocks=sum(PosteriorDistributionOfSpreadShocks>0,2)/NN;
TimeShocks=Time(13:T);
%
figure(3)
subplot(1,2,1)
plot(TimeShocks,FractionOfPositiveSpreadShocks,'k',Time,Regimes,'b','LineWidth',2)
axis([2004 Time(T) 0 1])
%
% Then we take a rolling mean, to smooth everything:
SmoothedFractionOfPositiveSpreadShocks=FractionOfPositiveSpreadShocks;
for tt=3:length(SmoothedFractionOfPositiveSpreadShocks)
    SmoothedFractionOfPositiveSpreadShocks(tt)=RollingMean(FractionOfPositiveSpreadShocks(tt-2:tt),3);
end
figure(3)
subplot(1,2,2)
plot(TimeShocks,SmoothedFractionOfPositiveSpreadShocks,'k',Time,Regimes,'b','LineWidth',2)
axis([2004 Time(T) 0 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('------------------------------------------------------------------')
disp('On Tuesday October 23 we ended here ....')
disp('------------------------------------------------------------------')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       Now let's perform several counterfactuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First we get, as we did last time, the structural shocks, which we are
% going to need in all of the counterfactuals we do ...
StructuralShocks=zeros(size(UUUU));
% 
ss=1;
while ss<=NN
    A0=AA00(:,:,ss);
    StructuralShocks(:,:,ss)=A0\UUUU(:,:,ss);
    ss=ss+1;
end
% Then we just look at them, to chect that there's nothing weird ...
% remember: you should always take a look at the series you work with and
% the different objects you generate in MATLAB: you typically know ex ante
% 'what they should look like', more or less, and you hav to check that
% they do exhibit strange features ...
xx=1;
while xx<=N
    figure(1)
    subplot(1,4,xx)
    plot(squeeze(StructuralShocks(xx,:,:)))
    xlim([1 size(StructuralShocks,2)])
    if xx==1
        ylabel('Identified shocks')
        title('Monetary policy shock')
    elseif xx==2
        title('Spread shock')
    elseif xx==3
        title('Demand non-policy shock')
    elseif xx==4
        title('Supply shock')
    else
    end
    xx=xx+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           I: Counterfactual in which we re-run the entire
%                          post-WWII U.S. history 'killing off' specific shocks
%
Inflation=RollingMean([zeros(11,1); CPIInflation],12);
XX=[FEDFUNDS Spread Inflation UnemploymentRate];
%
% This is the shcok we want to kill off: 
% 1 = monetary shock
% 2 = spread shock
% 3 = demand non policy shock
% 4 = supply
%
disp('---------------------------------------------')
disp('Remember to check that if you dont kill off')
disp('any shock you recover the original series ...')
disp('---------------------------------------------')
%
Shock=2;
%
CounterfactualXX=zeros(T,N,NN);
CounterfactualXXMinusActualX=zeros(T,N,NN);
%
ss=1;
while ss<=NN
    ss
    mu=MUMU(:,:,ss);
    bb=BBBB(:,:,ss);
    a0=AA00(:,:,ss);
    %
    CounterfactualX=X';
    %
    tt=LagOrder+1;
    while tt<=T
        % These are the identified shocks for month tt
        % and draw ss from the posterior dostribution:
        Shocks=StructuralShocks(:,tt-LagOrder,ss);
        % Here we kill off the shock we are not interested in:
        Shocks(Shock)=0;
        % Here we re-run history after having killed off the relevant shock:
        CounterfactualX(:,tt)=mu+bb*vec(fliplr(CounterfactualX(:,tt-LagOrder:tt-1)))+a0*Shocks;
        tt=tt+1;
    end
    % Here we store it:
    CounterfactualXX(:,:,ss)=CounterfactualX';
    CounterfactualXX(:,3,ss)=RollingMean([ones(11,1)*NaN; CounterfactualXX(:,3,ss)],12);
    CounterfactualXXMinusActualX(:,:,ss)=CounterfactualXX(:,:,ss)-XX;
    ss=ss+1;
end
% Getting the percentiles of the counterfactual series:
CounterfactualSeries=sort(CounterfactualXX,3);
PercentilesCounterfactualSeries=CounterfactualSeries(:,:,Percentiles);
%
% Getting the percentiles of the difference between counterfactual and actual series:
CounterfactualMinusActualSeries=sort(CounterfactualXXMinusActualX,3);
PercentilesCounterfactualMinusActualSeries=CounterfactualMinusActualSeries(:,:,Percentiles);
%
xx=1;
while xx<=N
    %
    Counterfactual=squeeze(PercentilesCounterfactualSeries(:,xx,:));
    %
    figure(1)
    subplot(2,N,xx)
    plot(Time,Counterfactual(:,1),'k',Time,Counterfactual(:,4:5),'r',Time,XX(:,xx),'b',Time,zeros(size(Time)),'b:')
    xlim([Time(1) Time(T)])
    if xx==1
        title('Federal funds rate')
    elseif xx==2
        title('Spread')
    elseif xx==3
        title('Inflation')
    elseif xx==4
        title('Unemployment rate')
    else
    end
    %
    CounterfactualMinusActual=squeeze(PercentilesCounterfactualMinusActualSeries(:,xx,:));
    %
    figure(1)
    subplot(2,N,N+xx)
    plot(Time,CounterfactualMinusActual(:,1),'k',Time,CounterfactualMinusActual(:,4:5),'r',Time,zeros(size(Time)),'b:')
    xlim([Time(1) Time(T)])
    %
    xx=xx+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     II: A counterfactual in which we re-run the entire post-WWII U.S. history
%                        by imposing a specific monetary policy rule over the entire sample
%
disp('----------------------------------------------------------------------------------------------')
disp('Now lets see what happens when you change the parameters')
disp('of the monetary policy rule, instead of just shutting off')
disp('the structural shocks one a  time ...')
disp('----------------------------------------------------------------------------------------------')
%
CounterfactualXX=zeros(T,N,NN);
Index=zeros(NN,1);
%
% This is the coefficient which parameterizes the increase in the
% aggressiveness against inflation of the monetary policy rule, that is:
% how much we make the reaction of the FED Funmds rate to inflation 'more
% aggressive' ...
K=0; 
%
ss=1;
while ss<=NN
    ss
    mu=MUMU(:,:,ss);
    bb=BBBB(:,:,ss);
    a0=AA00(:,:,ss);
    % First, we need to get the inverse of the matrix A0:
    inva0=inv(a0);
    % Then we get the structural form of the SVAR, by pre-multiplying
    % everything by the inverse of A0:
    structural_mu=inva0*mu;
    structural_bb=inva0*bb;
    original_structural_bb=structural_bb;
    % The structure of the structural form of the SVAR is:
    % 
    % inva0 * Y(t) = structural_mu + structural_bb * [lags of Y(t)] + StructuralShocks(t)
    %
    % so considering a more aggressive monetary rule against inflation
    % means 'increasing the coefficients on lagged inflation in the monetary policy rule'
    % The order of the variables is: 
    % 
    %                  X=[FEDFUNDS Spread CPIInflation UnemploymentRate]; 
    %
    % so for each of the AR matrices of the VAR we want to increase the
    % coefficient on lagged inflation in that matrix. 
    % First of all, we need to identify the monetary policy rule in the
    % structural form of the SVAR ... 
    % How do we do that? The monetary rule is uniqeuly identified by the
    % monetary policy shock: the order of the shocks is:
    %
    %                        [Monetary Spread Demand Supply]
    %
    % and so the monetary policy rule is the fiest equation of the
    % structural form of the SVAR, and it is givem by the following
    % objects:
    % On the left hand-side: inva0(1,:)
    % On the right hand-side: structural_mu(1,:) and structural_bb(1,:)
    %
    % Then, we need to (say) increase the coefficient on inflation in the
    % structural monetary rule, that is, in each AR matrix, the coefficient
    % on the third variable (inflation) in the first row ...
    %
    % To see how the structural monetery ruke changes, let's keep track of
    % the old values, and then let's compare them with the new values:
    % original_structural_bb(1,:)=structural_bb(1,:);
    %
    % So we run a loop:
    for xx=1:LagOrder
        structural_bb(1,(xx-1)*N+3)=K*structural_bb(1,(xx-1)*N+3);
    end
    % To see that we have truly changed the coefficients ion inflation that
    % we wanted to change, let's check:
    %                  original_structural_bb(1,:)-structural_bb(1,:)
    % and you see that we got it right: we indeed changed the coefficinents
    % on inflation in the structural monetary policy rule in the SVAR ...
    % Now, starting from the equation:
    %       inva0 * Y(t) = structural_mu + structural_bb * [lags of Y(t)] + a0 * StructuralShocks(t)
    % we multiply everything by a0, and we get the counterfactual mu and
    % bb:
    counterfactual_mu=a0*structural_mu;
    counterfactual_bb=a0*structural_bb;
    %
    % Important thing here: when you change the parameters of the monetary
    % policy rule, you change the structure of the VAR ...
    % Therefore, even if the VAR was stationary originally (because we
    % imposed that in estimation) we do not have any assurance that it is
    % still going to be stationary after we have changd the monetary policy
    % rule ...
    % So we need to check this, and we only do the counterfactual for those
    % draws for whcuih this is the case ...
    R=varroots(LagOrder,N,[mu counterfactual_bb]);
    MaxAbsRoots=max(abs(R));
    %
    if MaxAbsRoots<1
        CounterfactualX=X';
        %
        tt=LagOrder+1;
        while tt<=T
            Shocks=StructuralShocks(:,tt-LagOrder,ss);
            % Here we re-run history :
            CounterfactualX(:,tt)=counterfactual_mu+counterfactual_bb*vec(fliplr(CounterfactualX(:,tt-LagOrder:tt-1)))+a0*Shocks;
            tt=tt+1;
        end
        % Here we store it:
        CounterfactualXX(:,:,ss)=CounterfactualX';
        CounterfactualXX(:,3,ss)=RollingMean([ones(11,1)*NaN; CounterfactualXX(:,3,ss)],12);
        Index(ss)=1;
    end
    ss=ss+1;
end
% Let's start by looking at how many successful draws we had:
sum(Index)/NN
% 89.8 per cent is good: If we had considered a more drastic change in the monetary
% rule, we might have ended up with less draws ...
%
SuccessfulDraws=find(Index);
CounterfactualXX=CounterfactualXX(:,:,SuccessfulDraws);
% Getting the percentiles of the counterfactual series:
CounterfactualSeries=sort(CounterfactualXX,3);
Percentiles=fix(size(CounterfactualSeries,3)*[0.5 0.16 0.84 0.05 0.95]');
PercentilesCounterfactualSeries=CounterfactualSeries(:,:,Percentiles);
%
%
xx=1;
while xx<=N
    %
    Counterfactual=squeeze(PercentilesCounterfactualSeries(:,xx,:));
    %
    figure(1)
    subplot(1,N,xx)
    plot(Time,Counterfactual(:,1),'k',Time,Counterfactual(:,4:5),'r',Time,XX(:,xx),'b',Time,zeros(size(Time)),'b:')
    xlim([Time(1) Time(T)])
    if xx==1
        title('Federal funds rate')
    elseif xx==2
        title('Spread')
    elseif xx==3
        title('Inflation')
    elseif xx==4
        title('Unemployment rate')
    else
    end
    %
    xx=xx+1;
end
% What you see is pretty typical of policy counterfactuals: not much
% happens ...



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







