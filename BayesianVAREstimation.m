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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumberOfDraws=1000; % The number of draws from the posterior distribution
Intercept='Y'; % 'Y' if you want an intercept, 'N' otherwise
Stationary='Y'; % 'Y' if you want to impose stationarity, 'N' otherwise
%
disp('-------------------------') 
disp('Important point to notice') 
disp('-------------------------') 
%
% Since we are here working within a Bayesian context, we are in principle 
% allowed to impose any priors that we want on the economy, as long as we 
% can 'reasonably justify them' based on some meaningful argument ...
%
% This means that, if we want, we are perfectly entitled to impose the
% prior belief that the economy -- and therefore the VAR -- be stationary.
%
Trend=0; % This is the order of the deterministic time trend to be included in the VAR.
%          0 = no time trend; 1 = linear time trend; 2 = quadratic time trend
%
disp('----------------------------------------------------------------------') 
disp('Then, these are the priors (see discussion on slide 42 of October 18):') 
disp('----------------------------------------------------------------------') 
%
% B0,N0,S0,v0 = the Bayesian priors
%
disp('--------------------------------------------------------------------------------------------') 
disp('If they are not inputted in the code, they are set to the values found in Uhlig (1998, 2005)') 
disp('--------------------------------------------------------------------------------------------') 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LagOrder=12;
%
[B,Sigma,U,TimeTrend]=BayesianVARp(X,LagOrder,Intercept,Trend,NumberOfDraws,Stationary);
disp('--------------------------------------------------------------------------------------------') 
disp('After running the code, show how it works') 
disp('--------------------------------------------------------------------------------------------') 
%
% Here we are extracting the relevant objects from the output produced by
% the code BayesianVARp ...
%
if Trend>0
    CoefficientsOnTimeTrend=B(:,size(B,2),:);
    B=B(:,1:size(B,2)-1,:);
end
if Intercept=='Y'
    MU=B(:,1,:);
    B=B(:,2:size(B,2),:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('---------------------------------------------------') 
disp('Now we plot the results. Lets do it for example for') 
disp('the elements of the VAR covariance matrix Sigma ...') 
disp('---------------------------------------------------') 
%
figure(1)
for hh=1:N
    for kk=1:N
        PosteriorDistribution=squeeze(Sigma(hh,kk,:)); % Turning the 3-dimensional object Sigma(hh,kk,:) into a vector via the command squeeze
        PosteriorDistribution=sort(PosteriorDistribution); % Sorting the elements of the vector in ascending order
        LimInf=PosteriorDistribution(fix(0.01*NumberOfDraws));
        LimSup=PosteriorDistribution(fix(0.99*NumberOfDraws));
        subplot(N,N,(hh-1)*N+kk)
        hist(PosteriorDistribution,100,'k')
        xlim([LimInf LimSup])
    end
end
disp('--------------------------------------------------------')
disp('So you see that the entire logic is identical to what we') 
disp('saw in the Classical case, when we estimated the VAR via') 
disp('OLS, we bootstrapped it, and then we plotted the histogram') 
disp('of the bootstrapped distributions of the various elements')
disp('of the VAR ...') 
disp('--------------------------------------------------------')
disp('This means that everything we did then based on the reduced-form ') 
disp('VAR estimated via OLS -- computing projections, computing confidence') 
disp('intervals for the projections, etc. etc. -- can be done here too.') 
disp('--------------------------------------------------------')
disp('So lets compute a projection for the U.S. unemployment rate') 
disp('So lets compute a projection for the U.S. unemployment rate') 
disp('5 years ahead: the data we have in the spreadsheet were ') 
disp('downloaded 2 hours ago, so these are the data that policymakers') 
disp('and forecasters currently have to project the evolution of U.S.') 
disp('unemployment in the future ...') 
disp('Clearly, our VAR is extraordinarly simple, but it is nonetheless') 
disp('instructive to see what comes out of this forecasting exercise ...') 
disp('---------------------------------------------------') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
ForecastHorizon=5*12;
%
ForecastedX=zeros(1+ForecastHorizon,N,NumberOfDraws);
TimeForecast=(Time(T)):0.25:(Time(T)+0.25*ForecastHorizon);
%
xx=1;
while xx<=NumberOfDraws
    VARCovarianceMatrix=Sigma(:,:,xx);
    SqrtV=sqrtm(VARCovarianceMatrix);
    % SqrtV*SqrtV'-VARCovarianceMatrix is of an order of magnitude of 10 to
    % the minus 14, so it's fine ...
    % We use SqrtV in a moment to rescale the disturbances that we
    % randomly draw:
    Shocks=SqrtV*randn(N,LagOrder+1+ForecastHorizon);
    %
    forecastedx(:,1:LagOrder)=X(T-LagOrder+1:T,:)';
    for tt=LagOrder+1:LagOrder+1+ForecastHorizon
        % Here we stochastically simulate the VAR into the future:
        forecastedx(:,tt)=MU(:,:,xx)+B(:,:,xx)*vec(fliplr(forecastedx(:,tt-LagOrder:tt-1)))+Shocks(:,tt);
    end
    ForecastedX(:,:,xx)=forecastedx(:,LagOrder+1:LagOrder+1+ForecastHorizon)';
    xx=xx+1;
end
%
% Here we extract the projection for the unemployment rate:
% Remember the order of the variables: 
% X=[FEDFUNDS Spread CPIInflation UnemploymentRate];
ForecastedUnemploymentRate=squeeze(ForecastedX(:,4,:));
% Here we sort -- that is: we order -- the projections:
ForecastedUnemploymentRate=sort(ForecastedUnemploymentRate,2);
% Here we get the percentiles:
Percentiles=NumberOfDraws*[0.5 0.16 0.84 0.05 0.95 0.01 0.99]';
ForecastedUnemploymentRate=ForecastedUnemploymentRate(:,Percentiles);
%
% Then we plot:
%
Zero=zeros(size(Time));
ZeroForecast=zeros(size(TimeForecast));
%
figure(2)
subplot(1,2,1)
plot(Time,X(:,4),'k',TimeForecast,ForecastedUnemploymentRate(:,1),'b',TimeForecast,ForecastedUnemploymentRate(:,2:3),'r')
xlim([Time(1) TimeForecast(length(TimeForecast))])
title('Unemployment rate')

























