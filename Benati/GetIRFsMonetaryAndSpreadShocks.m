function [IRFsOnlyFirstShock,IRFsOnlySecondShock]=GetIRFsMonetaryAndSpreadShocks(B,A0,N,LagOrder,Horizon)
%
YY=zeros(N,LagOrder+1+Horizon);
YY(:,LagOrder+1)=A0*[1 zeros(1,N-1)]';
for tt=LagOrder+2:LagOrder+1+Horizon
    YY(:,tt)=B*vec(fliplr(YY(:,tt-LagOrder:tt-1)));
end
IRFsOnlyFirstShock=YY(:,LagOrder+1:LagOrder+1+Horizon)';
%
YY=zeros(N,LagOrder+1+Horizon);
YY(:,LagOrder+1)=A0*[0 1 zeros(1,N-2)]';
for tt=LagOrder+2:LagOrder+1+Horizon
    YY(:,tt)=B*vec(fliplr(YY(:,tt-LagOrder:tt-1)));
end
IRFsOnlySecondShock=YY(:,LagOrder+1:LagOrder+1+Horizon)';
