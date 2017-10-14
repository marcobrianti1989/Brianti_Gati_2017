function [IRFsOnlyFirstShock,IRFsOnlySecondShock,IRFsOnlyThirdShock,IRFsOnlyFourthShock]=GetIRFsFourShocks(B,A0,N,LagOrder,Horizon)
%
YY=zeros(N,LagOrder+1+Horizon);
YY(:,LagOrder+1)=A0*[1 zeros(1,N-1)]';
for tt=LagOrder+2:LagOrder+1+Horizon
    YY(:,tt)=B*vec(fliplr(YY(:,tt-LagOrder:tt-1)));
end
YY=YY';
YY(:,1)=cumsum(YY(:,1));
YY(:,3)=cumsum(YY(:,3));
IRFsOnlyFirstShock=YY(LagOrder+1:LagOrder+1+Horizon,:);
%
YY=zeros(N,LagOrder+1+Horizon);
YY(:,LagOrder+1)=A0*[0 1 zeros(1,N-2)]';
for tt=LagOrder+2:LagOrder+1+Horizon
    YY(:,tt)=B*vec(fliplr(YY(:,tt-LagOrder:tt-1)));
end
YY=YY';
YY(:,1)=cumsum(YY(:,1));
YY(:,3)=cumsum(YY(:,3));
IRFsOnlySecondShock=YY(LagOrder+1:LagOrder+1+Horizon,:);
%
YY=zeros(N,LagOrder+1+Horizon);
YY(:,LagOrder+1)=A0*[zeros(1,2) 1 zeros(1,N-3)]';
for tt=LagOrder+2:LagOrder+1+Horizon
    YY(:,tt)=B*vec(fliplr(YY(:,tt-LagOrder:tt-1)));
end
YY=YY';
YY(:,1)=cumsum(YY(:,1));
YY(:,3)=cumsum(YY(:,3));
IRFsOnlyThirdShock=YY(LagOrder+1:LagOrder+1+Horizon,:);
%
YY=zeros(N,LagOrder+1+Horizon);
YY(:,LagOrder+1)=A0*[zeros(1,3) 1 zeros(1,N-4)]';
for tt=LagOrder+2:LagOrder+1+Horizon
    YY(:,tt)=B*vec(fliplr(YY(:,tt-LagOrder:tt-1)));
end
YY=YY';
YY(:,1)=cumsum(YY(:,1));
YY(:,3)=cumsum(YY(:,3));
IRFsOnlyFourthShock=YY(LagOrder+1:LagOrder+1+Horizon,:);




