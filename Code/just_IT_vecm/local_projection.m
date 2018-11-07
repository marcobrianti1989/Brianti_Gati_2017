function [IRF,res,Rsquared,B_store,tuple_store,VarY] = local_projection(y,x,u,lags,H)

% H is the horizon of the irf
% y and x are (T,1); dependend and observable control variables
% u is exogenous regressor
% lags is how many lags of x and u we want to take in the regression
% regression: y(t+h) = alpha + B*u(t) + C1*u(t-1) + D1*y(t-1) + G1*x(t-1) + ...

for h = 1:H
      Y          = y(lags+h:end,:);
      X          = u(lags+1:end-h+1,:); %Plus 1 because when jj = 1 is contemporaneous
      if lags > 0 % which allows for controls
            for jj = 1:lags
                  if sum(sum(x.^2)) > 0 % Add into (big) X controls (small) x
                        X = [X, u(lags-jj+1:end-jj-h+1,:), y(lags-jj+1:end-jj-h+1,:), ...
                              x(lags-jj+1:end-jj-h+1,:)];
                  else % no controls 
                        X = [X, u(lags-jj+1:end-jj-h+1,:), y(lags-jj+1:end-jj-h+1,:)];
                  end
            end
      end
      % Add linear trend - [1:1:length(Y)]'
      trend              = [1:1:length(Y)]';      
      X                  = [X, ones(length(Y),1), trend];
      tuple_store{h}     = [Y X];    
      B                  = X'*X\(X'*Y); 
      IRF(h)             = B(1); 
      res{h}             = Y - X*B;
      Rsquared(h)        = 1 - var(res{h})/var(Y);
      B_store(:,h)       = B; 
      
      % Var(Y|all the shocks) - just eliminate the trend and the past
      Xcontrols          = X(:,2:end);
      Bcontrols          = Xcontrols'*Xcontrols\(Xcontrols'*Y); 
      VarY(h)            = var(Y - Xcontrols*Bcontrols);
      
end