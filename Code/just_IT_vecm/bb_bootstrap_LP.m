function [Yboot Xboot] = bb_bootstrap_LP(Y,X,nsimul,lags)
% no bias correction yet - although it should not change substantially 

% Inputs:
% Y which is the projected variable at horizon h. Length is T
% X is the matrix of regressor. Length is T
% nsimul = number of simulations, i.e. how many bootstrapped datasets to generate
% As block_size we use approximately half of the dataset
% CI is confidence interval. Say 95% for example

T                  = length(Y); 
q                  = size(X,2); 
l                  = floor((T-q)^(1/3)); %quite arbitrary 

% need to do for every h
% there are (T-l+1) (1+q)-tuples

Yboot = nan(floor(T/l)*l,nsimul);
Xboot = nan(floor(T/l)*l,size(X,2),nsimul);
for isimul = 1:nsimul
    for j  = 1:floor(T/l)
      draw                           = randi(T-l+1,1); %check
      Yboot(1+(j-1)*l:j*l,isimul)    = Y(draw:draw+l-1);
      Xboot(1+(j-1)*l:j*l,:,isimul)  = X(draw:draw+l-1,:);
    end
end

end