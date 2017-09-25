function [A, B, B_boot, B_boot_Kilian, mean_B_boot_check, shock_vec, A99, A1, A95, A5, A86, A16] = ...
      cholboot(dataset,lag_number,which_shock,total_extractions)

if lag_number == 0
      error('Number of lags cannot be zero')
end

%Taking the lag on the dataset
for ilag = 1:lag_number+1
      eval(['data_', num2str(ilag), ' = dataset(2+lag_number-ilag:end+1-ilag,:);'])
end

%Building the matrix of regressor (which are the lags of the left hand side)
reg = zeros(length(dataset)-lag_number,lag_number*size(dataset,2));
for ilag = 2:lag_number+1
      for ivar = 1:size(dataset,2)
            eval(['reg(:,ivar+size(dataset,2)*(ilag-2)) = data_',num2str(ilag),'(:,ivar);'])
      end
end

%OLS in the Reduced Form Vector Autoregression
B = (reg'*reg)^(-1)*(reg'*data_1);

%Evaluate the residuals (res = yhat - y)
res = data_1 - reg*B;

%Compute the variance-covariance matrix
sigma = res'*res;

%Static rotation matrix
A = chol(sigma,'lower');

% Proving that A is what we are looking for A'u'uA = res*res where u'u = I
zero1 = A*A' - sigma;
zero2 = sum(sum(zero1.^2));
if zero2 > 10^(-16)
      error('The rotation matrix is not correct. Check the code')
end

%Shock we want to focus on
shock_vec = zeros(1,size(dataset,2));
shock_vec(1,which_shock) = 1;

reg_boot = zeros(length(dataset)-2*lag_number,lag_number*size(dataset,2),total_extractions);
for i_repeat = 1:total_extractions
      %Random extraction of the residuals
<<<<<<< HEAD
      res_boot = res(randsample(size(res,1),1),:);      
      %Building many bootstrapped datasets
      reg_new = reg(1,:);
      for i_boot = 1:size(res,1)
      res_boot = res(randperm(size(res,1),1),:)
      yhat = reg_new*B;
      ystar = yhat + res_boot
      reg_new = [ystar reg(i_boot,size(dataset,2)+1:end)]
      dataset_boot(i_boot,:,i_repeat) = ystar;
      end
=======
      rand_sorter = randsample(size(res,1),size(res,1));
%       rand_sorter = randsample(size(res,1),size(res,1), true); % true = with replacement

      for i_s = 1:size(rand_sorter,1)
            res_boot(i_s,:,i_repeat) = res(rand_sorter(i_s),:); % <---
      end
      %Building many bootstrapped datasets
%       dataset_boot(:,:,i_repeat) = reg*B + res_boot(:,:,i_repeat);
      dataset_boot(:,:,i_repeat) = reg*B + res_boot(:,:,i_repeat)*A;
>>>>>>> b4e741ad45a0017dc35b52ee9d08497a4bc9748e
      for ilag = 1:lag_number+1
            eval(['data_boot', num2str(ilag),'(:,:,i_repeat) = dataset_boot(2+lag_number-ilag:end+1-ilag,:,i_repeat);'])
      end
      for ilag = 2:lag_number + 1
            for ivar = 1:size(dataset_boot,2)
                  eval(['reg_boot(:,ivar+size(dataset_boot,2)*(ilag-2),i_repeat) = data_boot',num2str(ilag),'(:,ivar,i_repeat);'])
            end
      end
      %OLS in the Reduced Form Vector Autoregression for each boot_dataset
      B_boot(:,:,i_repeat) = (reg_boot(:,:,i_repeat)'*reg_boot(:,:,i_repeat))^(-1)*...
            ((reg_boot(:,:,i_repeat))'*data_boot1(:,:,i_repeat));
      
      %Kilian Correction
      adjmeanB = mean(B_boot,3);
<<<<<<< HEAD
      B_boot(:,:,i_repeat) = 2*B_boot(:,:,i_repeat) - adjmeanB;

=======
      B_boot_Kilian(:,:,i_repeat) = 2*B_boot(:,:,i_repeat) - adjmeanB;
end
mean_B_boot_check = mean(B_boot,3);
for i_repeat = 1:total_extractions
>>>>>>> b4e741ad45a0017dc35b52ee9d08497a4bc9748e
      %Evaluate the residuals (res = yhat - y)
      res_boottwo(:,:,i_repeat) = ...
            data_boot1(:,:,i_repeat) - reg_boot(:,:,i_repeat)*B_boot(:,:,i_repeat);
      %Compute the variance-covariance matrix
      sigma_boot(:,:,i_repeat) = res_boottwo(:,:,i_repeat)'*res_boottwo(:,:,i_repeat);
      %Static rotation matrix
      A_boot(:,:,i_repeat) = chol(sigma_boot(:,:,i_repeat),'lower');
end
% Proving that A_boot is what we are looking for A'u'uA = res*res where u'u = I
zero1_boot(:,:,i_repeat) = A_boot(:,:,i_repeat)*A_boot(:,:,i_repeat)' - sigma_boot(:,:,i_repeat);
zero2_boot = sum(sum(sum(zero1_boot.^2)));
if zero2_boot > 10^(-10)
      error('The rotation matrix of the bootstrap is not correct. Check the code')
end

%Confidence Intervals
sort_A_boot = sort(A_boot,3);
A99 = sort_A_boot(:,:,ceil(total_extractions*0.99));
A1 = sort_A_boot(:,:,ceil(total_extractions*0.01));
A95 = sort_A_boot(:,:,ceil(total_extractions*0.95));
A5 = sort_A_boot(:,:,ceil(total_extractions*0.05));
A86 = sort_A_boot(:,:,ceil(total_extractions*0.86));
A16 = sort_A_boot(:,:,ceil(total_extractions*0.16));

B(1,:)
mean_B_boot_check = mean(B_boot,3);
mean_B_boot_check(1,:)
res_boot(1,:,1)
res

end













