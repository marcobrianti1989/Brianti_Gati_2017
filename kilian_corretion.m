function [beta_corrected, bias] = kilian_corretion(beta, beta_boot)
% beta is the matrix of regressor coefficients from the VAR. It is (nlga*nvar x nvar)
% beta_boot is the matrix of regressor coefficients from the first round bootstrapped VARs 

average_beta_boot    = mean(beta_boot,3); %Expectation over the number of simulations
beta_corrected       = 2*beta - average_beta_boot; %Remove the bias from beta
bias                 = sum(sum(abs(beta - average_beta_boot))); %useful to check with the second average_beta_boot

end


