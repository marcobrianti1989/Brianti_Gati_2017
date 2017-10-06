function [beta_corrected] = kilian_corretion(beta,beta_boot)
% beta is the matrix of regressor coefficients from the VAR. It is (nlga*nvar x nvar)
% beta_boot is the matrix of regressor coefficients from the first round bootstrapped VARs 

average_beta_boot    = mean(beta_boot,3); %Expectation over the number of simulations
beta_corrected       = 2*beta - average_beta_boot; %Remove the bias from beta

end


