function [A_boot, B_boot] = bootstrappedCIs(B, nburn, res, nsimul, which_correction, q, nvar, nlags)
% This is a first pass bootstrap that's not quite correct because it needs
% to use Ryan's ID strategy every time to recover "bootstrapped gammas" and
% generate IRFs from A*gamma.

% Generate bootstrapped data samples
dataset_boot = data_boot(B, nburn, res, nsimul, which_correction, q);

% Redo VAR nsimul times on the bootstrapped datasets
A_boot = zeros(nvar,nvar,nsimul);
B_boot = zeros(nvar*nlags+1,nvar,nsimul);
for i_simul = 1:nsimul
    [A_boot(:,:,i_simul), B_boot(:,:,i_simul), ~, ~] = ...
        sr_var(dataset_boot(:,:,i_simul), nlags);
end

% Kilian correction
[B_corrected,  bias] = kilian_corretion(B, B_boot);
dataset_boot_corrected = data_boot(B_corrected, nburn, res, nsimul, which_correction, q);
A_boot_corrected = zeros(nvar,nvar,nsimul);
B_boot_corrected = zeros(nvar*nlags+1,nvar,nsimul);
for i_simul = 1:nsimul
    [A_boot_corrected(:,:,i_simul), B_boot_corrected(:,:,i_simul), ~, ~] = ...
        sr_var(dataset_boot_corrected(:,:,i_simul), nlags);
end
B_boot_test = mean(B_boot_corrected,3); %It should be very close to B
bias_test = sum(sum(abs(B - B_boot_test)));
if bias < bias_test
    error('Kilian correction should decrease the bias of beta and mean(beta_boot).')
end

A_boot = A_boot_corrected;
B_boot = B_boot_corrected;