function [A_boot, B_boot, fake_impact_boot] = bootstrappedCIs(B, nburn, res, nsimul, which_correction, blocksize, nvar, ...
        nlags, pos_rel_prices, which_shocks, which_variable,H)
% This is a first pass bootstrap that's not quite correct because it needs
% to use Ryan's ID strategy every time to recover "bootstrapped gammas" and
% generate IRFs from A*gamma.

% Generate bootstrapped data samples
dataset_boot = data_boot(B, nburn, res, nsimul, which_correction, blocksize);

% Redo VAR nsimul times on the bootstrapped datasets
impact_boot    = zeros(nvar,length(which_shocks),nsimul);
B_boot         = zeros(nvar*nlags+1,nvar,nsimul);
% for i_simul = 1:nsimul
%     [A_boot(:,:,i_simul), B_boot(:,:,i_simul), ~, ~] = ...
%         sr_var(dataset_boot(:,:,i_simul), nlags);
% end

disp('Going into the 1st bootstrap loop...')
for i_simul = 1:nsimul
    [A_step1, B_boot(:,:,i_simul), ~, ~] = sr_var(dataset_boot(:,:,i_simul), nlags);
%     [impact_boot(:,:,i_simul),~,~,~,~,~] =  ...
%         Ryan_two_stepsID(which_variable,which_shocks,H,B_boot(:,:,i_simul),A_step1, pos_rel_prices);
end


% Kilian correction
[B_corrected,  bias]    = kilian_corretion(B, B_boot);
dataset_boot_corrected  = data_boot(B_corrected, nburn, res, nsimul, which_correction, blocksize);
impact_boot_corrected   = zeros(nvar,length(which_shocks),nsimul);
B_boot_corrected        = zeros(nvar*nlags+1,nvar,nsimul);
% for i_simul = 1:nsimul
%     [A_boot_corrected(:,:,i_simul), B_boot_corrected(:,:,i_simul), ~, ~] = ...
%         sr_var(dataset_boot_corrected(:,:,i_simul), nlags);
% end

disp('Going into the 2nd bootstrap loop...')
for i_simul = 1:nsimul
    disp(['Simulation ', num2str(i_simul), ' out of ', num2str(nsimul)])
    [A_corrected, B_boot_corrected(:,:,i_simul), ~, ~] = sr_var(dataset_boot_corrected(:,:,i_simul), nlags);
    [impact_boot_corrected(:,:,i_simul),~,~,~,~,~] = Ryan_two_stepsID(which_variable,which_shocks,H, ...
        B_boot_corrected(:,:,i_simul),A_corrected, pos_rel_prices);
end

disp('Finishing the bootstrap...')


B_boot_test = mean(B_boot_corrected,3); %It should be very close to B
bias_test = sum(sum(abs(B - B_boot_test)));
if bias < bias_test
    warning('Kilian correction should decrease the bias of beta and mean(beta_boot).')
end

A_boot = impact_boot_corrected;
B_boot = B_boot_corrected;

fake_impact_boot = zeros(nvar,nvar,nsimul);
for i_simul = 1:nsimul
    fake_impact_boot(:,which_shocks,i_simul) = impact_boot_corrected(:,:,i_simul);
end



