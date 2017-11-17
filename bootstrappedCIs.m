function [fake_A_boot, B_boot] = bootstrappedCIs(B, nburn, res, nsimul, which_correction, ...
        blocksize, nvar, nlags, pos_rel_prices, which_shocks, which_variable,H,which_ID,impact)
    % This is a first pass bootstrap that's not quite correct because it needs
    % to use Ryan's ID strategy every time to recover "bootstrapped gammas" and
    % generate IRFs from A*gamma.
    % which_ID chooses between the following identification strategies
    % 'barskysims', 'FEVmax_sr', 'Ryan_two_stepsID', and 'sr'
    % impact is gam from the point estimation. It is a (nvar,2)
    
    % Generate bootstrapped data samples
    dataset_boot = data_boot(B, nburn, res, nsimul, which_correction, blocksize);
    
    % Redo VAR nsimul times on the bootstrapped datasets
    impact_boot    = zeros(nvar,length(which_shocks),nsimul);
    B_boot         = zeros(nvar*nlags+1,nvar,nsimul);
    
    disp('Going into the bootstrap loop...')
    delt = 10^(-4); %parameter to decrease B_boot to make it stationary
    for i_simul = 1:nsimul
        [~, B_boot(:,:,i_simul), ~, ~] ...
            = sr_var(dataset_boot(:,:,i_simul), nlags);
        % don't need to check this b/c we know our VAR point estimate is
        % stationary
%         %Checking if the VAR is stationary
%         flag = test_stationarity(B_boot(:,:,i_simul)')
%         while flag == 1
%             B_boot(:,:,i_simul) ...
%                 = (1 - delt)*B_boot(:,:,i_simul);
%             flag = test_stationarity(B_boot(:,:,i_simul)');
%             disp('I am doing Vitos stationarization of B_boot')
%         end
    end
    
    % Kilian correction
    [B_corrected,  bias]    = kilian_corretion(B, B_boot);
    dataset_boot_corrected  = data_boot(B_corrected, nburn, res, nsimul, which_correction, blocksize);
    impact_boot_corrected   = zeros(nvar,length(which_shocks),nsimul);
    B_boot_corrected        = zeros(nvar*nlags+1,nvar,nsimul);
    
    switch which_ID
        
        case 'sr'
            for i_simul = 1:nsimul
                disp(['Simulation ', num2str(i_simul), ' out of ', num2str(nsimul)])
                [impact_boot_corrected(:,:,i_simul), B_boot_corrected(:,:,i_simul), ~, ~] = ...
                    sr_var(dataset_boot_corrected(:,:,i_simul), nlags);
            end
            
        case 'Ryan_two_stepsID'
            for i_simul = 1:nsimul
                disp(['Simulation ', num2str(i_simul), ' out of ', num2str(nsimul)])
                [A_corrected, B_boot_corrected(:,:,i_simul), ~, ~] = ...
                    sr_var(dataset_boot_corrected(:,:,i_simul), nlags);
                [impact_boot_corrected(:,:,i_simul),~,~,~,~,~] = Ryan_two_stepsID(which_variable,which_shocks,H, ...
                    B_boot_corrected(:,:,i_simul),A_corrected, pos_rel_prices);
               %
                
                
            end
            
        case 'FEVmax_sr'
            pos_mich_index = 3;
            disp(['For the restriction, the position of Mich Index should be ' num2str(pos_mich_index)])
            for i_simul = 1:nsimul
                disp(['Simulation ', num2str(i_simul), ' out of ', num2str(nsimul)])
                [A_corrected, B_boot_corrected(:,:,i_simul), ~, ~] = ...
                    sr_var(dataset_boot_corrected(:,:,i_simul), nlags);
                [impact_boot_corrected(:,:,i_simul),~,~,~,~,~] = FEVmax_sr_ID(which_variable,which_shocks,H, ...
                    B_boot_corrected(:,:,i_simul),A_corrected, pos_mich_index);
            end
            
        case 'barskysims'
            for i_simul = 1:nsimul
                disp(['Simulation ', num2str(i_simul), ' out of ', num2str(nsimul)])
                [A_corrected, B_boot_corrected(:,:,i_simul), ~, ~] = ...
                    sr_var(dataset_boot_corrected(:,:,i_simul), nlags);
                [impact_boot_corrected(:,:,i_simul),~,~,~,~,~] = ...
                    barskysims(which_variable,which_shock,H,B_boot_corrected(:,:,i_simul),A_corrected);
            end
            
    end
    
    disp('Finishing the bootstrap...')
    
    B_boot_test = nanmean(B_boot_corrected,3); %It should be very close to B
    bias_test = sum(sum(abs(B - B_boot_test)));
    if bias < bias_test
        warning('Kilian correction should decrease the bias of beta and mean(beta_boot).')
    end
    
    %If the ID is max something we obtain just gamma which is not the full
    %impact and then we create a fake squared matrix
    if strcmp(which_ID, 'sr')
        fake_A_boot = impact_boot_corrected;
    else
        A_boot = impact_boot_corrected;
        B_boot = B_boot_corrected;
        
        fake_A_boot = zeros(nvar,nvar,nsimul);
        for i_simul = 1:nsimul
            fake_A_boot(:,which_shocks,i_simul) = A_boot(:,:,i_simul);
        end
    end
    
end



