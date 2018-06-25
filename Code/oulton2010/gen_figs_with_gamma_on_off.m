% gen_figs_with_gamma_on_off - This file just imports the results of
% main_prog_spillover_news.m with the optimal parameter values, but once
% gamma is its optimal value of 0.58, once it is its minimum value of 0.14
% (should be 0).
clear

% Load IRFs for when gamma is optimal
load IRFs_gamma_opt
IRFs_gam_opt = IRFs_VAR;
clear IRFs_VAR

% Load IRFs for when gamma is shut off (set to its minimum value)
load IRFs_gamma_min
IRFs_gam_min = IRFs_VAR;

print_figs = 'yes';
plot2IRFs(IRFs_gam_min,IRFs_gam_opt,T,which_shock,shocknames, varnames_matching, '\gamma = 0.14', '\gamma = 0.58', print_figs, base_path, 'Varying_gamma')