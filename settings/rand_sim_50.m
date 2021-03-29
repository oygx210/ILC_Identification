function set = rand_sim_50( ~ )

    % number of iterations
    set.nz = 10;
    
    % first cycle of the model update
    set.gUpdateCycle = 2;
    
    % learning rate
    set.gamma = 0.7;
    
    set.enable_stalibity_stop = true;
    
    set.enable_cycle_error = true;
    
    % enable fixed or random reference
    set.ref_type = 'rand';
    
    % enable fixed or random model and system
    set.system_type = 'rand';
    
    set.enable_var_filter_time = false;
    
    % obsolete,fixed time was only required for plotting
    set.enable_fixed_time = false;
    set.tend = 2;
    set.nt = 200;    
        
    % enble plots
    set.enable_plots = false;
    
    % enable model update
    set.enable_error_update = true;

    % type of model update
    set.model_update_type = 'algebraic_gme';
    
    % type of control
    set.control_type = 'quadprog';
    
    % eliminate index y
    % noise settings
    set.enable_noise_y = false;
    set.mu_y = 0.5  *2*(rand-0.5);
    set.sigma_y = 0.01  *2*(rand-0.5);
    
    set.enable_model_chec = true; 
    set.res_threshold = 0.01;
    
end