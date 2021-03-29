function set = no_update_2( k_run )

    % number of iterations
    set.nz = 10;
    
    % first cycle of the model update
    set.gUpdateCycle = 2;
    
    % learning rate
    set.gamma = 0.7;
    
    set.enable_stalibity_stop = true;
    
    set.enable_cycle_error = true;
    
    % enable fixed or random reference
    set.ref_type = 'single';
    
    % enable fixed or random model and system
    set.system_type = 'single';
    
    set.enable_var_filter_time = true;
    
    set.enable_fixed_time = true;
    set.tend = 5;
    set.nt = 200;    
        
    % enble plots
    set.enable_plots = false;
    
    % enable model update
    if k_run == 1
        set.enable_error_update = true;
    else
        set.enable_error_update = false;  
    end
    
    % type of model update
    set.model_update_type = 'fmincon_gm';
    
    % type of control
    set.control_type = 'quadprog';
    
    % noise settings
    set.enable_noise_y = false;
    set.enable_noise_u = false;
    set.mu_u = 0.1;
    set.sigma_u = 0.01;
    set.mu_y = 0.5;
    set.sigma_y = 0.01;
    
    
    
end