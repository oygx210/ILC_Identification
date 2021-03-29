function set = noise_sim_50( k1 )

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
    set.enable_error_update = true;

    
    % type of model update
    if k1 <= 25
        set.model_update_type = 'tfest_gm_io';
    else
        set.model_update_type = 'tfest_gm';
    end
    
    % type of control
    set.control_type = 'quadprog';
    
    % noise settings
    set.enable_noise_y = true;
    set.mu_y = 0.5  *2*(rand-0.5);
    set.sigma_y = 0.01  *2*(rand-0.5);
    
end