function set = load_default_settings()

    % number of iterations
    set.nz = 5;
    
    % first cycle of the model update
    set.gUpdateCycle = 2;
    
    % learning rate
    set.gamma = 0.7;
    
    % enable fixed or random reference
    set.ref_type = 'transport';
    
    % enable fixed or random model and system
    set.system_type = 'transport';
    
    % enble plots
    set.enable_plots = false;
    
    % enable model update
    set.enable_error_update = true;
    
    % type of model update
    set.model_update_type = 'tfest_gm';
    
    % type of control
    set.control_type = 'quadprog';
    
    % noise settings
    set.enable_noise_y = false;
    set.enable_noise_u = false;
    set.mu_u = 0.1;
    set.sigma_u = 0.01;
    set.mu_y = 0.5;
    set.sigma_y = 0.01;
    
    set.enable_plot_vel_error = true;
    set.enable_cycle_error = false;
    
end