function set = default()

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
set.model_update_type = 'fmincon_gm';

% type of control
set.control_type = 'quadprog';

% noise settings
set.enable_noise = false;
set.mu = 0.5;
set.sigma = 0.01;

% enable stability checks
set.enable_model_check = true;
set.res_threshold = 0.01;

end