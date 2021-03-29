function set = rand_sim_50( ~ )

% load default settings and only override changes
set = default();


% enable fixed or random reference
set.ref_type = 'rand';


% enable fixed or random model and system
set.system_type = 'rand';


% obsolete,fixed time was only required for plotting
set.enable_fixed_time = false;
set.tend = 2;
set.nt = 200;


% type of model update
set.model_update_type = 'algebraic_gme';


end