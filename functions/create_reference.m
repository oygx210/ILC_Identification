function r = create_reference(time, ref_type, enable_var_filter_time)


switch ref_type
        
    case 'rand'
        rand_traj = rand;
        
        Range = [-1,1];
        Band = [0 0.05];
        prbs = idinput(2540,'prbs',Band,Range);
        short_prbs = prbs( (1:time.nt) + floor( (length(prbs) - time.nt)*rand_traj ));
        
        
        if enable_var_filter_time
            Tf = 0.1*rand + 0.05;
        else
            Tf = 0.05;
        end
        filter = tf(1,[Tf, 1])^7;
        
        r = lsim(filter, short_prbs, time.vec);
        
    case 'single'
        
        rand_traj = 0.1;
        
        Range = [-1,1];
        Band = [0 0.05];
        prbs = idinput(2540,'prbs',Band,Range);
        short_prbs = prbs( (1:time.nt) + floor( (length(prbs) - time.nt)*rand_traj ));
        
        Tf = 0.05;
        filter = tf(1,[Tf, 1])^7;
        
        r = lsim(filter, short_prbs, time.vec);
        
end
end