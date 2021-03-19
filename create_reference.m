function r = create_reference(time, ref_type)


switch ref_type
    case 'transport_v0'
        
        safe = 1; % safety factor
        
        accel_s_max = 10; % m/s^2
        vel_s_max = 2.5; % m/s
        
        accel_c_max = 5; % m/s^2
        vel_c_max = 1; % m/s
        
        accel = zeros(1, time.nt);
        
        % change velocity to first straight
        
        t_accel_1 = safe * vel_s_max / accel_s_max;
        nt_1 = ceil(t_accel_1 / time.dt);
        accel(0.1*time.nt + (1:nt_1)) = accel_s_max;
        
        % change velocity to first curve
        t_accel_2 = safe * (vel_s_max - vel_c_max) / accel_s_max;
        nt_2 = ceil(t_accel_2 / time.dt);
        accel(0.5*time.nt + (1:nt_2)) = -accel_s_max;
        
        % change velocity to standstill
        t_accel_3 = safe * vel_c_max / accel_c_max;
        nt_3 = ceil(t_accel_3 / time.dt);
        accel(0.7*time.nt + (1:nt_3)) = -accel_c_max;
        
        s = tf('s');
        int = 1/s;
        
        % vel = lsim(int, accel, time.vec);
        pos = lsim(int^2, accel, time.vec);
        
        %
        %     figure
        %     subplot(3,1,1)
        %     plot(time.vec, accel)
        %     subplot(3,1,2)
        %     plot(time.vec, vel)
        %     subplot(3,1,3)
        %     plot(time.vec, pos)
        
        r = pos;
        
    case 'transport'
                
        safe = 1; % safety factor
        
        accel_s_max = 10; % m/s^2
        vel_s_max = safe* 2.5; % m/s
        
        accel_c_max = 5; % m/s^2
        vel_c_max = safe* 1; % m/s
        
        accel = zeros(1, time.nt);
        
        
        % dist_ref = vel_s_max * dt1 + vel_c_max * dt2;
        % t1 = a * t2;
        % dist_ref = vel_s_max * a * dt2 + vel_c_max * dt2;
        % dist_ref = (vel_s_max * a + vel_c_max) * dt2;
        
        a = 1.5;
        dist_ref = 10;
        ind_t0 = 50;
        dt2 = dist_ref / (vel_s_max * a + vel_c_max);
        dt1 = a * dt2;
        ind_dt1 = ceil(dt1 / time.dt);
        ind_dt2 = ceil(dt2 / time.dt);
        
        
        % change velocity to first straight
        
        t_accel_1 = vel_s_max / accel_s_max;
        nt_1 = ceil(t_accel_1 / time.dt);
        accel(ind_t0 + (1:nt_1)) = accel_s_max;
        
        % change velocity to first curve
        t_accel_2 = (vel_s_max - vel_c_max) / accel_s_max;
        nt_2 = ceil(t_accel_2 / time.dt);
        accel(ind_t0 + ind_dt1 + (1:nt_2)) = -accel_s_max;
        
        % change velocity to standstill
        t_accel_3 = vel_c_max / accel_c_max;
        nt_3 = ceil(t_accel_3 / time.dt);
        accel(ind_t0 + ind_dt1 + ind_dt2 + (1:nt_3)) = -accel_c_max;
        
        s = tf('s');
        int = 1/s;
        
        vel = lsim(int, accel, time.vec);
        pos = lsim(int^2, accel, time.vec);
        
        
            figure
            subplot(3,1,1)
            plot(time.vec, accel)
            subplot(3,1,2)
            plot(time.vec, vel)
            subplot(3,1,3)
            plot(time.vec, pos)
        
        dist = pos(end)
        
        r = pos;
                
        
    case 'rand'
        rand_traj = rand;
        
        Range = [-1,1];
        Band = [0 0.05];
        prbs = idinput(2540,'prbs',Band,Range);
        short_prbs = prbs( (1:time.nt) + floor( (length(prbs) - time.nt)*rand_traj ));
        Tf = 0.05;
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