function out = master_sim_v3(varargin)

% Jonas Weigand <jonas.weigand@mv.uni-kl.de>
% Gajanan Kanagalingam <gajanan.kanagalingam@gmail.com>
% 08-02-2021





% load settings
if isempty(varargin)
    
    close all
    clc
    
    % load default settings
    set = load_default_settings();
    
    % and override
    set.enable_plots = true;
    
else
    set = varargin{1};
    
end


% load system, model and time
switch set.system_type
    
    case 'rand'
        [sys, time] = create_rand_system();
        
    case 'transport'
        [sys, time] = create_transport_system();
        
    case 'single'
        [sys, time] = create_single_system();
        
    otherwise
        disp('No system is chosen.')
        out = [];
        return
        
end

% create reference
r = create_reference(time, set.ref_type);



%% initialise variables

cycle_error_perfect     = zeros(set.nz, 1);
cycle_error_each_iter   = zeros(set.nz, 1); % all cases plus perfect measure
yIterAll                = zeros(set.nz, time.nt); % all cycles, all timesteps, with and without model update

g_m_all         = zeros(time.nt, 3); % initial, real, update
g_m_all(:, 1)   = sys.G_M(:,1);
g_m_all(:, 2)   = sys.G_S(:,1);


% the initial input equals the reference trajectory
u = r;

% for all cycles
for k2 = 1:set.nz
    
    % add noise to input
    if set.enable_noise_u
        u = u + set.mu_u + set.sigma_u*randn(size(u));
        
    end
    
    % calculate the "real" system
    y_LTI = sys.G_S * u;
    
    % add noise to output
    if set.enable_noise_y
        
        y_LTI = y_LTI + set.mu_y + set.sigma_y*randn(size(y_LTI));
        
    end
    
    % save norm of error
    cycle_error_each_iter(k2) = norm(r - y_LTI);
    
    % save y
    yIterAll(k2, :) = y_LTI;
    
    % save perfect error for comparision
    if k2 <= set.gUpdateCycle
        cycle_error_perfect(k2,3) = norm(r - y_LTI);
    else
        cycle_error_perfect(k2,3) = (1 - set.gamma)*cycle_error_perfect(k2-1,3);
    end
    
    
    % model update
    if set.enable_error_update && k2 == set.gUpdateCycle
        
        % get the current and previous error
        e_current   = r - yIterAll(k2,:)';
        e_prev      = r - yIterAll(k2-1,:)';
        
        % save data
        data.e_current = e_current;
        data.e_prev = e_prev;
        data.y_LTI = y_LTI;
        data.r = r;
        data.u = u;
        data.g_me_init = g_m_all(:, 1) .* e_prev;
        data.g_m_init = g_m_all(:, 1);
        data.gamma = set.gamma;
        
        % initialise g_m and g_me
        g_m_update = [];
        g_me_update = [];
        
        switch set.model_update_type
            
            case 'tfest_gme'
                
                g_me_update = update_tfest_gme(sys, time, data);
                
            case 'fmincon_gme'
                
                g_me_update = update_fmincon_gme(sys, time, data);
                
            case 'algebraic_gme'
                
                g_me_update = update_algebraic_gme(sys, data, time);
                
            case 'fmincon_gm_io'
                
                g_m_update = update_fmincon_gm_io(data, time);
                
            case 'fmincon_gm'
                
                g_m_update = update_fmincon_gm(sys, data, time);
                
            case 'tfest_gm_io'
                
                g_m_update = update_tfest_gm_io(data, time);
                
            case 'tfest_gm'
                
                g_m_update = update_tfest_gm(sys, data, time);
                
        end
        
        % perform update, depending on g_me or g_m type
        if not(isempty(g_m_update))
            
            % transfer update to a matrix
            sys.G_M = convert_from_vector_to_matrix( g_m_update, time );
            
            g_m_all(:, 3) = sys.G_M(:,1);
            
        elseif not(isempty(g_me_update))
            
            % transfer update to a matrix
            G_ME = convert_from_vector_to_matrix( g_me_update, time );
            
            % perform the update
            sys.G_M =  1/set.gamma*(sys.G_M - G_ME);
            
            g_m_all(:, 3) = sys.G_M(:,1);
            
        else
            disp('No model update applied.')
        end
    end
    
    % for each but the last cycle apply the update of the input trajectory
    if k2 < set.nz
        
        
        switch set.control_type
            
            % calulate the update using a quadratic form
            case 'quadprog'
                
                options = optimset('Display','off','Algorithm','interior-point-convex');
                
                % weighting matrix for quadratic optimization
                Q = eye(time.nt,time.nt);
                
                f = -2*(r - y_LTI)'*Q*sys.G_M;
                H = 2*sys.G_M'*Q*sys.G_M;

                
                
                if ~issymmetric(H)
                    H = (H' + H)/2;
                end
                
                u_update = quadprog(H,f, [],[],[],[],[],[],[], options);
                
                % or calculate the update using a nonlinear optimisation
            case 'fmincon'
                
                options = optimoptions('fmincon','Display',...
                    'off','Algorithm','sqp', ...
                    'MaxFunctionEvaluations', 1e5);
                
                costFun = @(u_use) cost_control( u_use , sys, time, r, y_LTI);
                
                u_update = fmincon( costFun, u , [],[],[],[],[],[],[],  options);
                
            case 'algebraic'
                
                
                W = eye(time.nt);
                I = eye(time.nt);
                
                tik_cand =10.^(-20:20);
                e_tik = inf(1,length(tik_cand));
                
                
                % calculate Tikhonov error
                for k3=1:1:length(tik_cand)
                    
                    M = (sys.G_M'*W*sys.G_M + tik_cand(k3)*I);
                    if rcond(M) > 10^-15
                        
                        % etk = norm(  r - Gm*u + (Gm'*Gm + I)^-1*Gm'*e )
                        e_tik(k3) = norm(r - sys.G_M*(u + M^-1*W*sys.G_M'*(r-y_LTI)));
                    end
                end
                
                [~, ind_tik]=min(e_tik);
                
                u_update=(sys.G_M'*W*sys.G_M+tik_cand(ind_tik)*I)^-1*W*sys.G_M'*(r-y_LTI);
                
                
            otherwise
                
                u_update = 0;
                disp('No control update applied.')
        end
        
        % apply the update to the trajectory
        u = u + set.gamma * u_update;
        
    end
end


%% Outputs
out.g_m_all = g_m_all;
out.time = time;
out.sys = sys;
out.yIterAll = yIterAll;
out.r = r;
out.set = set;
out.cycle_error_each_iter = cycle_error_each_iter;
out.cycle_error_perfect = cycle_error_perfect;



%% Visualisation
if set.enable_plots
    
    figure
    subplot(3,1,1)
    hold on
    plot(time.vec, yIterAll(1,:,1),'b');
    plot(time.vec, yIterAll(2,:,1),'m');
    plot(time.vec, yIterAll(end,:,1),'k');
    plot(time.vec, r,'r-.');
    legend('1st Cycle','2nd Cycle', 'Final Cycle', 'Reference', 'Location','eastoutside')
    xlabel('Time')
    ylabel('Output')
    subplot(3,1,2)
    hold on
    plot(time.vec, r - yIterAll(1,:)','b');
    plot(time.vec, r - yIterAll(2,:)','m');
    plot(time.vec, r - yIterAll(end,:)','k');
    legend('1st Error','2nd Error', 'Final Error', 'Location','eastoutside')
    xlabel('Time')
    ylabel('Error')
    subplot(3,1,3)
    hold on
    plot(time.vec, g_m_all(:,2),'m');
    plot(time.vec, g_m_all(:,1),'b');
    plot(time.vec, g_m_all(:,3),'k-.');
    legend('True', 'Initial', 'Updated', 'Location','eastoutside')
    xlabel('Time')
    ylabel('Impuse Response')
    
    figure
    grid on
    hold on
    plot(time.vec, g_m_all(:,2),'m');
    plot(time.vec, g_m_all(:,1),'b');
    plot(time.vec, g_m_all(:,3),'k-.');
    legend('True', 'Initial', 'Updated', 'Location','eastoutside')
    xlabel('Time')
    ylabel('Impuse Response')
    
    if set.enable_plot_vel_error
        figure
        subplot(2,1,1)
        hold on
        grid on
        pos1_err = r - yIterAll(1,:)';
        pos2_err = r - yIterAll(2,:)';
        pos3_err = r - yIterAll(end,:)';
        
            vel1 = [0; (pos1_err(2:end) - pos1_err(1:end-1))./time.dt];
            vel2 = [0; (pos2_err(2:end) - pos2_err(1:end-1))./time.dt];
            vel3 = [0; (pos3_err(2:end) - pos3_err(1:end-1))./time.dt];
            plot(time.vec, vel1,'b');
            % plot(time.vec, vel2,'m');
            plot(time.vec, vel3','k');
            legend('1st Vel Error', 'Final Vel Error', 'Location','eastoutside')
            xlabel('Time')
            ylabel('Vel Error')
                    subplot(2,1,2)
        hold on
        grid on
        pos1 = yIterAll(1,:)';
        pos2 = yIterAll(2,:)';
        pos3 = yIterAll(end,:)';
        
            vel1 = [0; (pos1(2:end) - pos1(1:end-1))./time.dt];
            vel2 = [0; (pos2(2:end) - pos2(1:end-1))./time.dt];
            vel3 = [0; (pos3(2:end) - pos3(1:end-1))./time.dt];
            vel_ref = [0; (r(2:end) - r(1:end-1))./time.dt];
            % plot(time.vec, vel1,'b');
            % plot(time.vec, vel2,'m');
            plot(time.vec, vel3','k');
            plot(time.vec, vel_ref','r');
            legend('1st Vel', 'Final Vel',  'Ref Vel', 'Location','eastoutside')
            xlabel('Time')
            ylabel('Vel')
    end

if set.enable_cycle_error

    figure
    hold on
    grid on
    plot(cycle_error_each_iter, 'r-*')
    plot(cycle_error_perfect, 'mo')
    ax = gca;
    ax.YScale = 'log';
    legend('Model Update', 'Model Init', 'Perfect', 'Location','southwest')
    xlabel('Cycle')
    ylabel('Error in L2-Norm')
end
end