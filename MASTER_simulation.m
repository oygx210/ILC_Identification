function out = MASTER_simulation(varargin)

% Jonas Weigand <jonas.weigand@mv.uni-kl.de>
% Gajanan Kanagalingam <gajanan.kanagalingam@gmail.com>
% 08-02-2021


addpath('functions', 'settings')


% load settings
if isempty(varargin)
    
    close all
    clc
    
    % load default settings
    set = default();
    
    % enable plots for single simulations
    set.enable_plots = true;
    
else
    
    % if function is called by MASTER_multirun, import settings
    set = varargin{1};
    
end


% load system, model and time
switch set.system_type
    
    case 'rand'
        [sys, time] = create_rand_system( set.enable_fixed_time );
        
    case 'single'
        [sys, time] = create_single_system();
        
    otherwise
        warning('No system is chosen.')
        out = [];
        return
        
end

% create reference
r = create_reference(time, set.ref_type, set.enable_var_filter_time);



%% initialise variables

% save the perfect cycle error for each iteration
cycle_error_perfect     = zeros(set.nz, 1);

% save the cycle error for each iteration
cycle_error_each_iter   = zeros(set.nz, 1);

% save the output for all timesteps and all cycles
y_iter_all              = zeros(set.nz, time.nt);


% impulse responses for initial model, real model, and updated model
g_m_all         = zeros(time.nt, 3);
g_m_all(:, 1)   = sys.G_M(:,1); % copy initial model
g_m_all(:, 2)   = sys.G_S(:,1); % copy real model

% save previous model
G_M_prev    = sys.G_M;

% the initial input equals the reference trajectory
u = r;

% boolean for stopping the update
update_stopped = false;



%% for all cycles
for k2 = 1:set.nz
    
    % calculate the "real" system
    y_LTI = sys.G_S * u;
    
    % add noise to output
    if set.enable_noise
        
        y_LTI = y_LTI + set.mu + set.sigma*randn(size(y_LTI));
        
    end
    
    % save norm of error
    cycle_error_each_iter(k2) = norm(r - y_LTI);
    
    % save y
    y_iter_all(k2, :) = y_LTI;
    
    % save perfect error for comparision
    if k2 <= set.gUpdateCycle
        cycle_error_perfect(k2,3) = norm(r - y_LTI);
    else
        cycle_error_perfect(k2,3) = (1 - set.gamma)*cycle_error_perfect(k2-1,3);
    end
    
    
    
    %% model update
    if set.enable_error_update && k2 == set.gUpdateCycle
        
        % get the current and previous error
        e_current   = r - y_iter_all(k2,:)';
        e_prev      = r - y_iter_all(k2-1,:)';
        
        % save data
        data.e_current  = e_current;
        data.e_prev     = e_prev;
        data.y_LTI      = y_LTI;
        data.r          = r;
        data.u          = u;
        data.g_me_init  = g_m_all(:, 1) .* e_prev;
        data.g_m_init   = g_m_all(:, 1);
        data.gamma      = set.gamma;
        
        % initialise g_m and g_me
        g_m_update      = [];
        g_me_update     = [];
        
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
            
            % convert g_m_vector to G_M_matrix
            G_M_update = convert_from_vector_to_matrix( g_m_update, time );
            
            % perform update and save in data structure
            sys.G_M = G_M_update;
            
        elseif not(isempty(g_me_update))
            
            % transfer update to a matrix
            G_ME = convert_from_vector_to_matrix( g_me_update, time );
            
            % perform update and save in data structure
            sys.G_M =  1/set.gamma*(sys.G_M - G_ME);
            
            
        else
            disp('No model update applied.')
        end
        
        % save impulse response of update
        g_m_all(:, 3) = sys.G_M(:,1);
        
        % stability checks - estimate the residual of the model update
        if not(isempty(g_m_update))
            
            res = norm( G_M_prev * (e_current - e_prev) + set.gamma * G_M_update * e_prev );
            
        elseif not(isempty(g_me_update))
            
            res = norm( sys.G_M * e_current - G_ME * e_prev );
            
        end
        
        if set.enable_model_check && res >= set.res_threshold
            sys.G_M = G_M_prev;
            update_stopped = true;
            disp('Stop update due to residual violation.')
        end
    end
    
    % for each but the last cycle apply the update of the input trajectory
    if k2 < set.nz
        
        switch set.control_type
            
            % calulate the update using a quadratic form
            case 'quadprog'
                
                % options = [];
                options = optimset('Display','off','Algorithm','interior-point-convex');
                
                % weighting matrix for quadratic optimization
                Q = eye(time.nt,time.nt);
                
                f = -2*(r - y_LTI)'*Q*sys.G_M;
                H = 2*sys.G_M'*Q*sys.G_M;
                
                % make sure the Hessian is symmetric
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
                
            otherwise
                
                u_update = 0;
                warning('No control update applied.')
        end
        
        % check if cycle error has increased
        if not(update_stopped) && set.enable_stalibity_stop && k2 > 2 && (cycle_error_each_iter(k2) - cycle_error_each_iter(k2-1)) >= 0
            update_stopped = true;
            disp('Stop update due to error increase.')
        end
        
        if update_stopped
            u_update = 0;
        end
        
        % apply the update to the trajectory
        u = u + set.gamma * u_update;
        
    end
end


%% Outputs
out.g_m_all = g_m_all;
out.time = time;
out.sys = sys;
out.yIterAll = y_iter_all;
out.r = r;
out.set = set;
out.cycle_error_each_iter = cycle_error_each_iter;
out.cycle_error_perfect = cycle_error_perfect;



%% Visualisation
if set.enable_plots
    
    figure
    subplot(3,1,1)
    hold on
    plot(time.vec, y_iter_all(1,:,1),'b');
    plot(time.vec, y_iter_all(2,:,1),'m');
    plot(time.vec, y_iter_all(end,:,1),'k');
    plot(time.vec, r,'r-.');
    legend('1st Cycle','2nd Cycle', 'Final Cycle', 'Reference', 'Location','eastoutside')
    xlabel('Time')
    ylabel('Output')
    subplot(3,1,2)
    hold on
    plot(time.vec, r - y_iter_all(1,:)','b');
    plot(time.vec, r - y_iter_all(2,:)','m');
    plot(time.vec, r - y_iter_all(end,:)','k');
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