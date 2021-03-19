function all_out = master_multistart_v2


% close all
clear all
clc


n_runs = 50;
all_out = cell(n_runs, 1);

set = load_default_settings();

for k1 = 1:n_runs
    
    disp(['Start run ',num2str(k1),' of ',num2str(n_runs),'.'])
    
    
    % override specific settings for different runs here
    if k1 <= n_runs/2
        set.model_update_type = 'tfest_gm';
    else
        set.model_update_type = 'tfest_gm_io';
    end
    
    % noise settings
    set.mu_y = 0.5*(2*rand-1);
    set.sigma_y = 0.01*(2*rand-1);
    
    
    % perform simulation
    all_out{k1} = master_sim_v3(set);
end




% %% Visualization
% figure
% for k1 = 1:n_runs
%     subplot(3,1,1)
%     hold on
%     plot(all_out{k1}.time.vec, all_out{k1}.yIterAll(1,:),'b');
%     plot(all_out{k1}.time.vec, all_out{k1}.yIterAll(2,:),'m');
%     plot(all_out{k1}.time.vec, all_out{k1}.yIterAll(end,:),'k');
%     plot(all_out{k1}.time.vec, all_out{k1}.r,'r-.');
%     legend('1st Cycle','2nd Cycle', 'Final Cycle', 'Reference', 'Location','eastoutside')
%     xlabel('Time')
%     ylabel('Output')
%     subplot(3,1,2)
%     hold on
%     plot(all_out{k1}.time.vec, all_out{k1}.r - all_out{k1}.yIterAll(1,:)','b');
%     plot(all_out{k1}.time.vec, all_out{k1}.r - all_out{k1}.yIterAll(2,:)','m');
%     plot(all_out{k1}.time.vec, all_out{k1}.r - all_out{k1}.yIterAll(end,:)','k');
%     legend('1st Error','2nd Error', 'Final Error', 'Location','eastoutside')
%     xlabel('Time')
%     ylabel('Error')
%     subplot(3,1,3)
%     hold on
%     plot(all_out{k1}.time.vec, all_out{k1}.g_m_all(:,2),'m');
%     plot(all_out{k1}.time.vec, all_out{k1}.g_m_all(:,1),'b');
%     plot(all_out{k1}.time.vec, all_out{k1}.g_m_all(:,3),'k-.');
%     legend('True', 'Initial', 'Updated', 'Location','eastoutside')
%     xlabel('Time')
%     ylabel('Impuse Response')
% end
% 
% figure
% for k1 = 1:n_runs
%     grid on
%     hold on
%     plot(all_out{k1}.time.vec, all_out{k1}.g_m_all(:,2),'m');
%     plot(all_out{k1}.time.vec, all_out{k1}.g_m_all(:,1),'b');
%     plot(all_out{k1}.time.vec, all_out{k1}.g_m_all(:,3),'k-.');
%     legend('True', 'Initial', 'Updated', 'Location','eastoutside')
%     xlabel('Time')
%     ylabel('Impuse Response')
% end
% 
% 
% figure
% for k1 = 1:n_runs
%     grid on
%     hold on
%     if k1 <= n_runs/2
%     plot(all_out{k1}.time.vec, all_out{k1}.g_m_all(:,2) - all_out{k1}.g_m_all(:,3),'r');
%     else
%         plot(all_out{k1}.time.vec, all_out{k1}.g_m_all(:,2) - all_out{k1}.g_m_all(:,3),'b');
%     end
%     % legend('Difference', 'Location','eastoutside')
%     xlabel('Time')
%     ylabel('Impuse Response Error')
% end
% 
% 
% figure
% hold on
% grid on
% for k1 = 1:n_runs
%     plot(all_out{k1}.cycle_error_each_iter, 'r-*')
%     plot(all_out{k1}.cycle_error_perfect, 'mo')
% end
% ax = gca;
% ax.YScale = 'log';
% legend('Model Update', 'Model Init', 'Perfect', 'Location','southwest')
% xlabel('Cycle')
% ylabel('Error in L2-Norm')





end