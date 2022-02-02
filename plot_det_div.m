function plot_det_div(exp_i,type,mi,model_pred,full_figure_flag, save_flag)
% Detection % if type ==1 or 3
%exp_i = 2;
%exp_i = 2;
%type = 1;
%full_figure_flag = 1;
%save_flag = 1;
%mi = 1;
%mi = 2;
load(['all_vars_exp_',num2str(exp_i),'_type_',num2str(type) ,'_mi_',num2str(mi), '.mat'])

set_plotting_params();

figure(1)
if full_figure_flag
    set(gcf, 'Position', [100 100 500 820])%set(gcf, 'Position', [100 100 500 920])
    guttera2 = guttera;
    marginsa2 = [0.14 0.14 0.08 0.08]; %marginsa; %left right bottom top
    
else
    set(gcf, 'Position', [100 100 500 190])
    guttera2 = guttera;
    marginsa2 = [0.14 0.14 0.2 0.2];
end




ticklengthsz = [0.025 0.025];


for cond = [1 3]
    if full_figure_flag
        
        % a)   set size
        %tight_subplot(5,2,1,0.5+cond/2, guttera2, marginsa2)
        tight_subplot(4,2,1,0.5+cond/2, guttera2, marginsa2)
        if model_pred
            h_p=fill([Nvec Nvec(end:-1:1)], [[squeeze(mean(p_pr_tp(:,cond,:),1))]'-[squeeze(std(p_pr_tp(:,cond,:),1))]'/sqrt(Nsubj)...
                [squeeze(mean(p_pr_tp(:,cond,end:-1:1),1))]'+[squeeze(std(p_pr_tp(:,cond,:,end:-1:1),1))]'/sqrt(Nsubj)],bluee_shade,'EdgeColor', 'None'); hold on;
            h_a=fill([Nvec Nvec(end:-1:1)], [[squeeze(mean(p_pr_ta(:,cond,:),1))]'-[squeeze(std(p_pr_ta(:,cond,:),1))]'/sqrt(Nsubj)...
                [squeeze(mean(p_pr_ta(:,cond,end:-1:1),1))]'+[squeeze(std(p_pr_ta(:,cond,:,end:-1:1),1))]'/sqrt(Nsubj)],redd_shade ,'EdgeColor', 'None'); hold on;
            
        end
        if model_pred %no line
            plot(Nvec, squeeze(mean(p_tp(:, cond, :),1)), 'Color', bluee, 'Linestyle', 'none'); hold on;
            errorbar(Nvec, squeeze(mean(p_tp(:, cond, :),1)),squeeze(std(p_tp(:, cond, :),1))/sqrt(Nsubj), 'o','MarkerSize', sz_dot, 'MarkerFaceColor', bluee, 'MarkerEdgeColor', bluee,'Color',bluee , 'LineWidth',1.3); hold on;
            plot(Nvec, squeeze(mean(p_ta(:, cond, :),1)), 'Color', redd, 'Linestyle', 'none'); hold on;
            errorbar(Nvec, squeeze(mean(p_ta(:, cond, :),1)),squeeze(std(p_ta(:, cond, :),1))/sqrt(Nsubj), 'o','MarkerSize', sz_dot, 'MarkerFaceColor', redd, 'MarkerEdgeColor',redd,'Color',redd, 'LineWidth',1.3); hold on;
        else % line
            plot(Nvec, squeeze(mean(p_tp(:, cond, :),1)), 'Color', bluee); hold on;
            errorbar(Nvec, squeeze(mean(p_tp(:, cond, :),1)),squeeze(std(p_tp(:, cond, :),1))/sqrt(Nsubj), 'o','MarkerSize', sz_dot, 'MarkerFaceColor', bluee, 'MarkerEdgeColor', bluee,'Color',bluee , 'LineWidth',1.3); hold on;
            plot(Nvec, squeeze(mean(p_ta(:, cond, :),1)), 'Color', redd); hold on;
            errorbar(Nvec, squeeze(mean(p_ta(:, cond, :),1)),squeeze(std(p_ta(:, cond, :),1))/sqrt(Nsubj), 'o','MarkerSize', sz_dot, 'MarkerFaceColor', redd, 'MarkerEdgeColor',redd,'Color',redd, 'LineWidth',1.3); hold on;
        end
        
        chance=plot(Nvec, 0.5*ones(1, length(Nvec)),'--k'); hold on;
        box off
        ylim([0.3 1])
        xlim([min(Nvec)-0.5 max(Nvec)+0.5]);
        box off
        set(gca, 'tickdir', 'out')
        set(gca, 'xtick', Nvec)
        set(gca, 'xticklabels', {'2', '3', '4', '6'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
        set(gca, 'ytick', [0:0.1:1])
        set(gca, 'yticklabels', {})
        set(gca,'ticklength', ticklengthsz)
        xlabel('Set size','FontName', 'Helvetica', 'FontSize', fontsz)
        if cond == 1
            ylabel('Proportion target present','FontName', 'Helvetica', 'FontSize', fontsz)
            set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
        end
        text(3, 1.16, text_labels{(1+cond)/2},  'FontName', 'Helvetica', 'FontSize', fontsz)
        
        
        %b) T-most similar of distractors
        
        for Nind = 1:4
            bincentersl=[squeeze(bincenterz(cond, Nind,1,:))]';
            
            %tight_subplot(5,2,2,0.5+cond/2, guttera2, marginsa2)
            tight_subplot(4,2,2,0.5+cond/2, guttera2, marginsa2)
            
            if model_pred
                hacc=fill([bincentersl bincentersl(end:-1:1)], [[squeeze(mean(accuracy_pred_msd_all(:,cond,Nind,:),1))]'-[squeeze(std(accuracy_pred_msd_all(:,cond,Nind,:),1))./sqrt(squeeze(sum(~isnan(accuracy_pred_msd_all(:,cond,Nind,:)),1)))]'...
                    [squeeze(mean(accuracy_pred_msd_all(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(accuracy_pred_msd_all(:,cond,Nind,end:-1:1),1))./sqrt(squeeze(sum(~isnan(accuracy_pred_msd_all(:,cond,Nind,:)),1)))]'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
                
            end
        end
        for Nind = 1:4
            bincentersl=[squeeze(bincenterz(cond, Nind,1,:))]';
            
            %tight_subplot(5,2,2,0.5+cond/2, guttera2, marginsa2)
            tight_subplot(4,2,2,0.5+cond/2, guttera2, marginsa2)
            if model_pred
                plot(bincentersl, squeeze(mean(accuracy_msd_all(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
                leg(Nind)=errorbar(bincentersl,squeeze(mean(accuracy_msd_all(:,cond,Nind,:),1)), squeeze(std(accuracy_msd_all(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none',  'Linewidth', 1.3);
            else
                plot(bincentersl, squeeze(mean(accuracy_msd_all(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
                leg(Nind)=errorbar(bincentersl,squeeze(mean(accuracy_msd_all(:,cond,Nind,:),1)), squeeze(std(accuracy_msd_all(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:),  'Linewidth', 1.3);
            end
            %plot(bincentersl,squeeze(mean(p_data_tp(:, cond, Nind, :),1)), '-o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',3*sz_dot); hold on;
            %plot(bincentersl,squeeze(mean(p_data_ta(:, cond, Nind, :),1)), '--o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',3*sz_dot); hold on;
            chance = plot(0:10:90, 0.5*ones(1, 10),'--k'); hold on;
            box off
            ylim([0.4 1])
            xlim([0 90])
            set(gca, 'xtick', [0:15:90])
            set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
            set(gca, 'ytick', [0:0.1:1])
            set(gca, 'yticklabels', {})
            if cond == 1
                set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
                ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', ticklengthsz)
            xlabel('min T-D difference (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
            
            
            if cond == 4
                llpa = legend([leg(1) leg(2) leg(3) leg(4) chance],'N = 2', 'N = 3','N = 4', 'N = 6','chance','FontName','Helvetica','FontSize', fontsz)
                legend boxoff
                set(llpa, 'FontName','Helvetica','FontSize', fontsz*0.9, 'Position', [0.9 0.84 0.04 0.04])
            end
            
            
            
            
            %c) T-mean of distractors
            
            binzz_meanE = 90/pi* binzz_meanE;  %since we doubled all our stimuli for modelling on the (-pi,pi) full space
            
            
            for Nind = 2:4
                
                binzz_mean_mean = [mean(squeeze(binzz_meanE(:,cond,Nind,:)),1)]';
                bincenterz_mean = [(binzz_mean_mean(1:end-1)+binzz_mean_mean(2:end))/2]';
                
                %tight_subplot(5,2,3,0.5+cond/2, guttera2, marginsa2)
                tight_subplot(4,2,3,0.5+cond/2, guttera2, marginsa2)
                
                if model_pred
                    h_m(Nind)=fill([bincenterz_mean bincenterz_mean(end:-1:1)], [[squeeze(mean(accuracy_pred_mean_all(:,cond,Nind,:),1))]'-[squeeze(std(accuracy_pred_mean_all(:,cond,Nind,:),1))]'./sqrt(squeeze(sum(~isnan(accuracy_pred_mean_all(:,cond,Nind,:)),1)))'...
                        [squeeze(mean(accuracy_pred_mean_all(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(accuracy_pred_mean_all(:,cond,Nind,end:-1:1),1))]'./sqrt(squeeze(sum(~isnan(accuracy_pred_mean_all(:,cond,Nind,:)),1)))'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
                    
                end
            end
            
            for Nind = 2:4
                
                binzz_mean_mean = [mean(squeeze(binzz_meanE(:,cond,Nind,:)),1)]';
                bincenterz_mean = [(binzz_mean_mean(1:end-1)+binzz_mean_mean(2:end))/2]';
                
                %tight_subplot(5,2,3,0.5+cond/2, guttera2, marginsa2)
                tight_subplot(4,2,3,0.5+cond/2, guttera2, marginsa2)
                if model_pred
                    plot(bincenterz_mean, squeeze(mean(accuracy_mean_all(:,cond,Nind,:),1)),'Color', colors_sz(Nind,:),'Linestyle', 'none','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
                    leg(Nind) = errorbar(bincenterz_mean,squeeze(mean(accuracy_mean_all(:,cond,Nind,:),1)), squeeze(std(accuracy_mean_all(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none',  'Linewidth', 1.3);
                else
                    plot(bincenterz_mean, squeeze(mean(accuracy_mean_all(:,cond,Nind,:),1)),'Color', colors_sz(Nind,:),'Linestyle', 'none','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
                    leg(Nind) = errorbar(bincenterz_mean,squeeze(mean(accuracy_mean_all(:,cond,Nind,:),1)), squeeze(std(accuracy_mean_all(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:),  'Linewidth', 1.3);
                end
                %plot(bincenterz_mean,squeeze(mean(p_data_tp_mean(:, cond, Nind, :),1)), '-o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',3*sz_dot); hold on;
                %plot(bincenterz_mean,squeeze(mean(p_data_ta_mean(:, cond, Nind, :),1)), '--o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',3*sz_dot); hold on;
                
                
                chance = plot(0:10:90, 0.5*ones(1, 10),'--k'); hold on;
                box off
                ylim([0.4 1])
                xlim([0 90])
                set(gca, 'xtick', [0:15:90])
                set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
                set(gca, 'ytick', [0:0.1:1])
                set(gca, 'yticklabels', {})
                set(gca, 'tickdir', 'out')
                set(gca,'ticklength', ticklengthsz)
                xlabel('T-D mean (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
                
                if cond == 1
                    set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
                    ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
                end
                
            end
            
            binzz_meanE = binzz_meanE *pi/90;
            
            
            % d) circ var of distractors
            for Nind = 2:4
                
                binzz_cvar_mean = [mean(squeeze(binzz_cvarE(:,cond,Nind,:)),1)]';
                bincenterz_cvar = [(binzz_cvar_mean(1:end-1)+binzz_cvar_mean(2:end))/2]';
                
                %tight_subplot(5,2,4,0.5+cond/2, guttera2, marginsa2)
                tight_subplot(4,2,4,0.5+cond/2, guttera2, marginsa2)
                if model_pred
                    h_v(Nind)=fill([bincenterz_cvar bincenterz_cvar(end:-1:1)], [[squeeze(nanmean(accuracy_pred_cvar_all(:,cond,Nind,:),1))]'-[squeeze(nanstd(accuracy_pred_cvar_all(:,cond,Nind,:),1))]'./sqrt(squeeze(sum(~isnan(accuracy_pred_cvar_all(:,cond,Nind,:)),1)))'...
                        [squeeze(nanmean(accuracy_pred_cvar_all(:,cond,Nind,end:-1:1),1))]'+[squeeze(nanstd(accuracy_pred_cvar_all(:,cond,Nind,end:-1:1),1))]'./sqrt(squeeze(sum(~isnan(accuracy_pred_cvar_all(:,cond,Nind,:)),1)))'],colors_sz_shade(Nind,:),'EdgeColor', 'None', 'Linestyle', 'none'); hold on;
                    
                end
                
            end
            for Nind = 2:4
                
                binzz_cvar_mean = [mean(squeeze(binzz_cvarE(:,cond,Nind,:)),1)]';
                bincenterz_cvar = [(binzz_cvar_mean(1:end-1)+binzz_cvar_mean(2:end))/2]';
                
                %tight_subplot(5,2,4,0.5+cond/2, guttera2, marginsa2)
                tight_subplot(4,2,4,0.5+cond/2, guttera2, marginsa2)
                if model_pred
                    plot(bincenterz_cvar, squeeze(nanmean(accuracy_cvar_all(:,cond,Nind,:),1)),  'o', 'Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
                    leg(Nind) = errorbar(bincenterz_cvar,squeeze(nanmean(accuracy_cvar_all(:,cond,Nind,:),1)), squeeze(nanstd(accuracy_cvar_all(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:),'Linestyle', 'none',  'Linewidth', 1.3);
                else
                    plot(bincenterz_cvar, squeeze(nanmean(accuracy_cvar_all(:,cond,Nind,:),1)),  'o', 'Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
                    leg(Nind) = errorbar(bincenterz_cvar,squeeze(nanmean(accuracy_cvar_all(:,cond,Nind,:),1)), squeeze(nanstd(accuracy_cvar_all(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:),  'Linewidth', 1.3);
                end
                chance = plot(0:0.10:0.90, ones(1,10)*1/Nvec(Nind),'--', 'Color',colors_sz(Nind,:)); hold on;
                box off
                ylim([0.4 1])
                xlim([0 1])
                set(gca, 'xtick',[0:0.2:1])
                set(gca, 'ytick', [0:0.1:1])
                set(gca, 'yticklabels', {})
                set(gca, 'tickdir', 'out')
                set(gca,'ticklength', ticklengthsz)
                xlabel('Distractor variance','FontName', 'Helvetica', 'FontSize', fontsz)
                set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
                if cond == 1
                    
                    ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
                end
            end
            
            %{
            % e) circ var of distractors, false alarms
            for Nind = 2:4
                
                binzz_cvar_mean = [mean(squeeze(binzz_cvarE(:,cond,Nind,:)),1)]';
                bincenterz_cvar = [(binzz_cvar_mean(1:end-1)+binzz_cvar_mean(2:end))/2]';
                
                tight_subplot(5,2,5,0.5+cond/2, guttera2, marginsa2)
                
                if model_pred
                    h_v(Nind)=fill([bincenterz_cvar bincenterz_cvar(end:-1:1)], [[squeeze(mean(p_pred_ta_cvar(:,cond,Nind,:),1))]'-[squeeze(std(p_pred_ta_cvar(:,cond,Nind,:),1))]'./sqrt(squeeze(sum(~isnan(p_pred_ta_cvar(:,cond,Nind,:)),1)))'...
                        [squeeze(nanmean(p_pred_ta_cvar(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(p_pred_ta_cvar(:,cond,Nind,end:-1:1),1))]'./sqrt(squeeze(sum(~isnan(p_pred_ta_cvar(:,cond,Nind,:)),1)))'],colors_sz_shade(Nind,:),'EdgeColor', 'None', 'Linestyle', 'none'); hold on;
                    
                end
                
            end
            for Nind = 2:4
                
                binzz_cvar_mean = [mean(squeeze(binzz_cvarE(:,cond,Nind,:)),1)]';
                bincenterz_cvar = [(binzz_cvar_mean(1:end-1)+binzz_cvar_mean(2:end))/2]';
                
                tight_subplot(5,2,5,0.5+cond/2, guttera2, marginsa2)
                if model_pred
                    plot(bincenterz_cvar, squeeze(mean(p_data_ta_cvar(:,cond,Nind,:),1)),  'o', 'Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
                    leg(Nind) = errorbar(bincenterz_cvar,squeeze(mean(p_data_ta_cvar(:,cond,Nind,:),1)), squeeze(std(p_data_ta_cvar(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:),'Linestyle', 'none');
                else
                    plot(bincenterz_cvar, squeeze(mean(p_data_ta_cvar(:,cond,Nind,:),1)),  'o', 'Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
                    leg(Nind) = errorbar(bincenterz_cvar,squeeze(mean(p_data_ta_cvar(:,cond,Nind,:),1)), squeeze(std(p_data_ta_cvar(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:));
                end
                %chance = plot(0:0.10:0.90, ones(1,10)*1/Nvec(Nind),'--', 'Color',colors_sz(Nind,:)); hold on;
                box off
                ylim([0.2 1])
                xlim([0 1])
                set(gca, 'xtick',[0:0.2:1])
                set(gca, 'ytick', [0:0.1:1])
                set(gca, 'yticklabels', {})
                set(gca, 'tickdir', 'out')
                set(gca,'ticklength', [0.04 0.04])
                xlabel('D variance','FontName', 'Helvetica', 'FontSize', fontsz)
                if cond == 1
                    set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
                    ylabel('False-alarm rate','FontName', 'Helvetica', 'FontSize', fontsz)
                end
            end
            %}
            
            
            
        end
    else
        %b) T-most similar of distractors
        
        for Nind = 1:4
            bincentersl=[squeeze(bincenterz(cond, Nind,1,:))]';
            
            tight_subplot(1,2,1,0.5+cond/2, guttera2, marginsa2)
            
            if model_pred
                hacc=fill([bincentersl bincentersl(end:-1:1)], [[squeeze(mean(accuracy_pred_msd_all(:,cond,Nind,:),1))]'-[squeeze(std(accuracy_pred_msd_all(:,cond,Nind,:),1))./sqrt(squeeze(sum(~isnan(accuracy_pred_msd_all(:,cond,Nind,:)),1)))]'...
                    [squeeze(mean(accuracy_pred_msd_all(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(accuracy_pred_msd_all(:,cond,Nind,end:-1:1),1))./sqrt(squeeze(sum(~isnan(accuracy_pred_msd_all(:,cond,Nind,:)),1)))]'],colors_sz_shade(Nind,:),'EdgeColor', 'None', 'Linestyle', 'none'); hold on;
                
            end
        end
        for Nind = 1:4
            bincentersl=[squeeze(bincenterz(cond, Nind,1,:))]';
            
            tight_subplot(1,2,1,0.5+cond/2, guttera2, marginsa2)
            if model_pred
                plot(bincentersl, squeeze(mean(accuracy_msd_all(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
                leg(Nind)=errorbar(bincentersl,squeeze(mean(accuracy_msd_all(:,cond,Nind,:),1)), squeeze(std(accuracy_msd_all(:,cond,Nind,:),1))/sqrt(Nsubj),'Color', colors_sz(Nind,:),'Linestyle', 'none',  'Linewidth', 1.3);
            else
                plot(bincentersl, squeeze(mean(accuracy_msd_all(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
                leg(Nind)=errorbar(bincentersl,squeeze(mean(accuracy_msd_all(:,cond,Nind,:),1)), squeeze(std(accuracy_msd_all(:,cond,Nind,:),1))/sqrt(Nsubj),'Color', colors_sz(Nind,:),  'Linewidth', 1.3);
            end
            chance = plot(0:10:90, 0.5*ones(1, 10),'--k'); hold on;
            
            box off
            ylim([0.4 1])
            xlim([0 90])
            set(gca, 'xtick', [0:15:90])
            set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
            set(gca, 'ytick', [0:0.1:1])
            set(gca, 'yticklabels', {})
            if cond == 1
                set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
                ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', ticklengthsz)
            xlabel('min T-D difference (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
            if Nind ==1
                text(30, 1.16, text_labels{(1+cond)/2},  'FontName', 'Helvetica', 'FontSize', fontsz)
            end
            if cond == 1 & Nind == 1
                text(-40, 1.13, 'Detection',  'FontName', 'Helvetica', 'FontSize', fontsz*1.2)
            end
            
            if cond == 4
                llpa = legend([leg(1) leg(2) leg(3) leg(4) chance],'N = 2', 'N = 3','N = 4', 'N = 6','chance','FontName','Helvetica','FontSize', fontsz)
                legend boxoff
                set(llpa, 'FontName','Helvetica','FontSize', fontsz*0.9, 'Position', [0.9 0.84 0.04 0.04])
            end
            
        end
        
    end
end

psname = ['DET_full_figure_',num2str(full_figure_flag),'_model_',num2str(model_pred),'_model_index_',num2str(mi),'_type_',num2str(type) ,'_exp_',num2str(exp_i), '.pdf']
if save_flag
    print_pdf(psname)
end
%% D & H direct test visualization
%close all;
figure(2)
set(gcf, 'Position', [100 100 500 820])
set_plotting_params();
colorsz_mean =  [255 165 0; 154 187 246;  204 221 123 ]./255;%[255 165 0; 154 217 246;  0 0 0 ]./255; %[154 217 246;  0 0 0; 255 165 0]./255;

guttera3 = [0.06 0.04];
marginsa3 = [0.14 0.14 0.08 0.08]; %marginsa; %left right bottom top


load('mean_min_cvar_hist_type_1.mat')


for cond = [1 3]
    for Nind = 2:4
        binz_mean_2_medians =  squeeze(binzz_mean_2_E(Nind,:));
        binz_cvar_3_medians =  squeeze(binzz_cvar_3_E(Nind,:));
        bincenters_mean_2_medians(Nind,:) = (binz_mean_2_medians(2:end)+ binz_mean_2_medians(1:end-1))/2;
        bincenters_cvar_3_medians(Nind,:) = (binz_cvar_3_medians(2:end)+ binz_cvar_3_medians(1:end-1))/2;
        
        %tight_subplot(6,4,1.5*cond-0.5,1, guttera3, marginsa3)
        % all combined
        
        
        tight_subplot(6,3,1.5*cond-0.5,Nind-1, guttera3, marginsa3)
        for jj = 1:2
            plot(bincenters_cvar_3_medians(Nind,:), mean(squeeze(1-error_rate_mean_cvar_all(:, cond, Nind, jj,:))), 'o','MarkerEdgeColor', colorsz_mean(jj,:), 'MarkerFaceColor', colorsz_mean(jj,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
            leg(jj)=errorbar(bincenters_cvar_3_medians(Nind,:),mean(squeeze(1-error_rate_mean_cvar_all(:, cond, Nind, jj,:))), std(squeeze(1-error_rate_mean_cvar_all(:, cond, Nind, jj,:)))/sqrt(Nsubj),'Color', colorsz_mean(jj,:),  'Linewidth', 1.3, 'CapSize', 0); hold on;
        end
        % jj = 2, but for ONLY 1 and 2 are the same values
        plot(bincenters_cvar_3_medians(Nind,:), mean(squeeze(1-error_rate_mean_cvar_all_ONLY(:, cond, Nind, jj,:))), 'o','MarkerEdgeColor', colorsz_mean(3,:), 'MarkerFaceColor', colorsz_mean(3,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
        leg_all= errorbar(bincenters_cvar_3_medians(Nind,:),mean(squeeze(1-error_rate_mean_cvar_all_ONLY(:, cond, Nind, jj,:))), std(squeeze(1-error_rate_mean_cvar_all_ONLY(:, cond, Nind, jj,:)))/sqrt(Nsubj),'Color', colorsz_mean(3,:),  'Linewidth', 1.6, 'CapSize', 0); hold on;
        xlim([0 1]); ylim([0.4 1]); box off
        set(gca, 'tickdir', 'out')
        
        
        if Nind == 2
            ylabel('Accuracy');
        end
        if Nind == 2 & cond == 1
            text(0.4, 1.1, 'N = 3')
            llpa = legend([leg(1) leg(2) leg_all],'T-D mean < median', 'T-D mean > median', 'all', 'FontName','Helvetica','FontSize', fontsz)
            legend boxoff
            set(llpa, 'FontName','Helvetica','FontSize', fontsz*0.9, 'Position', [0.17 0.94 0.04 0.04])
        end
        
        if Nind == 3 & cond == 1
            text(0.4, 1.1, 'N = 4')
        end
        if Nind == 4 & cond == 1
            text(0.4, 1.1, 'N = 6')
        end
        
        tight_subplot(6,3,1.5*cond+0.5,Nind-1, guttera3, marginsa3)
        for jj = 1:2
            plot(bincenters_cvar_3_medians(Nind,:), mean(squeeze(p_data_tp_mean_cvar(:, cond, Nind, jj,:))), 'o','MarkerEdgeColor', colorsz_mean(jj,:), 'MarkerFaceColor', colorsz_mean(jj,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
            leg(Nind)=errorbar(bincenters_cvar_3_medians(Nind,:), mean(squeeze(p_data_tp_mean_cvar(:, cond, Nind, jj,:))), std(squeeze(p_data_tp_mean_cvar(:, cond, Nind, jj,:)))/sqrt(Nsubj),'Color', colorsz_mean(jj,:),  'Linewidth', 1.3, 'CapSize', 0);
        end
        plot(bincenters_cvar_3_medians(Nind,:), mean(squeeze(p_data_tp_mean_cvar_ONLY(:, cond, Nind, jj,:))), 'o','MarkerEdgeColor', colorsz_mean(3,:), 'MarkerFaceColor', colorsz_mean(3,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
        errorbar(bincenters_cvar_3_medians(Nind,:),mean(squeeze(p_data_tp_mean_cvar_ONLY(:, cond, Nind, jj,:))), std(squeeze(p_data_tp_mean_cvar_ONLY(:, cond, Nind, jj,:)))/sqrt(Nsubj),'Color', colorsz_mean(3,:),  'Linewidth', 1.6, 'CapSize', 0);
        
        xlim([0 1]); ylim([0.4 1]); box off
        set(gca, 'tickdir', 'out')
        if Nind == 2
            ylabel('Hit rate');
        end
        
        tight_subplot(6,3,1.5*cond+1.5,Nind-1, guttera3, marginsa3)
        for jj = 1:2
            plot(bincenters_cvar_3_medians(Nind,:), mean(squeeze(p_data_ta_mean_cvar(:, cond, Nind, jj,:))), 'o','MarkerEdgeColor', colorsz_mean(jj,:), 'MarkerFaceColor', colorsz_mean(jj,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
            leg(Nind)=errorbar(bincenters_cvar_3_medians(Nind,:), mean(squeeze(p_data_ta_mean_cvar(:, cond, Nind, jj,:))), std(squeeze(p_data_ta_mean_cvar(:, cond, Nind, jj,:)))/sqrt(Nsubj),'Color', colorsz_mean(jj,:),  'Linewidth', 1.3, 'CapSize',0);
        end
        plot(bincenters_cvar_3_medians(Nind,:), mean(squeeze(p_data_ta_mean_cvar_ONLY(:, cond, Nind, jj,:))), 'o','MarkerEdgeColor', colorsz_mean(3,:), 'MarkerFaceColor', colorsz_mean(3,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
        errorbar(bincenters_cvar_3_medians(Nind,:), mean(squeeze(p_data_ta_mean_cvar_ONLY(:, cond, Nind, jj,:))), std(squeeze(p_data_ta_mean_cvar_ONLY(:, cond, Nind, jj,:)))/sqrt(Nsubj),'Color', colorsz_mean(3,:),  'Linewidth', 1.3, 'CapSize',0);
        
        xlim([0 1]); ylim([0 1]); box off
        set(gca, 'tickdir', 'out')
        if cond == 3
            xlabel('Distractor variance');
        end
        if Nind == 2
            ylabel('False-alarm rate');
        end
        
    end
    if save_flag
        psname = ['DD_variance_plot_mean_cvar_DET_exp_',num2str(exp_i),'.pdf']
        print_pdf(psname)
    end
end
%% HR AND FA breakdown
%close all;
figure(3)
set(gcf, 'Position', [100 100 500 820])
set_plotting_params();
colorsz_mean =  [255 165 0; 154 187 246;  204 221 123 ]./255;%[255 165 0; 154 217 246;  0 0 0 ]./255; %[154 217 246;  0 0 0; 255 165 0]./255;

guttera3 = [0.10 0.03]; %[0.06 0.04];
marginsa3 = [0.14 0.14 0.08 0.08]; %marginsa; %left right bottom top


model_pred = 1;


fontsz = 11; % vs 13
face_alpha_val = 0.7;

for cond = [1 3]
    % 6,2,1 HR MIN T-D DIFF
    for Nind = 1:4
        bincentersl=[squeeze(bincenterz(cond, Nind,1,:))]';
        
        tight_subplot(14,2,[ 1 2],0.5+cond/2, guttera3, marginsa2)
        
        if model_pred
            hacc(Nind)=fill([bincentersl bincentersl(end:-1:1)], [[squeeze(mean(p_pred_tp(:,cond,Nind,:),1))]'-[squeeze(std(p_pred_tp(:,cond,Nind,:),1))./sqrt(squeeze(sum(~isnan(p_pred_tp(:,cond,Nind,:)),1)))]'...
                [squeeze(mean(p_pred_tp(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(p_pred_tp(:,cond,Nind,end:-1:1),1))./sqrt(squeeze(sum(~isnan(p_pred_tp(:,cond,Nind,:)),1)))]'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
            set(hacc(Nind), 'FaceAlpha', face_alpha_val)
        end
    end
    for Nind = 1:4
        bincentersl=[squeeze(bincenterz(cond, Nind,1,:))]';
        
        
        tight_subplot(14,2,[1 2],0.5+cond/2, guttera3, marginsa2)
        if model_pred
            plot(bincentersl, squeeze(mean(p_data_tp(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
            leg(Nind)=errorbar(bincentersl,squeeze(mean(p_data_tp(:,cond,Nind,:),1)), squeeze(std(p_data_tp(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none',  'Linewidth', 1.3);
        else
            plot(bincentersl, squeeze(mean(p_data_tp(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
            leg(Nind)=errorbar(bincentersl,squeeze(mean(p_data_tp(:,cond,Nind,:),1)), squeeze(std(p_data_tp(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:),  'Linewidth', 1.3);
        end
        %chance = plot(0:10:90, 0.5*ones(1, 10),'--k'); hold on;
        box off
        ylim([0.3 1])
        xlim([0 90])
        set(gca, 'xtick', [0:15:90])
        set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
        set(gca, 'ytick', [0:0.1:1])
        set(gca, 'yticklabels', {})
        if cond == 1
            set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
            ylabel('Hit rate','FontName', 'Helvetica', 'FontSize', fontsz)
        end
        set(gca, 'tickdir', 'out')
        set(gca,'ticklength', ticklengthsz)
        %xlabel('min T-D difference (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
        if Nind == 1 & cond == 1
            text(30,1.4,'Perception', 'FontName', 'Helvetica', 'FontSize', fontsz)
            text(-27, 1.2, 'A', 'FontName', 'Helvetica', 'FontSize', fontsz)
        elseif Nind == 1 &  cond == 3
            text(30,1.4,'Memory', 'FontName', 'Helvetica', 'FontSize', fontsz)
        end
    end
    if cond == 1
        llpa = legend([leg(1) leg(2) leg(3) leg(4)],'N = 2', 'N = 3','N = 4', 'N = 6','FontName','Helvetica','FontSize', fontsz)
        legend boxoff
        set(llpa, 'FontName','Helvetica','FontSize', fontsz*0.9, 'Position', [0.9 0.84 0.04 0.04])
    end
    
    
    % 6,2,2 FA MIN T-D DIFF
    for Nind = 1:4
        bincentersl=[squeeze(bincenterz(cond, Nind,1,:))]';
        
        tight_subplot(14,2,[ 3 4],0.5+cond/2, guttera3, marginsa2)
        
        if model_pred
            hacc(Nind)=fill([bincentersl bincentersl(end:-1:1)], [[squeeze(mean(p_pred_ta(:,cond,Nind,:),1))]'-[squeeze(std(p_pred_ta(:,cond,Nind,:),1))./sqrt(squeeze(sum(~isnan(p_pred_ta(:,cond,Nind,:)),1)))]'...
                [squeeze(mean(p_pred_ta(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(p_pred_ta(:,cond,Nind,end:-1:1),1))./sqrt(squeeze(sum(~isnan(p_pred_ta(:,cond,Nind,:)),1)))]'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
            set(hacc(Nind), 'FaceAlpha', face_alpha_val);
        end
    end
    for Nind = 1:4
        bincentersl=[squeeze(bincenterz(cond, Nind,1,:))]';
        
        
        tight_subplot(14,2,[3 4],0.5+cond/2, guttera3, marginsa2)
        if model_pred
            plot(bincentersl, squeeze(mean(p_data_ta(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
            leg(Nind)=errorbar(bincentersl,squeeze(mean(p_data_ta(:,cond,Nind,:),1)), squeeze(std(p_data_ta(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none',  'Linewidth', 1.3);
        else
            plot(bincentersl, squeeze(mean(p_data_ta(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
            leg(Nind)=errorbar(bincentersl,squeeze(mean(p_data_ta(:,cond,Nind,:),1)), squeeze(std(p_data_ta(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:),  'Linewidth', 1.3);
        end
        %chance = plot(0:10:90, 0.5*ones(1, 10),'--k'); hold on;
        box off
        ylim([0.0 1])
        xlim([0 90])
        set(gca, 'xtick', [0:15:90])
        set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
        set(gca, 'ytick', [0:0.1:1])
        set(gca, 'yticklabels', {})
        if cond == 1
            set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
            ylabel('False-alarm rate','FontName', 'Helvetica', 'FontSize', fontsz)
        end
        set(gca, 'tickdir', 'out')
        set(gca,'ticklength', ticklengthsz)
        xlabel('min T-D difference (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
    end
    
    
    binzz_meanE = 90/pi* binzz_meanE;
    
    % 6,2,3 HR  T-D MEAN
    for Nind = 1:4
        binzz_mean_mean = [mean(squeeze(binzz_meanE(:,cond,Nind,:)),1)]';
        bincenterz_mean = [(binzz_mean_mean(1:end-1)+binzz_mean_mean(2:end))/2]';
        
        tight_subplot(14,2,[ 6 7],0.5+cond/2, guttera3, marginsa2)
        
        if model_pred
            hacc(Nind)=fill([bincenterz_mean bincenterz_mean(end:-1:1)], [[squeeze(mean(p_pred_tp_mean(:,cond,Nind,:),1))]'-[squeeze(std(p_pred_tp_mean(:,cond,Nind,:),1))./sqrt(squeeze(sum(~isnan(p_pred_tp_mean(:,cond,Nind,:)),1)))]'...
                [squeeze(mean(p_pred_tp_mean(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(p_pred_tp_mean(:,cond,Nind,end:-1:1),1))./sqrt(squeeze(sum(~isnan(p_pred_tp_mean(:,cond,Nind,:)),1)))]'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
            set(hacc(Nind), 'FaceAlpha', face_alpha_val)
        end
    end
    for Nind = 1:4
        binzz_mean_mean = [mean(squeeze(binzz_meanE(:,cond,Nind,:)),1)]';
        bincenterz_mean = [(binzz_mean_mean(1:end-1)+binzz_mean_mean(2:end))/2]';
        
        
        tight_subplot(14,2,[6 7],0.5+cond/2, guttera3, marginsa2)
        if model_pred
            plot(bincenterz_mean, squeeze(mean(p_data_tp_mean(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
            leg(Nind)=errorbar(bincenterz_mean,squeeze(mean(p_data_tp_mean(:,cond,Nind,:),1)), squeeze(std(p_data_tp_mean(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none',  'Linewidth', 1.3);
        else
            plot(bincenterz_mean, squeeze(mean(p_data_tp_mean(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
            leg(Nind)=errorbar(bincenterz_mean,squeeze(mean(p_data_tp_mean(:,cond,Nind,:),1)), squeeze(std(p_data_tp_mean(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:),  'Linewidth', 1.3);
        end
        %chance = plot(0:10:90, 0.5*ones(1, 10),'--k'); hold on;
        box off
        ylim([0.3 1])
        xlim([0 90])
        set(gca, 'xtick', [0:15:90])
        set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
        set(gca, 'ytick', [0:0.1:1])
        set(gca, 'yticklabels', {})
        if Nind == 1 & cond == 1
            text(-27, 1.17, 'B', 'FontName', 'Helvetica', 'FontSize', fontsz)
        end
        if cond == 1
            set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
            ylabel('Hit rate','FontName', 'Helvetica', 'FontSize', fontsz)
        end
        set(gca, 'tickdir', 'out')
        set(gca,'ticklength', ticklengthsz)
        %xlabel('T-D mean (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
    end
    
    % 6,2,4 FA T-D MEAN
    for Nind = 1:4
        binzz_mean_mean = [mean(squeeze(binzz_meanE(:,cond,Nind,:)),1)]';
        bincenterz_mean = [(binzz_mean_mean(1:end-1)+binzz_mean_mean(2:end))/2]';
        
        tight_subplot(14,2,[8 9],0.5+cond/2, guttera3, marginsa2)
        
        if model_pred
            hacc(Nind)=fill([bincenterz_mean bincenterz_mean(end:-1:1)], [[squeeze(mean(p_pred_ta_mean(:,cond,Nind,:),1))]'-[squeeze(std(p_pred_ta_mean(:,cond,Nind,:),1))./sqrt(squeeze(sum(~isnan(p_pred_ta_mean(:,cond,Nind,:)),1)))]'...
                [squeeze(mean(p_pred_ta_mean(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(p_pred_ta_mean(:,cond,Nind,end:-1:1),1))./sqrt(squeeze(sum(~isnan(p_pred_ta_mean(:,cond,Nind,:)),1)))]'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
            set(hacc(Nind), 'FaceAlpha', face_alpha_val)
        end
    end
    for Nind = 1:4
        binzz_mean_mean = [mean(squeeze(binzz_meanE(:,cond,Nind,:)),1)]';
        bincenterz_mean = [(binzz_mean_mean(1:end-1)+binzz_mean_mean(2:end))/2]';
        
        tight_subplot(14,2,[8 9],0.5+cond/2, guttera3, marginsa2)
        if model_pred
            plot(bincenterz_mean, squeeze(mean(p_data_ta_mean(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
            leg(Nind)=errorbar(bincenterz_mean,squeeze(mean(p_data_ta_mean(:,cond,Nind,:),1)), squeeze(std(p_data_ta_mean(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none',  'Linewidth', 1.3);
        else
            plot(bincenterz_mean, squeeze(mean(p_data_ta_mean(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
            leg(Nind)=errorbar(bincenterz_mean,squeeze(mean(p_data_ta_mean(:,cond,Nind,:),1)), squeeze(std(p_data_ta_mean(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:),  'Linewidth', 1.3);
        end
        %chance = plot(0:10:90, 0.5*ones(1, 10),'--k'); hold on;
        box off
        ylim([0.0 1])
        xlim([0 90])
        set(gca, 'xtick', [0:15:90])
        set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
        set(gca, 'ytick', [0:0.1:1])
        set(gca, 'yticklabels', {})
        if cond == 1
            set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
            ylabel('False-alarm rate','FontName', 'Helvetica', 'FontSize', fontsz)
        end
        set(gca, 'tickdir', 'out')
        set(gca,'ticklength', ticklengthsz)
        xlabel('T-D mean (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
    end
    
    binzz_meanE =  binzz_meanE * pi/90;
    
    % 6,2,5 HR CVAR
    for Nind = 2:4
        binzz_cvar_mean = [mean(squeeze(binzz_cvarE(:,cond,Nind,:)),1)]';
        bincenterz_cvar = [(binzz_cvar_mean(1:end-1)+binzz_cvar_mean(2:end))/2]';
        
        tight_subplot(14,2,[ 11 12],0.5+cond/2, guttera3, marginsa2)
        
        if model_pred
            hacc(Nind)=fill([bincenterz_cvar bincenterz_cvar(end:-1:1)], [[squeeze(mean(p_pred_tp_cvar(:,cond,Nind,:),1))]'-[squeeze(std(p_pred_tp_cvar(:,cond,Nind,:),1))./sqrt(squeeze(sum(~isnan(p_pred_tp_cvar(:,cond,Nind,:)),1)))]'...
                [squeeze(mean(p_pred_tp_cvar(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(p_pred_tp_cvar(:,cond,Nind,end:-1:1),1))./sqrt(squeeze(sum(~isnan(p_pred_tp_cvar(:,cond,Nind,:)),1)))]'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
            set(hacc(Nind), 'FaceAlpha', face_alpha_val)
        end
    end
    for Nind = 2:4
        
        binzz_cvar_mean = [mean(squeeze(binzz_cvarE(:,cond,Nind,:)),1)]';
        bincenterz_cvar = [(binzz_cvar_mean(1:end-1)+binzz_cvar_mean(2:end))/2]';
        
        tight_subplot(14,2,[11 12],0.5+cond/2, guttera3, marginsa2)
        if model_pred
            plot(bincenterz_cvar, squeeze(mean(p_data_tp_cvar(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
            leg(Nind)=errorbar(bincenterz_cvar,squeeze(mean(p_data_tp_cvar(:,cond,Nind,:),1)), squeeze(std(p_data_tp_cvar(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none',  'Linewidth', 1.3);
        else
            plot(bincenterz_cvar, squeeze(mean(p_data_tp_cvar(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
            leg(Nind)=errorbar(bincenterz_cvar,squeeze(mean(p_data_tp_cvar(:,cond,Nind,:),1)), squeeze(std(p_data_tp_cvar(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:),  'Linewidth', 1.3);
        end
        %chance = plot(0:10:90, 0.5*ones(1, 10),'--k'); hold on;
        box off
        ylim([0.3 1])
        xlim([0 1])
        set(gca, 'xtick',[0:0.2:1])
        set(gca, 'ytick', [0:0.1:1])
        set(gca, 'tickdir', 'out')
        set(gca, 'yticklabels', {})
        if cond == 1
            set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
            ylabel('Hit rate','FontName', 'Helvetica', 'FontSize', fontsz)
            text(-0.31, 1.17, 'C', 'FontName', 'Helvetica', 'FontSize', fontsz)
        end
        set(gca, 'tickdir', 'out')
        set(gca,'ticklength', ticklengthsz)
        %xlabel('distractor variance (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
    end
    
    
    % 6,2,6 FA CVAR
    for Nind = 2:4%4:-1:2%2:4
        binzz_cvar_mean = [mean(squeeze(binzz_cvarE(:,cond,Nind,:)),1)]';
        bincenterz_cvar = [(binzz_cvar_mean(1:end-1)+binzz_cvar_mean(2:end))/2]';
        
        tight_subplot(14,2,[13 14],0.5+cond/2, guttera3, marginsa2)
        
        if model_pred
            hacc(Nind)=fill([bincenterz_cvar bincenterz_cvar(end:-1:1)], [[squeeze(mean(p_pred_ta_cvar(:,cond,Nind,:),1))]'-[squeeze(std(p_pred_ta_mean(:,cond,Nind,:),1))./sqrt(squeeze(sum(~isnan(p_pred_ta_cvar(:,cond,Nind,:)),1)))]'...
                [squeeze(mean(p_pred_ta_cvar(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(p_pred_ta_cvar(:,cond,Nind,end:-1:1),1))./sqrt(squeeze(sum(~isnan(p_pred_ta_cvar(:,cond,Nind,:)),1)))]'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
            set(hacc(Nind), 'FaceAlpha', face_alpha_val)
        end
    end
    for Nind = 2:4 %2:4
        binzz_cvar_mean = [mean(squeeze(binzz_cvarE(:,cond,Nind,:)),1)]';
        bincenterz_cvar = [(binzz_cvar_mean(1:end-1)+binzz_cvar_mean(2:end))/2]';
        
        tight_subplot(14,2,[13 14],0.5+cond/2, guttera3, marginsa2)
        if model_pred
            plot(bincenterz_cvar, squeeze(mean(p_data_ta_cvar(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
            leg(Nind)=errorbar(bincenterz_cvar,squeeze(mean(p_data_ta_cvar(:,cond,Nind,:),1)), squeeze(std(p_data_ta_cvar(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none',  'Linewidth', 1.3);
        else
            plot(bincenterz_cvar, squeeze(mean(p_data_ta_cvar(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
            leg(Nind)=errorbar(bincenterz_cvar,squeeze(mean(p_data_ta_cvar(:,cond,Nind,:),1)), squeeze(std(p_data_ta_cvar(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:),  'Linewidth', 1.3);
        end
        %chance = plot(0:10:90, 0.5*ones(1, 10),'--k'); hold on;
        box off
        ylim([0.0 1])
        xlim([0 1])
        %set(gca, 'xtick', [0:15:90])
        %set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
        set(gca, 'ytick', [0:0.1:1])
        set(gca, 'yticklabels', {})
        set(gca, 'xtick',[0:0.2:1])
        set(gca, 'ytick', [0:0.1:1])
        set(gca, 'yticklabels', {})
        set(gca, 'tickdir', 'out')
        set(gca,'ticklength', ticklengthsz)
        xlabel('Distractor variance','FontName', 'Helvetica', 'FontSize', fontsz)
        
        if cond == 1
            set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
            ylabel('False-alarm rate','FontName', 'Helvetica', 'FontSize', fontsz)
        end
        set(gca, 'tickdir', 'out')
        set(gca,'ticklength', ticklengthsz)
        xlabel('distractor variance (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
    end
end
if save_flag
    psname = ['HR_FA_DET_exp_', num2str(exp_i), '_model_', num2str(mi), '.pdf']
    print_pdf(psname)
end

%% proportion correct with target value

%close all;
figure(4)
set(gcf, 'Position', [100 100 520 280])
set_plotting_params();
colorsz_mean =  [255 165 0; 154 187 246;  204 221 123 ]./255;%[255 165 0; 154 217 246;  0 0 0 ]./255; %[154 217 246;  0 0 0; 255 165 0]./255;

guttera3 = [0.07 0.03]; %[0.06 0.04];
marginsa3 = [0.12 0.14 0.22 0.12]; %marginsa; %left right bottom top


bincenterz_target_val = bincenterz_target_val * 180/pi;

model_pred = 1;
for cond = [1 3]
    
    for Nind = 1:4
        
        tight_subplot(1,2,[ 1],0.5+cond/2, guttera3, marginsa3)
        
        if model_pred
            hacc(Nind)=fill([bincenterz_target_val bincenterz_target_val(end:-1:1)], [[squeeze(mean(p_pred_tval_pc(:,cond,Nind,:),1))]'-[squeeze(std(p_pred_tval_pc(:,cond,Nind,:),1))./sqrt(squeeze(sum(~isnan(p_pred_tval_pc(:,cond,Nind,:)),1)))]'...
                [squeeze(mean(p_pred_tval_pc(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(p_pred_tval_pc(:,cond,Nind,end:-1:1),1))./sqrt(squeeze(sum(~isnan(p_pred_tval_pc(:,cond,Nind,:)),1)))]'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
            set(hacc(Nind), 'FaceAlpha', face_alpha_val)
        end
    end
    for Nind = 1:4
        
        
        tight_subplot(1,2,[1],0.5+cond/2, guttera3, marginsa3)
        %if model_pred
        %plot(bincenterz_target_val, squeeze(mean(p_data_tval_pc(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
        %leg(Nind)=errorbar(bincenterz_target_val,squeeze(mean(p_data_tval_pc(:,cond,Nind,:),1)), squeeze(std(p_data_tval_pc(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none',  'Linewidth', 1.3);
        %else
        plot(bincenterz_target_val, squeeze(mean(p_data_tval_pc(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
        leg(Nind)=errorbar(bincenterz_target_val,squeeze(mean(p_data_tval_pc(:,cond,Nind,:),1)), squeeze(std(p_data_tval_pc(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:),  'Linewidth', 1.3);
        %end
        %chance = plot(0:10:90, 0.5*ones(1, 10),'--k'); hold on;
        box off
        ylim([0.3 1])
        xlim([-20 160])
        %set(gca, 'xtick', [-180:30:180])
        %set(gca, 'xticklabels', {'-180', '', '-120', '', '-60', '', '0','', '60','', '120', '', '180' },'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
        set(gca, 'xtick', [-22.5:22.5:(135+22.5)])
        set(gca, 'xticklabels', {'-22.5', '0', '', '45', '', '90', '', '135',''},'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
        set(gca, 'ytick', [0:0.1:1])
        set(gca, 'yticklabels', {})
        if cond == 1
            set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
            ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
        end
        set(gca, 'tickdir', 'out')
        set(gca,'ticklength', ticklengthsz)
        xlabel('target value (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
        if Nind == 1 & cond == 1
            text(22.5,1.1,'Perception', 'FontName', 'Helvetica', 'FontSize', fontsz)
            %text(-90,1.1,'Perception', 'FontName', 'Helvetica', 'FontSize', fontsz)
            %text(-27, 1.2, 'A', 'FontName', 'Helvetica', 'FontSize', fontsz)
        elseif Nind == 1 &  cond == 3
            text(45,1.1,'Memory', 'FontName', 'Helvetica', 'FontSize', fontsz)
            %text(-90,1.1,'Memory', 'FontName', 'Helvetica', 'FontSize', fontsz)
        end
    end
    if cond == 1
        llpa = legend([leg(1) leg(2) leg(3) leg(4)],'N = 2', 'N = 3','N = 4', 'N = 6','FontName','Helvetica','FontSize', fontsz)
        legend boxoff
        set(llpa, 'FontName','Helvetica','FontSize', fontsz*0.9, 'Position', [0.9 0.84 0.04 0.04])
    end
    
end
if save_flag
    psname = ['Prop_corr_with_target_value_DET_NEWE_nbinzE_', num2str(length(bincenterz_target_val)), '_exp_', num2str(exp_i), '_model_', num2str(mi), '.pdf']
    print_pdf(psname)
end
%% 2-way repeated-measures ANOVA: set size and target orientation

%% Perception - cond = 1
cond = 1;
%squeeze(p_data_tval_pc(:,cond,Nind,:))
vals = [squeeze(p_data_tval_pc(:,cond,1,1)); squeeze(p_data_tval_pc(:,cond,1,2));  squeeze(p_data_tval_pc(:,cond,1,3)); squeeze(p_data_tval_pc(:,cond,1,4)) ;
    squeeze(p_data_tval_pc(:,cond,2,1)); squeeze(p_data_tval_pc(:,cond,2,2));  squeeze(p_data_tval_pc(:,cond,2,3)); squeeze(p_data_tval_pc(:,cond,2,4));
    squeeze(p_data_tval_pc(:,cond,3,1)); squeeze(p_data_tval_pc(:,cond,3,2));  squeeze(p_data_tval_pc(:,cond,3,3)); squeeze(p_data_tval_pc(:,cond,3,4));
    squeeze(p_data_tval_pc(:,cond,4,1)); squeeze(p_data_tval_pc(:,cond,4,2));  squeeze(p_data_tval_pc(:,cond,4,3)); squeeze(p_data_tval_pc(:,cond,4,4))];
% 4 orientation values, 4 set sizes
%Nsubj = 11;

Subjects = [[1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj]...
    [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj]];
group = [[ones(1,Nsubj) 2*ones(1,Nsubj)  3*ones(1,Nsubj) 4*ones(1,Nsubj)...
    ones(1,Nsubj) 2*ones(1,Nsubj)  3*ones(1,Nsubj) 4*ones(1,Nsubj)...
    ones(1,Nsubj) 2*ones(1,Nsubj)  3*ones(1,Nsubj) 4*ones(1,Nsubj)...
    ones(1,Nsubj) 2*ones(1,Nsubj)  3*ones(1,Nsubj) 4*ones(1,Nsubj)]'...
    [ones(1,(4*Nsubj)) 2*ones(1,(4*Nsubj)) 3*ones(1,(4*Nsubj)) 4*ones(1,(4*Nsubj))]'];

y = vals;
g1 = group(:,1);
g2 = group(:,2);


[p_full,table_full,stats_full] = anovan(y,{g1,g2, Subjects}, 'random',3,'varnames', {'target ori val', 'Set size', 'Subject'}, 'model', 'full');

%
cond = 3;
%squeeze(p_data_tval_pc(:,cond,Nind,:))
vals = [squeeze(p_data_tval_pc(:,cond,1,1)); squeeze(p_data_tval_pc(:,cond,1,2));  squeeze(p_data_tval_pc(:,cond,1,3)); squeeze(p_data_tval_pc(:,cond,1,4)) ;
    squeeze(p_data_tval_pc(:,cond,2,1)); squeeze(p_data_tval_pc(:,cond,2,2));  squeeze(p_data_tval_pc(:,cond,2,3)); squeeze(p_data_tval_pc(:,cond,2,4));
    squeeze(p_data_tval_pc(:,cond,3,1)); squeeze(p_data_tval_pc(:,cond,3,2));  squeeze(p_data_tval_pc(:,cond,3,3)); squeeze(p_data_tval_pc(:,cond,3,4));
    squeeze(p_data_tval_pc(:,cond,4,1)); squeeze(p_data_tval_pc(:,cond,4,2));  squeeze(p_data_tval_pc(:,cond,4,3)); squeeze(p_data_tval_pc(:,cond,4,4))];
% 4 orientation values, 4 set sizes
%Nsubj = 11;

Subjects = [[1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj]...
    [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj]];
group = [[ones(1,Nsubj) 2*ones(1,Nsubj)  3*ones(1,Nsubj) 4*ones(1,Nsubj)...
    ones(1,Nsubj) 2*ones(1,Nsubj)  3*ones(1,Nsubj) 4*ones(1,Nsubj)...
    ones(1,Nsubj) 2*ones(1,Nsubj)  3*ones(1,Nsubj) 4*ones(1,Nsubj)...
    ones(1,Nsubj) 2*ones(1,Nsubj)  3*ones(1,Nsubj) 4*ones(1,Nsubj)]'...
    [ones(1,(4*Nsubj)) 2*ones(1,(4*Nsubj)) 3*ones(1,(4*Nsubj)) 4*ones(1,(4*Nsubj))]'];

y = vals;
g1 = group(:,1);
g2 = group(:,2);


[p_full,table_full,stats_full] = anovan(y,{g1,g2, Subjects}, 'random',3,'varnames', {'target ori val', 'Set size', 'Subject'}, 'model', 'full');

%% group performance for 0 and 90, and respectively 45 and 135

cond = 1;
%squeeze(p_data_tval_pc(:,cond,Nind,:))
vals = [mean([squeeze(p_data_tval_pc(:,cond,1,1))';   squeeze(p_data_tval_pc(:,cond,1,3))'])'; mean([ squeeze(p_data_tval_pc(:,cond,1,2))'; squeeze(p_data_tval_pc(:,cond,1,4))'])' ;
    mean([squeeze(p_data_tval_pc(:,cond,2,1))';   squeeze(p_data_tval_pc(:,cond,2,3))'])'; mean([ squeeze(p_data_tval_pc(:,cond,2,2))'; squeeze(p_data_tval_pc(:,cond,2,4))'])' ;
    mean([squeeze(p_data_tval_pc(:,cond,3,1))';   squeeze(p_data_tval_pc(:,cond,3,3))'])'; mean([ squeeze(p_data_tval_pc(:,cond,3,2))'; squeeze(p_data_tval_pc(:,cond,3,4))'])' ;
    mean([squeeze(p_data_tval_pc(:,cond,4,1))';   squeeze(p_data_tval_pc(:,cond,4,3))'])'; mean([ squeeze(p_data_tval_pc(:,cond,4,2))'; squeeze(p_data_tval_pc(:,cond,4,4))'])' ;
    ];
% 4 orientation values, 4 set sizes
%Nsubj = 11;

Subjects = [[1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj]...
    ];
group = [[ones(1,Nsubj) 2*ones(1,Nsubj)  ...
    ones(1,Nsubj) 2*ones(1,Nsubj)  ...
    ones(1,Nsubj) 2*ones(1,Nsubj)  ...
    ones(1,Nsubj) 2*ones(1,Nsubj)  ]'...
    [ones(1,(2*Nsubj)) 2*ones(1,(2*Nsubj)) 3*ones(1,(2*Nsubj)) 4*ones(1,(2*Nsubj))]'];

y = vals;
g1 = group(:,1);
g2 = group(:,2);


[p_full,table_full,stats_full] = anovan(y,{g1,g2, Subjects}, 'random',3,'varnames', {'target ori val', 'Set size', 'Subject'}, 'model', 'full');

%
cond = 3;
%squeeze(p_data_tval_pc(:,cond,Nind,:))
vals = [mean([squeeze(p_data_tval_pc(:,cond,1,1))';   squeeze(p_data_tval_pc(:,cond,1,3))'])'; mean([ squeeze(p_data_tval_pc(:,cond,1,2))'; squeeze(p_data_tval_pc(:,cond,1,4))'])' ;
    mean([squeeze(p_data_tval_pc(:,cond,2,1))';   squeeze(p_data_tval_pc(:,cond,2,3))'])'; mean([ squeeze(p_data_tval_pc(:,cond,2,2))'; squeeze(p_data_tval_pc(:,cond,2,4))'])' ;
    mean([squeeze(p_data_tval_pc(:,cond,3,1))';   squeeze(p_data_tval_pc(:,cond,3,3))'])'; mean([ squeeze(p_data_tval_pc(:,cond,3,2))'; squeeze(p_data_tval_pc(:,cond,3,4))'])' ;
    mean([squeeze(p_data_tval_pc(:,cond,4,1))';   squeeze(p_data_tval_pc(:,cond,4,3))'])'; mean([ squeeze(p_data_tval_pc(:,cond,4,2))'; squeeze(p_data_tval_pc(:,cond,4,4))'])' ;
    ];
% 4 orientation values, 4 set sizes
%Nsubj = 11;

Subjects = [[1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj] [1:Nsubj]...
    ];
group = [[ones(1,Nsubj) 2*ones(1,Nsubj)  ...
    ones(1,Nsubj) 2*ones(1,Nsubj)  ...
    ones(1,Nsubj) 2*ones(1,Nsubj)  ...
    ones(1,Nsubj) 2*ones(1,Nsubj)  ]'...
    [ones(1,(2*Nsubj)) 2*ones(1,(2*Nsubj)) 3*ones(1,(2*Nsubj)) 4*ones(1,(2*Nsubj))]'];

y = vals;
g1 = group(:,1);
g2 = group(:,2);


[p_full,table_full,stats_full] = anovan(y,{g1,g2, Subjects}, 'random',3,'varnames', {'target ori val', 'Set size', 'Subject'}, 'model', 'full');




%% linear separability
% for this, we also need to load the type=2
%load(['all_vars_exp_',num2str(exp_i),'_type_',num2str(2) ,'_mi_',num2str(mi), '.mat'])
load(['lin_sep_pred_mi_',num2str(mi),'_exp_',num2str(exp_i),'_type_2.mat'])

acc_lin_sep_pred(:,2,:) = squeeze(acc_lin_sep_pred_type_2(:,2,:));
acc_lin_sep_pred(:,4,:) = squeeze(acc_lin_sep_pred_type_2(:,4,:));
acc_non_lin_sep_pred(:,2,:) = squeeze(acc_non_lin_sep_pred_type_2(:,2,:));
acc_non_lin_sep_pred(:,4,:) = squeeze(acc_non_lin_sep_pred_type_2(:,4,:));

%close all;
figure(5)
greyy = [128 128 0]/255;%[128 128 128]/255; % really olive
greyy_shade = (2*greyy*255+ [255 255 255])/(3*255);
blackk = [202 146 56]/255;%[0 0 0]/255;  % really orange
blackk_shade = (2*blackk*255+ [255 255 255])/(3*255);

sz_dot2 = 3*sz_dot;
for cond = 1:4%[1 3]
    
    % a)   set size
    %tight_subplot(5,2,1,0.5+cond/2, guttera2, marginsa2)
    if cond == 1 | cond == 3
        tight_subplot(4,2,1,0.5+cond/2, guttera2, marginsa2)
        chance=plot(Nvec, 0.5*ones(1, length(Nvec)),'--k'); hold on;
        ylim([0.49 0.8])%ylim([0.4 0.8])
        
    else
        tight_subplot(4,2,2,cond/2, guttera2, marginsa2)
        chance=plot(Nvec, 1./Nvec,'--k'); hold on;
        ylim([0.49 0.8])%ylim([0.3 0.8])
        
    end
    
    plot(Nvec(2:3), squeeze(mean(acc_non_lin_sep(:, cond, 2:3),1)), 'Color', greyy); hold on;
    lsep_no = errorbar(Nvec(2:3), squeeze(mean(acc_non_lin_sep(:, cond, 2:3),1)),squeeze(std(acc_non_lin_sep(:, cond, 2:3),1))/sqrt(Nsubj), 'o','MarkerSize', sz_dot2, 'MarkerFaceColor', greyy, 'MarkerEdgeColor', greyy,'Color',greyy , 'LineWidth',1.3, 'Capsize', 0); hold on;
    plot(Nvec(2:3), squeeze(mean(acc_non_lin_sep_pred(:, cond, 2:3),1)), 'Color', greyy_shade, 'Linewidth', 0.01); hold on;
    
    
    plot(Nvec(2:3), squeeze(mean(acc_lin_sep(:, cond, 2:3),1)), 'Color', blackk); hold on;
    lsep_yes = errorbar(Nvec(2:3), squeeze(mean(acc_lin_sep(:, cond, 2:3),1)),squeeze(std(acc_lin_sep(:, cond, 2:3),1))/sqrt(Nsubj), 'o','MarkerSize', sz_dot2, 'MarkerFaceColor', blackk, 'MarkerEdgeColor',blackk,'Color',blackk, 'LineWidth',1.3,'Capsize', 0); hold on;
    plot(Nvec(2:3), squeeze(mean(acc_lin_sep_pred(:, cond, 2:3),1)), 'Color', blackk_shade); hold on;
    
    if model_pred
        h_non_sep = fill([Nvec(2:3) Nvec(3:-1:2)], [[squeeze(mean(acc_non_lin_sep_pred(:,cond,2:3),1))]'-[squeeze(std(acc_non_lin_sep_pred(:,cond,2:3),1))./sqrt(squeeze(sum(~isnan(acc_non_lin_sep_pred(:,cond,2:3)),1)))]'...
            [squeeze(mean(acc_non_lin_sep_pred(:,cond,3:-1:2),1))]'+[squeeze(std(acc_non_lin_sep_pred(:,cond,3:-1:2),1))./sqrt(squeeze(sum(~isnan(acc_non_lin_sep_pred(:,cond,3:-1:2,:)),1)))]'],greyy_shade,'EdgeColor', 'None'); hold on;
        set(h_non_sep, 'FaceAlpha', face_alpha_val/1.5)
        
        h_sep = fill([Nvec(2:3) Nvec(3:-1:2)], [[squeeze(mean(acc_lin_sep_pred(:,cond,2:3),1))]'-[squeeze(std(acc_lin_sep_pred(:,cond,2:3),1))./sqrt(squeeze(sum(~isnan(acc_lin_sep_pred(:,cond,2:3)),1)))]'...
            [squeeze(mean(acc_lin_sep_pred(:,cond,3:-1:2),1))]'+[squeeze(std(acc_lin_sep_pred(:,cond,3:-1:2),1))./sqrt(squeeze(sum(~isnan(acc_lin_sep_pred(:,cond,3:-1:2,:)),1)))]'],blackk_shade,'EdgeColor', 'None'); hold on;
        set(h_sep, 'FaceAlpha', face_alpha_val/1.5)
    end
    
    
    if cond == 1
        lch =legend([chance], 'chance', 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz )
        legend boxoff
        set(lch, 'Location', 'Northeast')
    elseif cond == 2
        lllrev =legend([ lsep_yes lsep_no], 'Lin sep trials', 'Non lin sep trials','FontName', 'Helvetica', 'FontSize', 0.9*fontsz )
        legend boxoff
        %set(lllrev, 'Location', 'Southeast') % 0.1917    0.5393    0.2295    0.0774
        set(lllrev, 'Position', [0.19 0.33 0.22 0.07])
    end
    
    xlim([min(Nvec)-0.5 max(Nvec)+0.5]);
    box off
    set(gca, 'tickdir', 'out')
    set(gca, 'xtick', Nvec)
    set(gca, 'xticklabels', {'2', '3', '4', '6'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
    set(gca, 'ytick', [0:0.1:1])
    %set(gca, 'yticklabels', {})
    set(gca,'ticklength', ticklengthsz)
    xlabel('Set size','FontName', 'Helvetica', 'FontSize', fontsz)
    if cond == 1 | cond == 2
        ylabel('Accuracy','FontName', 'Helvetica', 'FontSize', fontsz)
    end
    set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
    if cond == 1
        %text(3, 1.16, text_labels{(1+cond)/2},  'FontName', 'Helvetica', 'FontSize', fontsz)
        text(3, 0.86, text_labels{(1+cond)/2},  'FontName', 'Helvetica', 'FontSize', fontsz)
    elseif cond == 3
        text(3, 0.86, text_labels{(1+cond)/2},  'FontName', 'Helvetica', 'FontSize', fontsz)
    end
    
    
end


if save_flag
    psname = ['Linear_sepEEE_mi_',num2str(mi), '_exp_', num2str(exp_i), '.pdf']
    print_pdf(psname)
end


