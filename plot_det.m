function plot_det(exp_i,type,mi,model_pred,full_figure_flag, save_flag)
% Detection % if type ==1 or 3

load(['all_vars_exp_',num2str(exp_i),'_type_',num2str(type) ,'_mi_',num2str(mi), '.mat'])

set_plotting_params();

figure
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
