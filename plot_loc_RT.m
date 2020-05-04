function plot_loc_RT(exp_i,type,mi,model_pred, model_predR,full_figure_flag, save_flag)
% Localization % if type ==2 or 3

%load('all_vars.mat')
load(['all_vars_exp_',num2str(exp_i),'_type_',num2str(type) ,'_mi_',num2str(mi), '.mat'])

set_plotting_params();


figure
if full_figure_flag
    set(gcf, 'Position', [100 100 500 920])
    guttera2 = guttera;
    marginsa2 = marginsa;
else
    set(gcf, 'Position', [100 100 500 190])
    guttera2 = guttera;
    marginsa2 = [0.14 0.14 0.2 0.2];
end



%ylim_min = 0.2 ; ylim_max = 0.7;
if exp_i == 1
    ylim_min = 0.4 ; ylim_max = 0.9;
elseif exp_i == 2
    ylim_min = 0.85 ; ylim_max = 1.4;
end

for cond = [2 4]
    
    if full_figure_flag
        %a) with set size
        tight_subplot(5,2,1,cond/2, guttera2, marginsa2)
        
        
        plot(Nvec, squeeze(mean(rt(:, cond, :),1)), 'Color', 'k', 'Linewidth', 1.3); hold on;
        errorbar(Nvec, squeeze(mean(rt(:, cond, :),1)),squeeze(std(rt(:, cond, :),1))/sqrt(Nsubj), 'o','MarkerSize', sz_dot, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color', 'k', 'LineWidth',1.3); hold on;
        
        box off
        ylim([ylim_min ylim_max])
        xlim([min(Nvec)-0.5 max(Nvec)+0.5]);
        box off
        set(gca, 'tickdir', 'out')
        set(gca, 'xtick', Nvec)
        set(gca, 'xticklabels', {'2','3', '4', '6'},'FontName', 'Helvetica', 'FontSize', fontsz)
        set(gca, 'ytick', [0:0.1:ylim_max])
        set(gca, 'yticklabels',{})
        set(gca, 'tickdir', 'out')
        set(gca,'ticklength', [0.04 0.04])
        xlabel('Set size','FontName', 'Helvetica', 'FontSize', fontsz)
        
        if cond == 2
            set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1','','1.2','','1.4','','1.6'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
            ylabel('Reaction time (sec)','FontName', 'Helvetica', 'FontSize', fontsz)
        end
        text(3, ylim_max*1.1, text_labels{cond/2},  'FontName', 'Helvetica', 'FontSize', fontsz)
        
        
        %b) T- MSD
        for Nind = 1:4
            bincentersl=[squeeze(bincenterz(cond, Nind,1,:))]';
            
            
            tight_subplot(5,2,2,cond/2, guttera2, marginsa2)
            
            
            plot(bincentersl, squeeze(mean(rt_data_w_msd(:,cond,Nind,:),1)), 'o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
            leg(Nind) = errorbar(bincentersl,squeeze(mean(rt_data_w_msd(:,cond,Nind,:),1)), squeeze(std(rt_data_w_msd(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linewidth', 1.3);
            
            
            box off
            ylim([ylim_min ylim_max])
            xlim([0 90])
            set(gca, 'xtick', [0:15:90])
            set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', fontsz)
            set(gca, 'ytick', [0:0.1:ylim_max])
            set(gca, 'yticklabels', {})
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', [0.04 0.04])
            if cond == 2
                set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1','','1.2','','1.4','','1.6'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
                ylabel('Reaction time (sec)','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            xlabel('min T-D difference (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
            
        end
        
        if cond == 4
            llpa = legend([leg(1) leg(2) leg(3) leg(4)],'N = 2', 'N = 3','N = 4', 'N = 6','FontName','Helvetica','FontSize', fontsz)
            legend boxoff
            set(llpa, 'FontName','Helvetica','FontSize', fontsz*0.9, 'Position', [0.9 0.84 0.04 0.04])
        end
        
        
        
        %c) T-mean of distractors
        binzz_meanE = 90/pi* binzz_meanE;  %since we doubled all our stimuli for modelling on the (-pi,pi) full space
        
        for Nind = 2:4
            
            binzz_mean_mean = [mean(squeeze(binzz_meanE(:,cond,Nind,:)),1)]';
            bincenterz_mean = [(binzz_mean_mean(1:end-1)+binzz_mean_mean(2:end))/2]';
            
            tight_subplot(5,2,3,cond/2, guttera2, marginsa2)
            
            plot(bincenterz_mean, squeeze(mean(rt_data_w_mean(:,cond,Nind,:),1)), 'o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
            leg(Nind) = errorbar(bincenterz_mean,squeeze(mean(rt_data_w_mean(:,cond,Nind,:),1)), squeeze(std(rt_data_w_mean(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linewidth', 1.3);
            
            
            box off
            ylim([ylim_min ylim_max])
            xlim([0 90])
            set(gca, 'xtick', [0:15:90])
            set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', fontsz)
            set(gca, 'ytick', [0:0.1:ylim_max])
            set(gca, 'yticklabels', {})
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', [0.04 0.04])
            
            xlabel('T-D mean (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
            
            if cond == 2
                set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1','','1.2','','1.4','','1.6'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
                ylabel('Reaction time (sec)','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            
        end
        
        
        binzz_meanE = pi/90* binzz_meanE;
        
        
        
        % d) circ var of distractors
        for Nind = 2:4
            
            binzz_cvar_mean = [mean(squeeze(binzz_cvarE(:,cond,Nind,:)),1)]';
            
            bincenterz_cvar = [(binzz_cvar_mean(1:end-1)+binzz_cvar_mean(2:end))/2]';
            
            tight_subplot(5,2,4,cond/2, guttera2, marginsa2)
            
            
            plot(bincenterz_cvar, squeeze(nanmean(rt_data_w_cvar(:,cond,Nind,:),1)), 'o','Color', colors_sz(Nind,:), 'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
            leg(Nind) = errorbar(bincenterz_cvar,squeeze(nanmean(rt_data_w_cvar(:,cond,Nind,:),1)), squeeze(nanstd(rt_data_w_cvar(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linewidth', 1.3);
            
            
            box off
            ylim([ylim_min ylim_max])
            xlim([0 1])
            set(gca, 'xtick',[0:0.2:1])
            set(gca, 'ytick', [0:0.1:ylim_max])
            set(gca, 'yticklabels', {})
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', [0.04 0.04])
            
            
            xlabel('Distractor variance','FontName', 'Helvetica', 'FontSize', fontsz)
            
            if cond == 2
                set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1','','1.2','','1.4','','1.6'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
                ylabel('Reaction time (sec)','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            
        end
        if cond == 4
            llpa = legend([leg(1) leg(2) leg(3) leg(4)],'N = 2', 'N = 3','N = 4', 'N = 6','FontName','Helvetica','FontSize', fontsz)
            legend boxoff
            set(llpa, 'FontName','Helvetica','FontSize', fontsz*0.9, 'Position', [0.9 0.84 0.04 0.04])
        end
        
        
        
    else
        
        %b) T- MSD
        for Nind = 1:4
            bincentersl=[squeeze(bincenterz(cond, Nind,1,:))]';
            
            
            tight_subplot(1,2,1,cond/2, guttera2, marginsa2)
            
            
            plot(bincentersl, squeeze(mean(rt_data_w_msd(:,cond,Nind,:),1)), 'o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
            leg(Nind) = errorbar(bincentersl,squeeze(mean(rt_data_w_msd(:,cond,Nind,:),1)), squeeze(std(rt_data_w_msd(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:));
            
            
            box off
            ylim([ylim_min ylim_max])
            xlim([0 90])
            set(gca, 'xtick', [0:15:90])
            set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', fontsz)
            set(gca, 'ytick', [0:0.1:ylim_max])
            set(gca, 'yticklabels', {})
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', [0.04 0.04])
            if cond == 2
                set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1','','1.2','','1.4','','1.6'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
                ylabel('Reaction time (sec)','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            xlabel('T - MSD orientation distance (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
            
        end
        if Nind == 1
            text(30, ylim_max*1.1, text_labels{(cond)/2},  'FontName', 'Helvetica', 'FontSize', fontsz)
        end
        if cond == 2
            text(-40, 1.23, 'Localization',  'FontName', 'Helvetica', 'FontSize', fontsz*1.2)
        end
        if cond == 4
            llpa = legend([leg(1) leg(2) leg(3) leg(4)],'N = 2', 'N = 3','N = 4', 'N = 6','FontName','Helvetica','FontSize', fontsz)
            legend boxoff
            set(llpa, 'FontName','Helvetica','FontSize', fontsz*0.9, 'Position', [0.9 0.84 0.04 0.04])
        end
        
        
        
        
    end
    
    psname = ['RT_LOC_full_figure_',num2str(full_figure_flag),'_type_',num2str(type) ,'_exp_',num2str(exp_i), '.pdf']
    if save_flag
        print_pdf(psname)
    end
    
    
end