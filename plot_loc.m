function plot_loc(exp_i,type,mi,model_pred, model_predR,full_figure_flag, save_flag)
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


if model_predR
    p_resp_ori_diff_rank_simim = squeeze(nanmean(p_resp_ori_diff_rank_simi,5));
    %p_resp_loc_diff_rank_simim = squeeze(mean(p_resp_loc_diff_rank_simi,5));
end

% construct chance resp distrib
chance_ori = nan(4,nbinz_jj);
chance_loc = nan(4,nbinz_jj);
for Nind = 2:4
    
    chance_ori(Nind,1:loc_vec(Nind)) = repmat(1/(Nvec(Nind)-1), 1, loc_vec(Nind));
    
end

if exp_i == 1
    chance_loc(2,1:loc_vec(2)) = [2/3 1/3];
    chance_loc(3,1:loc_vec(3)) = [3/6 2/6 1/6];
    chance_loc(4,1:loc_vec(4)) = [2/5 2/5 1/5];
elseif exp_i == 2
    chance_loc(2,1:loc_vec(2)) = [2/3 1/3];
    chance_loc(3,1:loc_vec(3)) = [3/6 2/6 1/6];
    chance_loc(4,1:loc_vec(4)) = [5/15 4/15 3/15 2/15 1/15];
end



if exp_i == 1
    p_resp_max = 0.35;
    x_mem_label = -0.6;
elseif exp_i == 2
    p_resp_max = 0.32;
    x_mem_label = -2.2;
end


ticklengthsz = [0.025 0.025];

for cond = [2 4]
    
    if full_figure_flag
        %a) with set size
        tight_subplot(5,2,1,cond/2, guttera2, marginsa2)
        if model_pred
            h_p=fill([Nvec Nvec(end:-1:1)], [[squeeze(mean(p_pr_c(:,cond,:),1))]'-[squeeze(std(p_pr_c(:,cond,:),1))]'/sqrt(Nsubj)...
                [squeeze(mean(p_pr_c(:,cond,end:-1:1),1))]'+[squeeze(std(p_pr_c(:,cond,:,end:-1:1),1))]'/sqrt(Nsubj)],col_gray,'EdgeColor', 'None'); hold on;
        end
        if model_pred
            plot(Nvec, squeeze(mean(accuracy(:, cond, :),1)), 'Color', 'k', 'Linestyle', 'none'); hold on;
            errorbar(Nvec, squeeze(mean(accuracy(:, cond, :),1)),squeeze(std(accuracy(:, cond, :),1))/sqrt(Nsubj), 'o','MarkerSize', sz_dot, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color', 'k', 'LineWidth',1.3, 'Linestyle', 'none'); hold on;
        else
            plot(Nvec, squeeze(mean(accuracy(:, cond, :),1)), 'Color', 'k'); hold on;
            errorbar(Nvec, squeeze(mean(accuracy(:, cond, :),1)),squeeze(std(accuracy(:, cond, :),1))/sqrt(Nsubj), 'o','MarkerSize', sz_dot, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color', 'k', 'LineWidth',1.3); hold on;
        end
        
        chance = plot(Nvec, 1./Nvec,'--k', 'Linewidth', 1.3); hold on;
        
        box off
        ylim([0 1])
        xlim([min(Nvec)-0.5 max(Nvec)+0.5]);
        box off
        set(gca, 'tickdir', 'out')
        set(gca, 'xtick', Nvec)
        set(gca, 'xticklabels', {'2','3', '4', '6'},'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
        set(gca, 'ytick', [0:0.1:1])
        set(gca, 'yticklabels',{})
        set(gca, 'tickdir', 'out')
        set(gca,'ticklength', ticklengthsz )
        xlabel('Set size','FontName', 'Helvetica', 'FontSize', fontsz)
        
        if cond == 2
            set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
            ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
        end
        text(3, 1.16, text_labels{cond/2},  'FontName', 'Helvetica', 'FontSize', fontsz)
        
        
        %b) T- MSD
        for Nind = 1:4
            bincentersl=[squeeze(bincenterz(cond, Nind,1,:))]';
            
            
            tight_subplot(5,2,2,cond/2, guttera2, marginsa2)
            
            if model_pred
                h_p=fill([bincentersl bincentersl(end:-1:1)], [[squeeze(mean(p_pred_w(:,cond,Nind,:),1))]'-[squeeze(std(p_pred_w(:,cond,Nind,:),1))]'./sqrt(squeeze(sum(~isnan(p_pred_w(:,cond,Nind,:)),1)))'...
                    [squeeze(mean(p_pred_w(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(p_pred_w(:,cond,Nind,end:-1:1),1))]'./sqrt(squeeze(sum(~isnan(p_pred_w(:,cond,Nind,:)),1)))'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
            end
            if model_pred
                plot(bincentersl, squeeze(mean(p_data_w(:,cond,Nind,:),1)), 'o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
                legg(Nind) = errorbar(bincentersl,squeeze(mean(p_data_w(:,cond,Nind,:),1)), squeeze(std(p_data_w(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none', 'Linewidth',1.3);
            else
                plot(bincentersl, squeeze(mean(p_data_w(:,cond,Nind,:),1)), 'o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
                legg(Nind) = errorbar(bincentersl,squeeze(mean(p_data_w(:,cond,Nind,:),1)), squeeze(std(p_data_w(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linewidth',1.3);
            end
            plot(0:10:90, ones(1,10)*1/Nvec(Nind),'--', 'Color',colors_sz(Nind,:), 'Linewidth', 1.3); hold on;
            box off
            ylim([0 1])
            xlim([0 90])
            set(gca, 'xtick', [0:15:90])
            set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
            set(gca, 'ytick', [0:0.1:1])
            set(gca, 'yticklabels', {})
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', ticklengthsz)
            if cond == 2
                set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
                ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            xlabel('min T-D difference (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
            
        end
        
        if cond == 4
            llpa = legend([legg(1) legg(2) legg(3) legg(4) chance],'N = 2', 'N = 3','N = 4', 'N = 6','chance','FontName','Helvetica','FontSize', fontsz)
            legend boxoff
            set(llpa, 'FontName','Helvetica','FontSize', fontsz*0.9, 'Position', [0.9 0.84 0.04 0.04])
        end
        
        
        
        %c) T-mean of distractors
        binzz_meanE = 90/pi* binzz_meanE;  %since we doubled all our stimuli for modelling on the (-pi,pi) full space
        
        for Nind = 2:4
            
            binzz_mean_mean = [mean(squeeze(binzz_meanE(:,cond,Nind,:)),1)]';
            bincenterz_mean = [(binzz_mean_mean(1:end-1)+binzz_mean_mean(2:end))/2]';
            
            tight_subplot(5,2,3,cond/2, guttera2, marginsa2)
            if model_pred
                h_p=fill([bincenterz_mean bincenterz_mean(end:-1:1)], [[squeeze(mean(p_pred_w_mean(:,cond,Nind,:),1))]'-[squeeze(std(p_pred_w_mean(:,cond,Nind,:),1))]'./sqrt(squeeze(sum(~isnan(p_pred_w_mean(:,cond,Nind,:)),1)))'...
                    [squeeze(mean(p_pred_w_mean(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(p_pred_w_mean(:,cond,Nind,end:-1:1),1))]'./sqrt(squeeze(sum(~isnan(p_pred_w_mean(:,cond,Nind,:)),1)))'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
            end
            if model_pred
                plot(bincenterz_mean, squeeze(mean(p_data_w_mean(:,cond,Nind,:),1)), 'o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
                leg(Nind) = errorbar(bincenterz_mean,squeeze(mean(p_data_w_mean(:,cond,Nind,:),1)), squeeze(std(p_data_w_mean(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none', 'Linewidth',1.3);
            else
                plot(bincenterz_mean, squeeze(mean(p_data_w_mean(:,cond,Nind,:),1)), 'o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
                leg(Nind) = errorbar(bincenterz_mean,squeeze(mean(p_data_w_mean(:,cond,Nind,:),1)), squeeze(std(p_data_w_mean(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linewidth',1.3);
            end
            plot(0:10:90, ones(1,10)*1/Nvec(Nind),'--', 'Color',colors_sz(Nind,:), 'Linewidth', 1.3); hold on;
            box off
            ylim([0 1])
            xlim([0 90])
            set(gca, 'xtick', [0:15:90])
            set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
            set(gca, 'ytick', [0:0.1:1])
            set(gca, 'yticklabels', {})
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', ticklengthsz)
            
            xlabel('T-D mean (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
            
            if cond == 2
                set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
                ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            
        end
        
        
        binzz_meanE = pi/90* binzz_meanE;
        
        
        
        % d) circ var of distractors
        for Nind = 2:4
            
            binzz_cvar_mean = [mean(squeeze(binzz_cvarE(:,cond,Nind,:)),1)]';
            
            bincenterz_cvar = [(binzz_cvar_mean(1:end-1)+binzz_cvar_mean(2:end))/2]';
            
            tight_subplot(5,2,4,cond/2, guttera2, marginsa2)
            
            if model_pred
                h_p=fill([bincenterz_cvar bincenterz_cvar(end:-1:1)], [[squeeze(nanmean(p_pred_w_cvar(:,cond,Nind,:),1))]'-[squeeze(nanstd(p_pred_w_cvar(:,cond,Nind,:),1))]'./sqrt(squeeze(sum(~isnan(p_pred_w_cvar(:,cond,Nind,:)),1)))'...
                    [squeeze(nanmean(p_pred_w_cvar(:,cond,Nind,end:-1:1),1))]'+[squeeze(nanstd(p_pred_w_cvar(:,cond,Nind,end:-1:1),1))]'./sqrt(squeeze(sum(~isnan(p_pred_w_cvar(:,cond,Nind,:)),1)))'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
            end
            if model_pred
                plot(bincenterz_cvar, squeeze(nanmean(p_data_w_cvar(:,cond,Nind,:),1)), 'o','Color', colors_sz(Nind,:), 'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
                leg(Nind) = errorbar(bincenterz_cvar,squeeze(nanmean(p_data_w_cvar(:,cond,Nind,:),1)), squeeze(nanstd(p_data_w_cvar(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none', 'Linewidth',1.3);
            else
                plot(bincenterz_cvar, squeeze(nanmean(p_data_w_cvar(:,cond,Nind,:),1)), 'o','Color', colors_sz(Nind,:), 'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
                leg(Nind) = errorbar(bincenterz_cvar,squeeze(nanmean(p_data_w_cvar(:,cond,Nind,:),1)), squeeze(nanstd(p_data_w_cvar(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linewidth',1.3);
            end
            plot(0:0.10:0.90, ones(1,10)*1/Nvec(Nind),'--', 'Color',colors_sz(Nind,:), 'Linewidth', 1.3); hold on;
            box off
            ylim([0 1])
            xlim([0 1])
            set(gca, 'xtick',[0:0.2:1])
            set(gca, 'ytick', [0:0.1:1])
            set(gca, 'yticklabels', {})
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', ticklengthsz)
            
            
            xlabel('Distractor variance','FontName', 'Helvetica', 'FontSize', fontsz)
            
            if cond == 2
                set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
                ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            
        end
        if cond == 4
            llpa = legend([legg(1) legg(2) legg(3) legg(4) chance],'N = 2', 'N = 3','N = 4', 'N = 6','chance','FontName','Helvetica','FontSize', fontsz)
            legend boxoff
            set(llpa, 'FontName','Helvetica','FontSize', fontsz*0.9, 'Position', [0.9 0.84 0.04 0.04])
        end
        
        
        
        %e) prop responses
        guttera3=[0.04 0.09]; %horiz vert
        marginsa3=[0.08 0.08 0.06 0.09]; %left right bottom top
        for Nind = 2:4
            
            tight_subplot(5,6, 5, 3*(cond/2-1) + Nind-1, guttera3, marginsa3)
            
            if model_predR
                h_p=fill([1:1:Nvec(Nind)-1 (Nvec(Nind)-1):-1:1], [[[squeeze(mean(p_resp_ori_diff_rank_simim(:,cond,Nind,1:Nvec(Nind)-1),1))]-[squeeze(std(p_resp_ori_diff_rank_simim(:,cond,Nind,1:Nvec(Nind)-1),1))]/sqrt(Nsubj)]'...
                    [[squeeze(mean(p_resp_ori_diff_rank_simim(:,cond,Nind,(Nvec(Nind)-1):-1:1),1))]+[squeeze(std(p_resp_ori_diff_rank_simim(:,cond,Nind,Nvec(Nind)-1:-1:1),1))]/sqrt(Nsubj)]'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
                %bar(1:1:(max(loc_vec)+0), squeeze(mean(p_resp_ori_diff_rank_simim(:,cond,Nind,1:(max(loc_vec)+0)),1)), 'EdgeColor', colors_sz(Nind,:), 'FaceColor', colors_sz_shade(Nind,:), 'BarWidth',0.8); hold on;
            end
            
            
            bar(1:1:Nvec(Nind)-1, squeeze(mean(p_resp_ori_diff_rank(:,cond,Nind,1:Nvec(Nind)-1),1)), 'EdgeColor', colors_sz(Nind,:),'LineWidth', 0.6,  'FaceColor', 'None', 'BarWidth',0.46); hold on;
            errorbar(1:1:Nvec(Nind)-1,squeeze(mean(p_resp_ori_diff_rank(:,cond,Nind,1:Nvec(Nind)-1),1)), squeeze(std(p_resp_ori_diff_rank(:,cond,Nind,1:Nvec(Nind)-1),1))/sqrt(Nsubj),'Color', colors_sz(Nind,:), 'LineStyle', 'none', 'Linewidth',1.3);
            if Nind>1
                plot(1:0.2:Nvec(Nind)-1, (1-mean(p_c(:,cond, Nind),1))*1/(Nvec(Nind)-1), '.','Color',  'k','MarkerSize',4,  'Linewidth', 1.3); hold on;
            end
            
            
            box off
            ylim([0 p_resp_max])
            xlim([0.4 max(Nvec)-1+0.6])
            set(gca, 'xtick', (0:1:max(Nvec)-1))
            set(gca, 'ytick', [0:0.1:0.4])
            if Nind>2
                set(gca, 'yticklabels', []);
            end
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', ticklengthsz)
            
            
            if cond ==2 & Nind == 4
                %xlabel('Rank of T - R orientation distance','FontName', 'Helvetica', 'FontSize', fontsz)
                xlabel('Rank of similarity of response to the target','FontName', 'Helvetica', 'FontSize', fontsz)
            elseif Nind == 2
                %ylabel('Proportion response','FontName', 'Helvetica', 'FontSize', fontsz)
                ylabel('Frequency of occurence','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            
        end
        
        
    else
        
        %b) T- MSD
        for Nind = 1:4
            bincentersl=[squeeze(bincenterz(cond, Nind,1,:))]';
            
            
            tight_subplot(1,2,1,cond/2, guttera2, marginsa2)
            
            if model_pred
                h_p=fill([bincentersl bincentersl(end:-1:1)], [[squeeze(mean(p_pred_w(:,cond,Nind,:),1))]'-[squeeze(std(p_pred_w(:,cond,Nind,:),1))]'./sqrt(squeeze(sum(~isnan(p_pred_w(:,cond,Nind,:)),1)))'...
                    [squeeze(mean(p_pred_w(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(p_pred_w(:,cond,Nind,end:-1:1),1))]'./sqrt(squeeze(sum(~isnan(p_pred_w(:,cond,Nind,:)),1)))'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
                
            end
            if model_pred
                plot(bincentersl, squeeze(mean(p_data_w(:,cond,Nind,:),1)), 'o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot,'Linestyle', 'none'); hold on;
                leg(Nind) = errorbar(bincentersl,squeeze(mean(p_data_w(:,cond,Nind,:),1)), squeeze(std(p_data_w(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none', 'Linewidth',1.3);
            else
                plot(bincentersl, squeeze(mean(p_data_w(:,cond,Nind,:),1)), 'o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot); hold on;
                leg(Nind) = errorbar(bincentersl,squeeze(mean(p_data_w(:,cond,Nind,:),1)), squeeze(std(p_data_w(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linewidth',1.3);
            end
            chance =  plot(0:10:90, ones(1,10)*1/Nvec(Nind),'--', 'Color',colors_sz(Nind,:), 'Linewidth', 1.3); hold on;
            box off
            ylim([0 1])
            xlim([0 90])
            set(gca, 'xtick', [0:15:90])
            set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
            set(gca, 'ytick', [0:0.1:1])
            set(gca, 'yticklabels', {})
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', ticklengthsz)
            if cond == 2
                set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
                ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            xlabel('min T-D difference (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
            
        end
        if Nind == 1
            text(30, 1.16, text_labels{(cond)/2},  'FontName', 'Helvetica', 'FontSize', fontsz)
        end
        if cond == 2
            text(-40, 1.23, 'Localization',  'FontName', 'Helvetica', 'FontSize', fontsz*1.2)
        end
        if cond == 4
            llpa = legend([leg(1) leg(2) leg(3) leg(4) chance],'N = 2', 'N = 3','N = 4', 'N = 6','chance','FontName','Helvetica','FontSize', fontsz)
            legend boxoff
            set(llpa, 'FontName','Helvetica','FontSize', fontsz*0.9, 'Position', [0.9 0.84 0.04 0.04])
        end
        
    end
    
    
end

psname = ['LOC_full_figure_',num2str(full_figure_flag),'_model_',num2str(model_pred),'_model_index_',num2str(mi),'_type_',num2str(type) ,'_exp_',num2str(exp_i), '.pdf']
if save_flag
    print_pdf(psname)
end


end