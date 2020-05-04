clear all; close all;
exp_i = 1;
load(['mdlALL_exp', num2str(exp_i),'.mat'])

if exp_i == 1
    Nsubj = 11;
elseif exp_i == 2
    Nsubj = 6;
end
Ncond = 4;

indi_rm = {[1],[1 2], [1 3], [1 4], [1 2 3 4], [1 3 4], [1 2 4], [1 2 3]};
coeff_len = [2 4 4 4 16 8 8 8] - 1; %  substract the intercept
var_indi_rm = {'N','min','mean', 'var','acc'};

regr_models_beta = nan(Nsubj, Ncond,length(indi_rm),max(coeff_len));
regr_models_beta_se = nan(Nsubj, Ncond,length(indi_rm),max(coeff_len));
regr_models_tStat = nan(Nsubj, Ncond,length(indi_rm),max(coeff_len));
regr_models_pValue = nan(Nsubj, Ncond,length(indi_rm), max(coeff_len));

regr_models_LL = nan(Nsubj, Ncond,length(indi_rm));
regr_models_AIC = nan(Nsubj, Ncond,length(indi_rm));
regr_models_BIC = nan(Nsubj, Ncond,length(indi_rm));
regr_models_Rsquared_adjusted = nan(Nsubj, Ncond,length(indi_rm));


for sbjid = 1: size(mdlALL,1)
    for cond = 1: size(mdlALL,2)
        for rmi = 1: size(mdlALL,3)
            regr_models_beta(sbjid, cond,rmi, 1:coeff_len(rmi)) = mdlALL{sbjid,cond,rmi}.Coefficients.Estimate(2:end); % do
            %not include initial slope
            regr_models_beta_se(sbjid, cond,rmi, 1:coeff_len(rmi)) = mdlALL{sbjid,cond,rmi}.Coefficients.SE(2:end);
            regr_models_beta_tStat(sbjid, cond,rmi, 1:coeff_len(rmi)) = mdlALL{sbjid,cond,rmi}.Coefficients.tStat(2:end);
            regr_models_beta_pValue(sbjid, cond,rmi, 1:coeff_len(rmi)) = mdlALL{sbjid,cond,rmi}.Coefficients.pValue(2:end);
            
            regr_models_LL(sbjid,cond,rmi) = mdlALL{sbjid,cond,rmi}.LogLikelihood;
            regr_models_AIC(sbjid,cond,rmi) = mdlALL{sbjid,cond,rmi}.ModelCriterion.AIC;
            regr_models_BIC(sbjid,cond,rmi) = mdlALL{sbjid,cond,rmi}.ModelCriterion.BIC;
            regr_models_Rsquared_adjusted(sbjid, cond, rmi) = mdlALL{sbjid,cond,rmi}.Rsquared.Adjusted;
        end
    end
end
%%
close all;
figure
eb_w = 0.12;%errorbar width
eb_t = 0.24; %errorbar line thickness
fontsz = 13;
sz_dot = 1.2;

cond_order = [2 4 1 3];
guttera = [0.06 0.09];
marginsa = [0.08 0.06 0.12  0.11]; %left right bottom top
log_regr_col = [255 255 255; 77 0 75; 140 107 177; 158 188 218]/255;

set(gcf, 'Position', [100 100 550 340])
diff_regr_models_AIC = NaN(size(regr_models_AIC));
diff_regr_models_AIC(:,:,1:4) = bsxfun(@minus,regr_models_AIC(:,:,1), regr_models_AIC(:,:,1:4) );
diff_regr_models_AIC(:,:,5:8) = bsxfun(@minus, regr_models_AIC(:,:,5:8), regr_models_AIC(:,:,5));


diff_regr_models_BIC = NaN(size(regr_models_BIC));
diff_regr_models_BIC(:,:,1:4) = bsxfun(@minus,regr_models_BIC(:,:,1), regr_models_BIC(:,:,1:4) );
diff_regr_models_BIC(:,:,5:8) = bsxfun(@minus, regr_models_BIC(:,:,5:8), regr_models_BIC(:,:,5));


for cond = 1 : Ncond
    
    tight_subplot(4,4,1,cond, guttera, marginsa)
    for rmi = 1 :4% length(indi_rm)
        bar(rmi, mean(diff_regr_models_AIC(:,cond_order(cond),rmi),1), std(diff_regr_models_AIC(:,cond_order(cond),rmi),1), 'FaceColor', log_regr_col(rmi,:),'EdgeColor',log_regr_col(rmi,:), 'BarWidth', 0.5); hold on;
        ei(rmi) = errorbar(rmi,  mean(diff_regr_models_AIC(:,cond_order(cond),rmi),1),std(diff_regr_models_AIC(:,cond_order(cond),rmi),1)/sqrt(Nsubj),...
            'o','MarkerSize', sz_dot, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color', 'k', 'LineWidth',1.3); hold on;
        errorbarT(ei(rmi),eb_w,eb_t)
        box off
        
    end
    xlim([0.5 4.5])
    ylim([-50 70])
    set(gca, 'tickdir', 'out')
    set(gca,'ticklength', [0.025 0.025] )
    set(gca, 'xtick', [1:1:4])
    set(gca, 'xticklabels', {'1','2','3','4','5','6','7','8'},'FontName','Helvetica', 'FontSize', 0.7*fontsz);
    set(gca, 'ytick', [-50:30:70])
    
    
    if cond == 1
        %ylabel( '\Delta AIC', 'FontName','Helvetica', 'FontSize', 0.9*fontsz)
        ylabel( 'KID AIC', 'FontName','Helvetica', 'FontSize', 0.9*fontsz)
    end
    
    tight_subplot(4,4,2,cond, guttera, marginsa)
    for rmi = 1 :4% length(indi_rm)
        bar(rmi, mean(diff_regr_models_BIC(:,cond_order(cond),rmi),1), std(diff_regr_models_BIC(:,cond_order(cond),rmi),1), 'FaceColor', log_regr_col(rmi,:),'EdgeColor',log_regr_col(rmi,:), 'BarWidth', 0.5); hold on;
        ei(rmi) = errorbar(rmi,  mean(diff_regr_models_BIC(:,cond_order(cond),rmi),1),std(diff_regr_models_BIC(:,cond_order(cond),rmi),1)/sqrt(Nsubj),...
            'o','MarkerSize', sz_dot, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color', 'k', 'LineWidth',1.3); hold on;
        
        errorbarT(ei(rmi),eb_w,eb_t)
        box off
        
    end
    xlim([0.5 4.5])
    ylim([-50 70])
    set(gca, 'tickdir', 'out')
    set(gca,'ticklength', [0.025 0.025] )
    set(gca, 'xtick', [1:1:4])
    set(gca, 'xticklabels', {'1','2','3','4','5','6','7','8'},'FontName','Helvetica', 'FontSize', 0.7*fontsz);
    set(gca, 'ytick', [-50:30:70])
    if cond == 1
        ylabel( 'KID BIC', 'FontName','Helvetica', 'FontSize', 0.9*fontsz)
    end
    
    
    tight_subplot(4,4,3,cond, guttera, marginsa)
    for rmi = 5 : length(indi_rm)
        bar(rmi, mean(diff_regr_models_AIC(:,cond_order(cond),rmi),1), std(diff_regr_models_AIC(:,cond_order(cond),rmi),1), 'FaceColor', log_regr_col(rmi-4,:),'EdgeColor',log_regr_col(rmi-4,:), 'BarWidth', 0.5); hold on;
        ei(rmi) = errorbar(rmi,  mean(diff_regr_models_AIC(:,cond_order(cond),rmi),1),std(diff_regr_models_AIC(:,cond_order(cond),rmi),1)/sqrt(Nsubj),...
            'o','MarkerSize', sz_dot, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color', 'k', 'LineWidth',1.3); hold on;
        errorbarT(ei(rmi),eb_w,eb_t)
        box off
        
        
    end
    xlim([4.5 8.5])
    ylim([-50 70])
    set(gca, 'tickdir', 'out')
    set(gca,'ticklength', [0.025 0.025] )
    set(gca, 'xtick', [5:1:8])
    set(gca, 'xticklabels', {'1','2','3','4'},'FontName','Helvetica', 'FontSize', 0.7*fontsz);
    set(gca, 'ytick', [-50:30:70])
    if cond == 1
        ylabel( 'KOD AIC', 'FontName','Helvetica', 'FontSize', 0.9*fontsz)
    end
    
    tight_subplot(4,4,4,cond, guttera, marginsa)
    for rmi = 5 : length(indi_rm)
        bar(rmi, mean(diff_regr_models_BIC(:,cond_order(cond),rmi),1), std(diff_regr_models_BIC(:,cond_order(cond),rmi),1), 'FaceColor', log_regr_col(rmi-4,:),'EdgeColor',log_regr_col(rmi-4,:), 'BarWidth', 0.5); hold on;
        ei(rmi) = errorbar(rmi,  mean(diff_regr_models_BIC(:,cond_order(cond),rmi),1),std(diff_regr_models_BIC(:,cond_order(cond),rmi),1)/sqrt(Nsubj),...
            'o','MarkerSize', sz_dot, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color', 'k', 'LineWidth',1.3); hold on;
        
        errorbarT(ei(rmi),eb_w,eb_t)
        box off
        
    end
    xlim([4.5 8.5])
    ylim([-50 70])
    set(gca, 'tickdir', 'out')
    set(gca,'ticklength', [0.025 0.025] )
    set(gca, 'xtick', [5:1:8])
    set(gca, 'xticklabels', {'1','2','3','4'},'FontName','Helvetica', 'FontSize', 0.7*fontsz);
    set(gca, 'ytick', [-50:30:70])
    
    xlabel('Regressor', 'FontName','Helvetica', 'FontSize', 0.9*fontsz)
    
    if cond == 1
        ylabel( 'KOD BIC', 'FontName','Helvetica', 'FontSize', 0.9*fontsz)
    end
    
end
%%
psname =  ['Regr_models_KID_KOD_comparison_Exp',num2str(exp_i),'.pdf']
%print_pdf(psname)


%%


sample=[];
sample2=[];
nboot = 100000;

for cond = 1: 4
    for rmi = [2:4 6:8]
        
        d_AIC_sum(cond_order(cond),rmi) = sum(squeeze(diff_regr_models_AIC(:,cond_order(cond),rmi)));
        d_BIC_sum(cond_order(cond),rmi) = sum(squeeze(diff_regr_models_BIC(:,cond_order(cond),rmi)));
        for kk=1:nboot
            sample=randsample(diff_regr_models_AIC(:,cond_order(cond),rmi),Nsubj,1);
            d_AIC_sums(kk,cond_order(cond), rmi) = sum(sample);
            
            sample2=randsample(diff_regr_models_BIC(:,cond_order(cond),rmi),Nsubj,1);
            d_BIC_sums(kk,cond_order(cond), rmi) = sum(sample2);
        end
    end
end
%%
ci_bnd_low = 0.025;
ci_bnd_high = 0.975;
for cond = 1:4
    for rmi = [2:4 6:8]
        bci_aic(cond,rmi,1:2) = [quantile(squeeze(d_AIC_sums(:,cond,rmi)),ci_bnd_low); quantile(squeeze(d_AIC_sums(:,cond,rmi)),ci_bnd_high)];
        bci_bic(cond,rmi,1:2) = [quantile(squeeze(d_BIC_sums(:,cond,rmi)),ci_bnd_low); quantile(squeeze(d_BIC_sums(:,cond,rmi)),ci_bnd_high)];
    end
end

%%
figure
set(gcf, 'Position', [100 100 450 490])
if exp_i == 1
    ylim_min1 = -200; ylim_max1 = 800;
    ylim_min2 = -500; ylim_max2 = 190;
    y_step = 300;
elseif exp_i == 2
    ylim_min1 = -100; ylim_max1 = 400;
    ylim_min2 = -350; ylim_max2 = 50;
    y_step = 200;
end

guttera = [0.06 0.05];%[0.08 0.07];
marginsa = [0.12 0.06 0.12  0.11];
for cond = 1 : Ncond
    
    tight_subplot(4,4,1,cond, guttera, marginsa)
    for rmi = 1 :4% length(indi_rm)
        hb(rmi)=bar(rmi, d_AIC_sum(cond_order(cond),rmi), 'FaceColor', log_regr_col(rmi,:),'EdgeColor',log_regr_col(rmi,:), 'BarWidth', 0.5); hold on;
        ei(rmi) = errorbar(rmi, d_AIC_sum(cond_order(cond),rmi), d_AIC_sum(cond_order(cond),rmi)-bci_aic(cond_order(cond),rmi,1),bci_aic(cond_order(cond),rmi,2)-d_AIC_sum(cond_order(cond),rmi),...
            'o','MarkerSize', sz_dot, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color', 'k', 'LineWidth',1.3); hold on;
        
        errorbarT(ei(rmi),eb_w,eb_t)
        box off
        
    end
    xlim([0.5 4.5])
    ylim([ylim_min1 ylim_max1])
    set(gca, 'tickdir', 'out')
    set(gca,'ticklength', [0.025 0.025] )
    set(gca, 'xtick', [1:1:4])
    set(gca, 'xticklabels', [])
    set(gca, 'ytick', [ylim_min1:y_step:ylim_max1])
   
    if cond == 1
        ylabel( 'KID AIC', 'FontName','Helvetica', 'FontSize', 0.9*fontsz)
    end
    if cond == 4
        ll=legend([hb(2) hb(3) hb(4)],'min T-D difference', 'T-D mean', 'distractor variance')
        set(ll, 'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
        set(ll, 'Position',[0.66 0.92 0.24 0.06])
    end
    
    tight_subplot(4,4,2,cond, guttera, marginsa)
    for rmi = 1 :4% length(indi_rm)
        bar(rmi, d_BIC_sum(cond_order(cond),rmi), 'FaceColor', log_regr_col(rmi,:),'EdgeColor',log_regr_col(rmi,:), 'BarWidth', 0.5); hold on;
        ei(rmi) = errorbar(rmi, d_BIC_sum(cond_order(cond),rmi) ,d_BIC_sum(cond_order(cond),rmi)-bci_bic(cond_order(cond),rmi,1),bci_bic(cond_order(cond),rmi,2)-d_BIC_sum(cond_order(cond),rmi),...
            'o','MarkerSize', sz_dot, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color', 'k', 'LineWidth',1.3); hold on;
        
        errorbarT(ei(rmi),eb_w,eb_t)
        box off
        
    end
    xlim([0.5 4.5])
    ylim([ylim_min1 ylim_max1])
    set(gca, 'tickdir', 'out')
    set(gca,'ticklength', [0.025 0.025] )
    set(gca, 'xtick', [1:1:4])
    set(gca, 'xticklabels', [])
    set(gca, 'ytick', [ylim_min1:y_step:ylim_max1])
   
    if cond == 1
        ylabel( 'KID BIC', 'FontName','Helvetica', 'FontSize', 0.9*fontsz)
    end
    
    
    tight_subplot(4,4,3,cond, guttera, marginsa)
    for rmi = 5 : length(indi_rm)
        bar(rmi, d_AIC_sum(cond_order(cond),rmi)  , 'FaceColor', log_regr_col(rmi-4,:),'EdgeColor',log_regr_col(rmi-4,:), 'BarWidth', 0.5); hold on;
        ei(rmi) = errorbar(rmi, d_AIC_sum(cond_order(cond),rmi), d_AIC_sum(cond_order(cond),rmi)-bci_aic(cond_order(cond),rmi,1),bci_aic(cond_order(cond),rmi,2)-d_AIC_sum(cond_order(cond),rmi),...
            'o','MarkerSize', sz_dot, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color', 'k', 'LineWidth',1.3); hold on;
        
        errorbarT(ei(rmi),eb_w,eb_t)
        box off
        
        
    end
    xlim([4.5 8.5])
    ylim([ylim_min2 ylim_max2])
    set(gca, 'tickdir', 'out')
    set(gca,'ticklength', [0.025 0.025] )
    set(gca, 'xtick', [5:1:8])
    set(gca, 'xticklabels', [])
    set(gca, 'ytick', [ylim_min2:y_step:ylim_max2])

    if cond == 1
        ylabel( 'KOD AIC', 'FontName','Helvetica', 'FontSize', 0.9*fontsz)
    end
    
    tight_subplot(4,4,4,cond, guttera, marginsa)
    for rmi = 5 : length(indi_rm)
        bar(rmi, d_BIC_sum(cond_order(cond),rmi) , 'FaceColor', log_regr_col(rmi-4,:),'EdgeColor',log_regr_col(rmi-4,:), 'BarWidth', 0.5); hold on;
        ei(rmi) = errorbar(rmi, d_BIC_sum(cond_order(cond),rmi), d_BIC_sum(cond_order(cond),rmi)-bci_bic(cond_order(cond),rmi,1),bci_bic(cond_order(cond),rmi,2)- d_BIC_sum(cond_order(cond),rmi),...
            'o','MarkerSize', sz_dot, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color', 'k', 'LineWidth',1.3); hold on;
        
        errorbarT(ei(rmi),eb_w,eb_t)
        box off
        
    end
    xlim([4.5 8.5])
    ylim([ylim_min2 ylim_max2])
    set(gca, 'tickdir', 'out')
    set(gca,'ticklength', [0.025 0.025] )
    set(gca, 'xtick', [5:1:8])
    set(gca, 'xticklabels', [])
    set(gca, 'xticklabels', [])
    set(gca, 'ytick', [ylim_min2:y_step:ylim_max2])

    if cond == 1
        ylabel( 'KOD BIC', 'FontName','Helvetica', 'FontSize', 0.9*fontsz)
    end
    
end
%%
psname =  ['Regr_models_KID_KOD_E',num2str(exp_i),'.pdf']
%print_pdf(psname)