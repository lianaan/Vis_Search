clear all; close all;
Nvec = [2 3 4 6 10 14];
Ntrials = 100000;
% in loc, T is always in the display
% in det, not always
kappa_prior = 0;  % kappa of the distribution from which to draw stimulus orientations
% Josh: conc parameter 1.5
%kappa_prior = 1.5;
%loc
figure
set(gcf, 'Position', [100 100 200 480])
set_plotting_params();
colors_szz = [colors_sz; [128 128 128]/256; [28 28 28]/256];
fontsz = 12;
guttera = [0.04 0.09];
marginsa = [0.19 0.1 0.14 0.1]; %left right bottom top


for Nind = 2: length(Nvec)
    stimvec = nan(Ntrials, (Nvec(Nind)));
    for Nti = 1:Ntrials
        stimvec(Nti,1: Nvec(Nind)) = circ_vmrnd(zeros(1,1),kappa_prior,Nvec(Nind)); %* 2;%/pi*90 ;
        targetidx(Nti) = randi(Nvec(Nind));
        targetval(Nti) = stimvec(Nti,targetidx(Nti));
        
        
        vall = [abs(bsxfun(@minus,stimvec(Nti,:)',targetval(Nti)))]'; %btwn (0,2*pi)
        vall(vall>pi) = 2*pi-vall(vall>pi);
        min_ori_diff(Nind,Nti) = min(vall(vall>0)); % btwn 0 and pi
        
        
        vall_mean = circ_mean(setdiff(stimvec(Nti,find(~isnan(stimvec(Nti,:)))), targetval(Nti))');
        vall_mean = abs(bsxfun(@minus,vall_mean, targetval(Nti)));
        vall_mean(vall_mean>pi) = 2*pi- vall_mean(vall_mean>pi);
        mean_ori_diff(Nind,Nti) = vall_mean;% in the interval (0,pi)
        
        vall_cvar_vec = (setdiff(stimvec(Nti,find(~isnan(stimvec(Nti,:)))), targetval(Nti)));
        vall_cvar_vec(vall_cvar_vec>pi) = 2*pi-vall_cvar_vec(vall_cvar_vec>pi);
        vall_cvar = circ_var(vall_cvar_vec,[],[],[]);
        cvar_ori_diff(Nind,Nti) = vall_cvar; % in the interval (0,1)
    end
    
    
    edges = linspace(0, pi, 100);
    %edges = linspace(0, pi, Ntrials);
    bin_size = edges(2)-edges(1);
    
    tight_subplot(3,1,1, 1,guttera, marginsa)
    [N1,BIN1] = histc(min_ori_diff(Nind,:), edges); hold on;
    %plot(180/pi* 1/2 * edges, N1/Ntrials/bin_size, 'Color', colors_szz(Nind,:), 'Linewidth',1.3); hold on;
    plot(180/pi* 1/2 * edges, N1/Ntrials, 'Color', colors_szz(Nind,:), 'Linewidth',1.3); hold on;
    xlim([0 90])
    %xlim([0 pi])
    %ylim([0 1.02*max(N1/Ntrials)])
    %ylim([0 15000])
    ylim([0 0.15])
    set(gca, 'tickdir', 'out')
    box off
    %title(['N=', num2str(Nvec(Nind))])
    set(gca, 'xtick', [0:15:90])
    set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
    if Nind == 2
        xlabel('min T-D difference (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
        %ylabel('Counts','FontName', 'Helvetica', 'FontSize', fontsz)
        %ylabel('Proportion','FontName', 'Helvetica', 'FontSize', fontsz)
        ylabel('Probability density','FontName', 'Helvetica', 'FontSize', fontsz)
    end
    
    tight_subplot(3,1,2, 1,guttera, marginsa)
    [N2,BIN2] = histc(mean_ori_diff(Nind,:), edges); hold on;
    plot(180/pi* 1/2 * edges, N2/Ntrials, 'Color', colors_szz(Nind,:), 'Linewidth',1.3); hold on;
    xlim([0 90])
    %xlim([0 pi])
    %ylim([0 1.15*max(N2/Ntrials)])
    %ylim([0 15000])
    ylim([0 0.15])
    set(gca, 'tickdir', 'out')
    box off
    set(gca, 'xtick', [0:15:90])
    set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
    if Nind == 2
        xlabel(' T-D mean (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
        %ylabel('counts','FontName', 'Helvetica', 'FontSize', fontsz)
        %ylabel('Proportion','FontName', 'Helvetica', 'FontSize', fontsz)
        ylabel('Probability density','FontName', 'Helvetica', 'FontSize', fontsz)
    end
    
    tight_subplot(3,1,3, 1,guttera, marginsa)
    [N3,BIN3] = histc(cvar_ori_diff(Nind,:), edges); hold on;
    plot(edges, N3/Ntrials, 'Color', colors_szz(Nind,:), 'Linewidth',1.3); hold on;
    xlim([0 1])
    %ylim([0 1.15*max(N3/Ntrials)])
    %ylim([0 15000])
    ylim([0 0.15])
    set(gca, 'tickdir', 'out')
    box off
    set(gca, 'xtick',[0:0.2:1])
    if Nind == 2
        xlabel('Distractor variance','FontName', 'Helvetica', 'FontSize', fontsz)
        %ylabel('counts','FontName', 'Helvetica', 'FontSize', fontsz)
        %ylabel('Proportion','FontName', 'Helvetica', 'FontSize', fontsz)
        ylabel('Probability density','FontName', 'Helvetica', 'FontSize', fontsz)
    end
    
    [corr12(Nind) pval12(Nind)] = circ_corrcc(min_ori_diff(Nind,:)', mean_ori_diff(Nind,:)');%, 'type', 'Spearman');
    
    [corr13(Nind) pval13(Nind)] = circ_corrcc(min_ori_diff(Nind,:)', cvar_ori_diff(Nind,:)');%, 'type', 'Spearman');
    
    [corr23(Nind) pval23(Nind)] = circ_corrcc(mean_ori_diff(Nind,:)', cvar_ori_diff(Nind,:)');%, 'type', 'Spearman');
    
    %if Nind == 2
    text_loc_ref = max(N3);
    %end
    %{
    text(0,-0.8*text_loc_ref, ['\rho_{12} = ', num2str(corr12(Nind), '%.2f'),', p = ', num2str(pval12(Nind), '%.2f')], 'FontSize', 9)
    text(0,-1*text_loc_ref, ['\rho_{13} = ', num2str(corr13(Nind), '%.2f'),', p = ', num2str(pval13(Nind), '%.2f')], 'FontSize', 9)
    text(0,-1.2*text_loc_ref, ['\rho_{23} = ', num2str(corr23(Nind), '%.2f'),', p = ', num2str(pval23(Nind), '%.2f')], 'FontSize', 9)
    %}
end
%%
psname = 'sum_stats_distr_loc_U_100000E4.pdf'
%print_pdf(psname)

