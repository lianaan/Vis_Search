clear all; close all;

Nvec = [2 3 4 6];
colors_sz = [217 30 0 ;  51 164 0;  242 135 65; 130 44 169]./255; % red, green, orange, purple
fontsz = 13;

fcn_quantile = @(theta,N) 1 - (1-theta/pi).^N;
fcn_bin_placement = @(quantile_val,N) pi *(1-(1-quantile_val).^(1/N));

nbinz = 6;
theta_range = linspace(0, pi,2*nbinz);
figure
set(gcf, 'Position', [100 100 400 200])
for Nind = 1: length(Nvec)
    N = Nvec(Nind);
    plot(theta_range * 180/pi* 1/2, fcn_quantile(theta_range,N), '-',  'Color',colors_sz(Nind,:),'Linewidth', 1.3); hold on;
    plot(fcn_bin_placement(linspace(0, 1,nbinz),N)* 180/pi * 1/2,  linspace(0, 1,nbinz) ,'o',  'MarkerFaceColor',colors_sz(Nind,:), 'MarkerEdgeColor',colors_sz(Nind,:)); hold on;
    set(gca, 'xtick', [0:15:90])
    set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', 0.9*fontsz)
    
    for ij = 4%1:nbinz
         plot(fcn_bin_placement(1/(nbinz-1)*ij,N)* 180/pi* 1/2 * ones(1,10), linspace(0,1/(nbinz-1)*ij,10), '--','Color',colors_sz(Nind,:)); hold on;
    end
    
end
for ij = 1: nbinz
    plot( [0:15:90],1/(nbinz-1)*ij* ones(1,length([0:15:90])),'--k'); hold on;
end
xlim([0 90])
ylim([0 1.01])
set(gca, 'tickdir', 'out')
box off
xlabel('min T-D difference (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
ylabel('Cumulative probability', 'FontName', 'Helvetica', 'FontSize', fontsz)

psname = 'ideal_binning.pdf';
%print_pdf(psname)
            