clear all; close all;
addpath('CircStat2012a')

Nvec = 2; 
Npts = 1200;
kappa_prior = 0;  % kappa of the distribution from which to draw stimulus orientations

Jset = exp([ 3 3; 3 1; 1 1]); 

figure
set(gcf, 'Position', [100 100 520 340])

fontsz = 11;
guttera = [0.08 0.13];
marginsa = [0.1 0.1 0.1 0.14]; %left right bottom top

for cond = 1:2
    for indir = 1:3
        clear dec_var1; clear dec_var2; clear dec_var_glob; clear dec_var_loc1;
        
        J         = repmat(Jset(indir,:), Npts,1);%gamrnd(J1bar(Nind)/tau, tau, length(idx), N, Nsamp); % Precisions: Nidx x N x Nsamp
        kappa     = fisher2kappa(J); % Concentration parameter of von Mises
        
        
        if  cond == 1
            
            tight_subplot(2,3,3-cond,indir, guttera, marginsa)
            
            xt1 = linspace(-pi, pi, Npts);
            xt2 = linspace(-pi, pi, Npts);
            for k1 = 1 : Npts
                dec_var1(k1) =  -log(besseli(0,kappa(k1,1),1)) - kappa(k1,1) + kappa(k1,1).*cos(xt1(k1));
                for k2 = 1 : Npts
                    dec_var2(k2) =  -log(besseli(0,kappa(k2,2),1)) - kappa(k2,2) + kappa(k2,2).*cos(xt2(k2));
                    dec_var_glob(k1,k2) = squeeze(log(mean(exp([dec_var1(k1) dec_var2(k2)]),2))); %+ log(pp/(1-pp)));
                end
            end
            imagesc(xt1/2,xt2/2, (dec_var_glob > 0)'); hold on;
            colormap(gray)
            set(gca, 'yDir', 'Normal')
            box off
            set(gca, 'tickdir', 'out')
            xlim([-pi/2 pi/2])
            ylim([-pi/2 pi/2])
            set(gca, 'xtick',[-pi/2: pi/6:pi/2])
            set(gca, 'xticklabels',{'-90', '-60', '-30', '0', '30', '60', '90'},  'FontName', 'Helvetica', 'FontSize', 0.8*fontsz)
            set(gca, 'ytick',[-pi/2: pi/6:pi/2])
            set(gca, 'yticklabels',{'-90', '-60', '-30', '0', '30', '60', '90'},  'FontName', 'Helvetica', 'FontSize', 0.8*fontsz)
            xlabel('measurement 1 - target ', 'FontName', 'Helvetica', 'FontSize', fontsz)
            if indir == 1
                ylabel('measurement 2 - target ', 'FontName', 'Helvetica', 'FontSize', fontsz)
                title('Decision boundary for target present vs absent','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            %}
        elseif cond == 2
            tight_subplot(2,3,3-cond,indir, guttera, marginsa)
            xt1 = linspace(-pi, pi, Npts);
            xt2 = linspace(-pi, pi, Npts);
            for k1 = 1 : Npts
                dec_var1(k1) = -log(besseli(0,kappa(k1,1),1)) - kappa(k1,1) + kappa(k1,1).*cos(xt1(k1));
                for k2 = 1 : Npts
                    dec_var2(k2) = -log(besseli(0,kappa(k2,2),1)) - kappa(k2,2) + kappa(k2,2).*cos(xt2(k2));
                    dec_var_loc1(k1,k2) = double(dec_var1(k1) > dec_var2(k2)); %+ log(pp/(1-pp)));
                end
            end
            
            imagesc(xt1/2,xt2/2,  (dec_var_loc1)') ; hold on;
            colormap(gray)
            set(gca, 'yDir', 'Normal')
            box off
            set(gca, 'tickdir', 'out')
            xlim([-pi/2 pi/2])
            ylim([-pi/2 pi/2])
            set(gca, 'xtick',[-pi/2: pi/6:pi/2])
            set(gca, 'xticklabels',{'-90', '-60', '-30', '0', '30', '60', '90'},  'FontName', 'Helvetica', 'FontSize', 0.8*fontsz)
            set(gca, 'ytick',[-pi/2: pi/6:pi/2])
            set(gca, 'yticklabels',{'-90', '-60', '-30', '0', '30', '60', '90'},  'FontName', 'Helvetica', 'FontSize', 0.8*fontsz)
            text( -1.7, 2.2, ['log J1 = ', num2str(log(Jset(indir,1)), '%.1f'), ', log J2 = ', num2str(log(Jset(indir,2)), '%.1f')])
            
            if indir == 1
                ylabel('measurement 2 - target ', 'FontName', 'Helvetica', 'FontSize', fontsz)
                title('Decision boundary for location 1 vs 2','FontName', 'Helvetica', 'FontSize', fontsz)
            end
        end
    end
end
%%
psname = 'Decision_boundaries.pdf'
%print_pdf(psname)

