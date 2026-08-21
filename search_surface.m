%clear all;
close all;


T_NT_sim = 90/pi*abs(squeeze(mean_ori_diff(Nind, :)));   %since we doubled all our stimuli for modelling on the (-pi,pi) full space
NT_NT_sim = squeeze(cvar_ori_diff(Nind,:));

nbinz = 7;
binz_1 = [];
binz_2 = [];


binz_1 = [linspace(min(T_NT_sim), max(T_NT_sim), nbinz+1) ];
binz_posE_1= (binz_1(2:end)+binz_1(1:end-1))/2;
binz_2 = [linspace(min(NT_NT_sim), max(NT_NT_sim), nbinz+1)];
binz_posE_2= (binz_2(2:end)+binz_2(1:end-1))/2;

accuracy_d_binned = nan(nbinz, nbinz);
accuracy_l_binned = nan(nbinz, nbinz);
if type == 1
    
    indi_fa_all = (Cdata_d == 0 & Resp_d' == 1);
    indi_miss_all = (Cdata_d == 1 & Resp_d' == 0);
    indi_hit_all = (Cdata_d == 1 & Resp_d' == 1);
    
    indi_target_pres = Cdata_d == 1;
    indi_target_abs = Cdata_d == 0;
end
for ii = 1:nbinz
    indi_1 = T_NT_sim > binz_1(ii) & T_NT_sim<= binz_1(ii+1);
    for jj = 1: nbinz
        indi_2 = NT_NT_sim > binz_2(jj) & NT_NT_sim<= binz_2(jj+1);
        n_points(ii,jj) = sum(indi_1 & indi_2);
        if type == 1
            accuracy_d_binned(ii,jj) = mean(Cdata_d(indi_1 & indi_2) == Resp_d(indi_1 & indi_2)');
            fa_rate_binned(ii,jj) = sum(indi_1 & indi_2 & indi_fa_all)/sum(indi_target_abs & indi_1 & indi_2);
            miss_rate_binned(ii,jj) = sum(indi_1 & indi_2 & indi_miss_all)/sum(indi_target_pres & indi_1 & indi_2);
            hit_rate_binned(ii,jj) = sum(indi_1 & indi_2 & indi_hit_all)/sum(indi_target_pres & indi_1 & indi_2);
        elseif type == 2
            accuracy_l_binned(ii,jj) = mean(Ldata_l(indi_1 & indi_2) == Resp_l(indi_1 & indi_2)');
        end
    end
end

error_rate_d_binned = 1-accuracy_d_binned;
error_rate_l_binned = 1-accuracy_l_binned;



%%
figure

%for paper
measure = fa_rate_binned; climz = [0 1];
%measure = hit_rate_binned; climz = [0 1];
%measure = 1 - accuracy_d_binned; climz = [0 1];

%measure =  1 - accuracy_l_binned; climz = [0 1];

imagesc(binz_posE_1, binz_posE_2, measure', climz); 
colormap(gray); colorbar
set(gca, 'Ydir', 'normal')
colorbar
box off
set(gca, 'tickdir', 'out')
xlabel('T-D mean (  )')
ylabel('distractor variance')

%psname = 'det_fa_rate.pdf'
%psname = 'det_hit_rate.pdf'
%psname = 'det_error_rate.pdf'

%psname = 'loc_error_rate.pdf'
print_pdf(psname)



%% calc the distributions of distances from circle 3 and circle 4 in Fig 11 from paper

figure;

%circle 3
ii = 1; jj = 1;
indi_1_3 = T_NT_sim > binz_1(ii) & T_NT_sim<= binz_1(ii+1);
indi_2_3 = NT_NT_sim > binz_2(jj) & NT_NT_sim<= binz_2(jj+1);
mask3 = indi_1_3 & indi_2_3;

vec_3 = circ_dist(sdata_d(mask3,:), target_d(mask3)');

target_present_3 = Cdata_d(mask3) == 1;
target_idx_3 = targetidx(mask3);

for k = 1:size(vec_3,1)
    if target_present_3(k)
        vec_3(k,target_idx_3(k)) = NaN;
    end
end

distances_3 = vec_3(:);
distances_3 = distances_3(~isnan(distances_3));

subplot(1,2,1)
hist(distances_3); xlim([-pi pi]); ylim([0 10000])


%circle 4
ii = 1; jj = 7;
indi_1_4 = T_NT_sim > binz_1(ii) & T_NT_sim<= binz_1(ii+1);
indi_2_4 = NT_NT_sim > binz_2(jj) & NT_NT_sim<= binz_2(jj+1);
mask4 = indi_1_4 & indi_2_4;

vec_4 = circ_dist(sdata_d(mask4,:), target_d(mask4)');

target_present_4 = Cdata_d(mask4) == 1;
target_idx_4 = targetidx(mask4);

for k = 1:size(vec_4,1)
    if target_present_4(k)
        vec_4(k,target_idx_4(k)) = NaN;
    end
end

distances_4 = vec_4(:);
distances_4 = distances_4(~isnan(distances_4));

subplot(1,2,2)
hist(distances_4); xlim([-pi pi]); ylim([0 10000])


psname = 'sim_data_for_fig_11_EF.pdf'
print_pdf(psname)