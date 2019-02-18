clear all; close all;

exp_i = 1; % exp 1: stimulus spacing = 60 dva, exp 2: 30 dva
load(['alldata_exp',num2str(exp_i),'.mat'])
mi = 1; % 1: wo decision noise, 2: decision noise
type = 3; %1: detection, 2: localization, 3: joint


model_pred = 1; % flag if we want to also calculate and plot model predictions
model_predR = 0; % for model of specific localization responses, beyond prop corr
Nsamp = 200;

Nsubj = size(alldata,1);
Ncond = size(alldata,2);
Ntrials = length(alldata(1).data.N);
Nvec = unique(alldata(1).data.N)';


nbinz = 6;
nbinz2 = 6;
nbinz3 = 3;
if exp_i == 1
    nbinz_jj = 3;
elseif exp_i == 2
    nbinz_jj = 5;  % for exp 2, distances can go up to 5
end

switch exp_i
    case 1
        loc_vec = [0 2 3 3];
        loc_vec_all = [60 120 180];
        step_loc = 60;
    case 2
        loc_vec = [0 2 3 5];
        loc_vec_all = [30 60 90 120 150];
        step_loc = 30;
end


for Nb = 1:nbinz
    binz_types(Nb,1:(nbinz+1)) = pi*(1-(1-(([1:1:(nbinz+1)]-1)/nbinz)).^(1/Nb));
end

for i = 1:4 %cond
    for j = 1:4 %set size
        if mod(i,2) == 1 %if cond odd
            binz(i,j,1,1:(nbinz+1)) = binz_types(Nvec(j)-1,1:(nbinz+1)); %if target present
            binz(i,j,2,1:(nbinz+1)) = binz_types(Nvec(j),1:(nbinz+1)); %if target present
        else %if cond even
            binz(i,j,1,1:(nbinz+1)) = binz_types(Nvec(j)-1,1:(nbinz+1)); %target always present for loc
            binz(i,j,2,1:(nbinz+1)) = binz_types(Nvec(j)-1,1:(nbinz+1)); %same
        end 
        bincenters(i,j,1,1:nbinz) = (squeeze(binz(i,j,1,2:end))+squeeze(binz(i,j,1,1:end-1)))/2;
        bincenters(i,j,2,1:nbinz) = (squeeze(binz(i,j,2,2:end))+squeeze(binz(i,j,2,1:end-1)))/2; 
    end
end

bincenters = bincenters/2;  %since we doubled all our stimuli for modelling on the (-pi,pi) full space
bincenters = 180*bincenters/pi;


all_min_ori_diff = [];
all_mean_ori_diff = [];
error_rate = nan(Nsubj,Ncond,length(Nvec));


colormat = [77 0 75; 129 15 124; 136 65 157; 140 107 177; ...
    140 150 198; 158 188 218; 191 211 230; 224 236 244; 247 252 253]/256;
eb_w = 0.4;%errorbar width
eb_t = 0.5; %errorbar line thickness

text_labels = {'Perception', 'Memory'};


n_models = 1;

p_c = nan(Nsubj, Ncond, length(Nvec));
preds = nan(Nsubj, n_models, Ncond, Ntrials);

p_resp_ori_diff_rank = nan(Nsubj, Ncond, length(Nvec), nbinz_jj+1);
p_resp_loc_diff_rank = nan(Nsubj, Ncond, length(Nvec), nbinz_jj+1);
p_resp_ori_loc_diff_rank = nan(Nsubj, Ncond, length(Nvec), nbinz_jj+1, nbinz_jj+1);

binzz_meanE = nan(Nsubj, Ncond, length(Nvec), nbinz2+1);
binzz_cstdE = nan(Nsubj, Ncond, length(Nvec), nbinz2+1);


curr_dir = pwd;


for sbjid = 1:Nsubj
    
    sbjid
    
    
    % read in the output of modelling
    if model_pred
        
        model_fits_name = ['/model_pred/']
        
        dirname = [curr_dir,model_fits_name];
        filez = dir([curr_dir,model_fits_name]);
        flst = {filez.name};
        
        cd(dirname)
        
        load(['model_pred_exp_', num2str(exp_i),'_model_',num2str(mi), '_type_', num2str(type),'_sbj_',num2str(sbjid),'.mat'])
       
        predz(mi,1,:) = pred(1,:);
        predz(mi,2,:) = pred(2,:);
        predz(mi,3,:) = pred(3,:);
        predz(mi,4,:) = pred(4,:);
        
        
        if type == 1 % det
            preds(sbjid,mi,1,:) = pred(1,:);
            preds(sbjid,mi,3,:) = pred(3,:);
        elseif type == 2 % loc
            preds(sbjid,mi,2,:) = pred(2,:);
            preds(sbjid,mi,4,:) = pred(4,:);
        elseif type == 3 %both
            preds(sbjid,mi,1,:) = pred(1,:);
            preds(sbjid,mi,3,:) = pred(3,:);
            
            preds(sbjid,mi,2,:) = pred(2,:);
            preds(sbjid,mi,4,:) = pred(4,:);
        end
        
    end
    
    
    for cond = 1:Ncond
        
        data = alldata(sbjid,cond).data;
        
        if mod(cond,2) == 1    % Detection trials
            
            idx_pres = data.target_pres == 1;
            idx_abs  = data.target_pres == 0;
            
            idx_corr  = (data.target_pres == data.response);
            
            if model_pred
                ptp = squeeze(preds(sbjid,mi,cond,:));
                %ptp=Predict_all(params(sbjid,cond,:), {data},1);
                
            end
            
            % characterize prop target present and absent relative to
            % 3 summary statistics related to the distractors distribution
            % 1) min ori diff of the target to the distractors
            % 2) mean ori diff of the distractors
            % 3) circular standard deviation of the distractors
            
            min_ori_diff = nan(Ntrials,1);
            mean_ori_diff = nan(Ntrials,1);
            cstd_ori_diff = nan(Ntrials,1);
            
            % split by set size
            for Nind = 1:length(Nvec)
                N = Nvec(Nind);
                idx_N = data.N == N;
                
                p_tp(sbjid, cond, Nind) = sum(data.response(idx_pres & idx_N)==1)/sum(idx_pres & idx_N);
                p_ta(sbjid, cond, Nind) = sum(data.response(idx_abs & idx_N)==1)/sum(idx_abs & idx_N);
                
                n_tp(sbjid,cond,Nind) = sum(idx_pres & idx_N);
                n_ta(sbjid,cond,Nind) = sum(idx_abs & idx_N);
                
                error_rate(sbjid, cond, Nind) = ((1-p_tp(sbjid, cond, Nind))*sum(idx_pres & idx_N) + p_ta(sbjid, cond, Nind)* sum(idx_abs & idx_N))/ sum(idx_N);
                error_rate_mi(sbjid, cond, Nind) = ((1-p_tp(sbjid, cond, Nind))*sum(idx_pres & idx_N)/ sum( idx_N))/2;
                error_rate_fa(sbjid, cond, Nind) = (p_ta(sbjid, cond, Nind)* sum(idx_abs & idx_N)/ sum( idx_N))/2;
                
                if model_pred
                    p_pr_tp(sbjid, cond, Nind) = mean(ptp(idx_pres & idx_N));
                    p_pr_ta(sbjid, cond, Nind) = mean(ptp(idx_abs & idx_N));
                    error_rate_pred(sbjid, cond, Nind) = ((1-p_pr_tp(sbjid, cond, Nind))*sum(idx_pres & idx_N) + p_pr_ta(sbjid, cond, Nind)* sum(idx_abs & idx_N))/ sum(idx_N);
                end
                
                
                idx_target_abs = find(data.N==N & data.target_pres==0);
                for j = 1:length(idx_target_abs)
                    idx = idx_target_abs(j);
                    vall = [abs(bsxfun(@minus,data.stims(idx,:),data.target_val(idx)))]'; % in the interval (0,2*pi)
                    vall(vall>pi) = 2*pi-vall(vall>pi);
                    min_ori_diff(idx) = min(vall);% in the interval (0,pi)
                    
                    vall_mean = mean(setdiff(data.stims(idx,find(~isnan(data.stims(idx,:)))), data.target_val(idx)));
                    vall_mean = abs(bsxfun(@minus,vall_mean, data.target_val(idx)));
                    vall_mean(vall_mean>pi) = 2*pi- vall_mean(vall_mean>pi);
                    mean_ori_diff(idx) = vall_mean;% in the interval (0,pi)
                    
                    vall_cstd_vec = (setdiff(data.stims(idx,find(~isnan(data.stims(idx,:)))), data.target_val(idx)));
                    vall_cstd_vec(vall_cstd_vec>pi) = 2*pi-vall_cstd_vec(vall_cstd_vec>pi);
                    vall_cstd = circ_var(vall_cstd_vec,[],[],[]);
                    cstd_ori_diff(idx) = vall_cstd; % in the interval (0,1)
                end
                
                
                idx_target_pres = find(data.N ==N & data.target_pres==1);
                for j = 1:length(idx_target_pres)
                    idx = idx_target_pres(j);
                    vall = [abs(bsxfun(@minus,setdiff(data.stims(idx,:),data.target_val(idx)),data.target_val(idx)))]'; %btwn (0,2*pi)
                    vall(vall>pi) = 2*pi-vall(vall>pi);
                    min_ori_diff(idx) = min(vall);
                    
                    vall_mean = mean(setdiff(data.stims(idx,find(~isnan(data.stims(idx,:)))), data.target_val(idx)));
                    vall_mean = abs(bsxfun(@minus,vall_mean, data.target_val(idx)));
                    vall_mean(vall_mean>pi) = 2*pi- vall_mean(vall_mean>pi);
                    mean_ori_diff(idx) = vall_mean;
                    
                    
                    vall_cstd_vec = (setdiff(data.stims(idx,find(~isnan(data.stims(idx,:)))), data.target_val(idx)));
                    vall_cstd_vec(vall_cstd_vec>pi)= 2*pi-vall_cstd_vec(vall_cstd_vec>pi);
                    vall_cstd = circ_var(vall_cstd_vec,[],[],[]);
                    cstd_ori_diff(idx) = vall_cstd;
                end
                
                % split by min ori diff
                for j = 1:nbinz
                    idx_bin_p  = min_ori_diff<binz(cond,Nind,1,j+1) & min_ori_diff>binz(cond,Nind,1,j);
                    idx_bin_a  = min_ori_diff<binz(cond,Nind,2,j+1) & min_ori_diff>binz(cond,Nind,2,j);
                    
                    
                    idx_p = idx_N & idx_bin_p & idx_pres;
                    p_data_tp(sbjid, cond, Nind, j) = sum(data.response(idx_p)==1)/sum(idx_p);
                    n_data_tp(sbjid, cond, Nind, j) = sum(idx_p);
                    
                    idx_a = idx_N & idx_bin_a & idx_abs;
                    p_data_ta(sbjid, cond, Nind, j) = sum(data.response(idx_a)==1)/sum(idx_a);
                    n_data_ta(sbjid, cond, Nind, j) = sum(idx_a);
                    
                    error_rate_all(sbjid, cond, Nind,j)=((1-p_data_tp(sbjid, cond, Nind,j))*sum(idx_p) + p_data_ta(sbjid, cond, Nind,j)* sum(idx_a))/ (sum(idx_p)+sum(idx_a));
                    
                    rt_data_w(sbjid, cond, Nind, j) = median(data.reaction_time(idx_N & idx_bin_p));
                    
                    if model_pred
                        p_pred_tp(sbjid, cond, Nind, j) = mean(ptp(idx_p));
                        p_pred_ta(sbjid, cond, Nind, j) = mean(ptp(idx_a));
                        error_rate_pred_all(sbjid, cond, Nind,j) = ((1-p_pred_tp(sbjid, cond, Nind,j))*sum(idx_p) + p_pred_ta(sbjid, cond, Nind,j)* sum(idx_a))/ (sum(idx_p)+sum(idx_a));
                    end
                    
                    
                end
                
                binzz_mean = [];
                binzz_cstd = [];
                for j = 1: nbinz2
                    binzz_mean(j) = quantile(mean_ori_diff(data.N==N),j/(nbinz2));
                    binzz_cstd(j) = quantile(cstd_ori_diff(data.N==N),j/(nbinz2));
                end
                binzz_mean = [min(mean_ori_diff(data.N==N)) binzz_mean];
                binzz_cstd = [min(cstd_ori_diff(data.N==N)) binzz_cstd];
                
                binzz_meanE(sbjid,cond,Nind,:) = binzz_mean;
                binzz_cstdE(sbjid,cond,Nind,:) = binzz_cstd;
                
                %split by mean of distractors and cstd of distractors
                for j = 1:nbinz2
                    
                    idx_bin_mean  = mean_ori_diff<binzz_mean(j+1) & mean_ori_diff>binzz_mean(j);
                    idx_bin_cstd  = cstd_ori_diff<binzz_cstd(j+1) & cstd_ori_diff>binzz_cstd(j);
                    
                    idxb_mean_p = idx_pres & idx_N & idx_bin_mean;
                    idxb_mean_a = idx_abs & idx_N & idx_bin_mean;
                    
                    idxb_cstd_p = idx_pres & idx_N & idx_bin_cstd;
                    idxb_cstd_a = idx_abs & idx_N & idx_bin_cstd;
                    
                    p_data_tp_mean(sbjid, cond, Nind, j) = sum(data.response(idxb_mean_p)==1)/sum(idxb_mean_p);
                    p_data_ta_mean(sbjid, cond, Nind, j) = sum(data.response(idxb_mean_a)==1)/sum(idxb_mean_a);
                    
                    p_data_tp_cstd(sbjid, cond, Nind, j) = sum(data.response(idxb_cstd_p)==1)/sum(idxb_cstd_p);
                    p_data_ta_cstd(sbjid, cond, Nind, j) = sum(data.response(idxb_cstd_a)==1)/sum(idxb_cstd_a);
                    
                    n_data_w_mean(sbjid, cond, Nind, j) = sum(idx_N & idx_bin_mean);
                    n_data_w_cstd(sbjid, cond, Nind, j) = sum(idx_N & idx_bin_cstd);
                    
                    error_rate_mean_all(sbjid, cond, Nind,j) = ((1-p_data_tp_mean(sbjid, cond, Nind,j))*sum(idxb_mean_p) + p_data_ta_mean(sbjid, cond, Nind,j)* sum(idxb_mean_a))/ (sum(idxb_mean_p)+sum(idxb_mean_a));
                    error_rate_cstd_all(sbjid, cond, Nind,j) = ((1-p_data_tp_cstd(sbjid, cond, Nind,j))*sum(idxb_cstd_p) + p_data_ta_cstd(sbjid, cond, Nind,j)* sum(idxb_cstd_a))/ (sum(idxb_cstd_p)+sum(idxb_cstd_a));
                    
                    
                    if model_pred
                        p_pred_tp_mean(sbjid, cond, Nind, j) = mean(ptp(idxb_mean_p));
                        p_pred_ta_mean(sbjid, cond, Nind, j) = mean(ptp(idxb_mean_a));
                        
                        error_rate_mean_pred_all(sbjid, cond, Nind,j)=((1-p_pred_tp_mean(sbjid, cond, Nind,j))*sum(idxb_mean_p) + p_pred_ta_mean(sbjid, cond, Nind,j)* sum(idxb_mean_a))/ (sum(idxb_mean_p)+sum(idxb_mean_a));
                        
                        p_pred_tp_cstd(sbjid, cond, Nind, j) = mean(ptp(idxb_cstd_p));
                        p_pred_ta_cstd(sbjid, cond, Nind, j) = mean(ptp(idxb_cstd_a));
                        
                        error_rate_cstd_pred_all(sbjid, cond, Nind,j) = ((1-p_pred_tp_cstd(sbjid, cond, Nind,j))*sum(idxb_cstd_p) + p_pred_ta_cstd(sbjid, cond, Nind,j)* sum(idxb_cstd_a))/ (sum(idxb_cstd_p)+sum(idxb_cstd_a));
                    end
                    
                end
                
                
                
            end
            
            
            all_min_ori_diff = [all_min_ori_diff min_ori_diff];
        else %cond even, localization data
            
            if model_pred
                pc = squeeze(preds(sbjid,mi,cond,:));
                %pc=Predict_all(params(sbjid,cond,:), {data},2);
                if model_predR
                    [pcr ddr] = Predict_allR(params(sbjid, mi,cond/2, 1:5), 0, {data}, type);
                end
            end
            
            idx_corr  = (data.target_loc == data.response) ;
            idx_incorr  = (data.target_loc ~= data.response) ;
            
            
            min_ori_diff = nan(Ntrials,1);
            mean_ori_diff = nan(Ntrials,1);
            cstd_ori_diff = nan(Ntrials,1);
            
            
            loc_diff = nan(Ntrials,1);
            resp_ori_diff_rank = nan(Ntrials,1);
            resp_loc_diff_rank = nan(Ntrials,1);
            
            for Nind = 1:length(Nvec)
                
                N = Nvec(Nind);
                idx_N = data.N == N;
                
                p_c(sbjid, cond, Nind) = sum(idx_corr & idx_N)/sum(idx_N);
                %rt_c(sbjid, cond, Nind) = median(data.reaction_time(idx_N));
                %rt_c_corr(sbjid, cond, Nind) = median(data.reaction_time(idx_corr & idx_N));
                %rt_c_incorr(sbjid, cond, Nind) = median(data.reaction_time(idx_incorr & idx_N));
                
                n_c(sbjid,cond,Nind) = sum(idx_N);
                
                error_rate(sbjid,cond,Nind) = 1-p_c(sbjid, cond, Nind);
                if model_pred
                    p_pr_c(sbjid, cond, Nind) = mean(pc(idx_corr & idx_N));
                end
                
                
                idx_N_val = find(data.N==N);
                
                for j=1:length(idx_N_val)
                    idx=idx_N_val(j);
                    vall = [abs(bsxfun(@minus,data.stims(idx,:),data.target_val(idx)))]'; %btwn (0,2*pi)
                    vall(vall>pi) = 2*pi-vall(vall>pi);
                    min_ori_diff(idx) = min(vall(vall>0)); %min(vall); %btwn (0,pi)
                    
                    vall2 = abs(data.stims(idx,data.response(idx)) - data.target_val(idx));
                    if vall2 > pi
                        vall2 = 2*pi-vall2;
                    end
                    vall_sorted = sort(vall);
                    
                    vall_mean = mean(setdiff(data.stims(idx,find(~isnan(data.stims(idx,:)))), data.target_val(idx)));
                    vall_mean = abs(bsxfun(@minus,vall_mean, data.target_val(idx)));
                    vall_mean(vall_mean>pi) = 2*pi- vall_mean(vall_mean>pi);
                    mean_ori_diff(idx) = vall_mean;
                    
                    
                    vall_cstd_vec = (setdiff(data.stims(idx,find(~isnan(data.stims(idx,:)))), data.target_val(idx)));
                    vall_cstd_vec(vall_cstd_vec>pi) = 2*pi-vall_cstd_vec(vall_cstd_vec>pi);
                    vall_cstd = circ_var(vall_cstd_vec,[],[],[]);
                    cstd_ori_diff(idx) = vall_cstd;
                    
                    if model_predR
                        for simi = 1: Nsamp
                            data_stims_idx = data.stims(idx,~isnan(data.stims(idx,:)));
                            vall3 = abs(data_stims_idx(ddr(idx,simi)) - data.target_val(idx));
                            if vall3 > pi
                                vall3 = 2*pi-vall3;
                            end
                            resp_ori_diff_rank_simi(idx,simi) = find(vall3==vall_sorted)-1;
                            resp_loc_diff_rank_simi(idx,simi) = abs(find(vall==0) - find(vall==vall3));
                            
                            if exp_i == 1
                                
                                if resp_ori_diff_rank_simi(idx, simi) >= 4
                                    resp_ori_diff_rank_simi(idx, simi) = 6-resp_ori_diff_rank_simi(idx, simi);
                                end
                                if resp_loc_diff_rank_simi(idx, simi) >= 4
                                    resp_loc_diff_rank_simi(idx, simi) = 6-resp_loc_diff_rank_simi(idx, simi);
                                end
                                
                            elseif exp_i == 2
                                
                                if resp_ori_diff_rank_simi(idx, simi) >= 7
                                    resp_ori_diff_rank_simi(idx, simi) = 12-resp_ori_diff_rank_simi(idx, simi);
                                end
                                if resp_loc_diff_rank_simi(idx, simi) >= 7
                                    resp_loc_diff_rank_simi(idx, simi) = 12-resp_loc_diff_rank_simi(idx, simi);
                                end
                                
                            end
                            
                        end
                    end
                    
                    resp_ori_diff_rank(idx) = find(vall2==vall_sorted)-1;
                    resp_loc_diff_rank(idx) = abs(find(vall==0) - find(vall==vall2));
                    
                    min_loc=find(vall == min_ori_diff(idx));
                    loc_diff(idx) = abs(min_loc-data.target_loc(idx));
                    
                    if exp_i == 1
                        
                        if loc_diff(idx) >= 4
                            loc_diff(idx) = 6-loc_diff(idx);
                        end
                        if resp_ori_diff_rank(idx) >= 4
                            resp_ori_diff_rank(idx) = 6-resp_ori_diff_rank(idx);
                        end
                        if resp_loc_diff_rank(idx) >= 4
                            resp_loc_diff_rank(idx) = 6-resp_loc_diff_rank(idx);
                        end
                        
                    elseif exp_i == 2
                        
                        if loc_diff(idx) >= 7
                            loc_diff(idx) = 12 - loc_diff(idx);
                        end
                        if resp_ori_diff_rank(idx) >= 7
                            resp_ori_diff_rank(idx) = 12-resp_ori_diff_rank(idx);
                        end
                        if resp_loc_diff_rank(idx) >= 7
                            resp_loc_diff_rank(idx) = 12-resp_loc_diff_rank(idx);
                        end
                        
                    end
                    
                end
                
                %location only, not also noise
                for jj = 1:nbinz_jj
                    idx_loci = [(loc_diff==jj)]; %1,2,3 or 1,2,3,4,5
                    idxl = idx_corr & idx_N  & idx_loci;
                    idxl_inc = idx_incorr & idx_N & idx_loci;
                    p_data_l(sbjid, cond, Nind, jj) = sum(idxl)/sum(idx_N &  idx_loci);
                    n_data_l(sbjid, cond, Nind, jj) = sum(idx_N & idx_loci);
                    if model_pred
                        p_pred_l(sbjid, cond, Nind, jj) = mean(pc(idxl));
                    end

                
                    %rt_data_l(sbjid,cond, Nind, jj) = median(data.reaction_time(idx_N  & idx_loci));
                    %rt_data_l_corr(sbjid,cond, Nind, jj) = median(data.reaction_time(idxl));%median(data.reaction_time(idx));
                    %rt_data_l_incorr(sbjid, cond, Nind, jj) = median(data.reaction_time(idxl_inc));%median(data.reaction_time(idx_inc));
                    
                    if jj <= loc_vec(Nind)
                        idx_loco = [(resp_ori_diff_rank==jj)];
                        idx_lo = idx_loco & idx_N;
                        p_resp_ori_diff_rank(sbjid, cond, Nind, jj) = sum(idx_lo)/length(idx_N_val);
                        
                        
                        idx_locl = [(resp_loc_diff_rank==jj)];
                        idx_ll = idx_locl & idx_N;
                        p_resp_loc_diff_rank(sbjid, cond, Nind,jj) = sum(idx_ll)/length(idx_N_val);
                    end
                    
                    if model_predR
                        for simi = 1:Nsamp
                            if jj <= loc_vec(Nind)
                                idx_loco = [(resp_ori_diff_rank_simi(:,simi)==jj)];
                                idx_lo = idx_loco & idx_N;
                                p_resp_ori_diff_rank_simi(sbjid, cond, Nind, jj, simi) = sum(idx_lo)/length(idx_N_val);
                                
                                
                                idx_locl = [(resp_loc_diff_rank_simi(:,simi)==jj)];
                                idx_ll = idx_locl & idx_N;
                                p_resp_loc_diff_rank_simi(sbjid, cond, Nind,jj, simi) = sum(idx_ll)/length(idx_N_val);
                            end
                        end
                    end
                    
                end
                
                if Nind>1
                    for iii = 1:loc_vec(Nind)
                        for jjj = 1:loc_vec(Nind)
                            idx_locoi = [(resp_ori_diff_rank==iii)];
                            idx_loclj = [(resp_loc_diff_rank==jjj)];
                            idx_lol = idx_locoi & idx_loclj & idx_N;
                            p_resp_ori_loc_diff_rank(sbjid, cond, Nind, iii, jjj) = sum(idx_lol)/length(idx_N_val);
                        end
                    end
                end
                
                for j = 1:nbinz
                    idx_bin  = min_ori_diff<binz(cond,Nind,1,j+1) & min_ori_diff>binz(cond,Nind,1,j);
                    idxb = idx_corr & idx_N & idx_bin;
                    idx_bin_rt = idx_bin;
                    
                    p_data_w(sbjid, cond, Nind, j) = sum(idxb)/sum(idx_N & idx_bin);
                    n_data_w(sbjid, cond, Nind, j) = sum(idx_N & idx_bin);
                    
                    rt_data_w(sbjid, cond, Nind, j) = median(data.reaction_time(idx_N & idx_bin));
                    %rt_data_w_corr(sbjid, cond, Nind, j) = median(data.reaction_time(idx_corr & idx_N & idx_bin));
                    %rt_data_w_incorr(sbjid, cond, Nind, j) = median(data.reaction_time(idx_incorr & idx_N & idx_bin));
                    
                    if model_pred
                        p_pred_w(sbjid, cond, Nind, j) = mean(pc(idxb));
                    end
                    
                    for jj = 1:nbinz_jj % noise and loc
                        idx_loci = [(loc_diff==jj)];
                        idxx = idx_corr & idx_N & idx_bin & idx_loci;
                        idxx_inc = idx_incorr & idx_N & idx_bin & idx_loci;
                        
                        p_data_wl(sbjid, cond, Nind, j, jj) = sum(idxx)/sum(idx_N & idx_bin & idx_loci);
                        n_data_wl(sbjid, cond, Nind, j, jj) = sum(idx_N & idx_bin & idx_loci);
                        
                        %rt_data_wl(sbjid,cond, Nind, j, jj) =  median(data.reaction_time(idx_N & idx_bin & idx_loci));
                        %rt_data_wl_corr(sbjid,cond, Nind, j, jj) = median(data.reaction_time(idxx));%median(data.reaction_time(idx));
                        %rt_data_wl_incorr(sbjid, cond, Nind, j, jj) = median(data.reaction_time(idxx_inc));%median(data.reaction_time(idx_inc));
                    end
                    
                end
                
                binzz_mean = [];
                binzz_cstd = [];
                for j = 1: nbinz2
                    binzz_mean(j) = quantile(mean_ori_diff(data.N==N),j/(nbinz2));
                    binzz_cstd(j) = quantile(cstd_ori_diff(data.N==N),j/(nbinz2));
                end
                binzz_mean = [min(mean_ori_diff(data.N==N)) binzz_mean];
                binzz_cstd = [min(cstd_ori_diff(data.N==N)) binzz_cstd];
                
                binzz_meanE(sbjid,cond,Nind,:) = binzz_mean;
                binzz_cstdE(sbjid,cond,Nind,:) = binzz_cstd;
                for j = 1:nbinz2
                    idx_bin_mean  = mean_ori_diff<binzz_mean(j+1) & mean_ori_diff>binzz_mean(j);
                    idx_bin_cstd  = cstd_ori_diff<binzz_cstd(j+1) & cstd_ori_diff>binzz_cstd(j);
                    
                    idxb_mean = idx_corr & idx_N & idx_bin_mean;
                    idxb_cstd = idx_corr & idx_N & idx_bin_cstd;
                    
                    p_data_w_mean(sbjid, cond, Nind, j) = sum(idxb_mean)/sum(idx_N & idx_bin_mean);
                    p_data_w_cstd(sbjid, cond, Nind, j) = sum(idxb_cstd)/sum(idx_N & idx_bin_cstd);
                    
                    n_data_w_mean(sbjid, cond, Nind, j) = sum(idx_N & idx_bin_mean);
                    n_data_w_cstd(sbjid, cond, Nind, j) = sum(idx_N & idx_bin_cstd);
                    if model_pred
                        %p_pred_w(sbjid, cond, Nind, j) = mean(pc(idxb));
                        p_pred_w_mean(sbjid, cond, Nind, j) = mean(pc(idxb_mean));
                        p_pred_w_cstd(sbjid, cond, Nind, j) = mean(pc(idxb_cstd));
                    end
                end
                
            end
            
        end
        
        
        all_min_ori_diff = [all_min_ori_diff min_ori_diff];
        all_mean_ori_diff = [all_mean_ori_diff mean_ori_diff];
        
        %end
    end
end


if model_pred
    cd ..
    close all;
end

%%
savefilename_pars_nll = ['params_nll_model_',num2str(mi),'_type_',num2str(type) ,'_exp_',num2str(exp_i), '.mat']
%save(savefilename_pars_nll,'params','nll_sbj','-mat')
%% comprehensive figure


for i=1:4 %cond
    for j=1:4 %set size
        if mod(i,2)==1 %if cond odd
            binz(i,j,1,1:(nbinz+1)) = binz_types(Nvec(j)-1,1:(nbinz+1)); %if target present
            binz(i,j,2,1:(nbinz+1)) = binz_types(Nvec(j),1:(nbinz+1)); %if target present
        else %if cond even
            binz(i,j,1,1:(nbinz+1)) = binz_types(Nvec(j)-1,1:(nbinz+1)); %target always present for loc
            binz(i,j,2,1:(nbinz+1)) = binz_types(Nvec(j)-1,1:(nbinz+1)); %same
        end
        
        bincenters(i,j,1,1:nbinz) = (squeeze(binz(i,j,1,2:end))+squeeze(binz(i,j,1,1:end-1)))/2;
        bincenters(i,j,2,1:nbinz) = (squeeze(binz(i,j,2,2:end))+squeeze(binz(i,j,2,1:end-1)))/2;
        
    end
end
bincenters = bincenters/2;  %since we doubled all our stimuli for modelling on the (-pi,pi) full space
bincenters = 180*bincenters/pi;

binzz_meanE = binzz_meanE/2;  %since we doubled all our stimuli for modelling on the (-pi,pi) full space
binzz_meanE = 180* binzz_meanE/pi;


%%

bincenterz = bincenters;

fontsz=10;
sz_dot=1.2;

color1=[58 67 192]/255;  %blue present
color2=[149 16 16]/255;  %red absent

color_gray=[210 210 210]/255;
color1_shade=(color1+color_gray*1.5)/2;
color2_shade=(color2+color_gray*1.5)/2;

color3=[57 63 0]/255;
color3_shade=(color3+color_gray*1.5)/2;


accuracy=bsxfun(@minus,1,error_rate);
accuracy_all = bsxfun(@minus,1,error_rate_all);
accuracy_mean_all = bsxfun(@minus,1,error_rate_mean_all);
accuracy_cstd_all = bsxfun(@minus,1,error_rate_cstd_all);

if model_pred
    accuracy_pred = bsxfun(@minus,1,error_rate_pred);
    accuracy_pred_all = bsxfun(@minus,1,error_rate_pred_all);
    accuracy_pred_mean_all = bsxfun(@minus,1, error_rate_mean_pred_all);
    accuracy_pred_cstd_all = bsxfun(@minus,1, error_rate_cstd_pred_all);
end

colors_sz = [217 30 0 ;  51 164 0;  242 135 65; 130 44 169]./255;
% red, green, orange, purple
colors_sz_nuances(:,1,1:3) = [142 164 73; 51 164 0; 14 164 108]; % green nuances
colors_sz_nuances(:,2,1:3) = [242 144 118; 242 135 65; 242 34 26]; % orange nuances
colors_sz_nuances(:,3,1:3) = [119 97 169; 130 44 169; 169 10 132]; % purple nuances
colors_sz_nuances = colors_sz_nuances./255;

col_gray = [168 168 168]/255;
color_gray=[210 210 210]/255;
colors_sz_shade =(colors_sz+1.5*repmat(color_gray,4,1))/2;
colors_sz_shade(colors_sz_shade>0.99)= 0.99;

redd = [0.9047    0.1918    0.1988];
bluee = [0.2941    0.5447    0.7494];
redd_shade=(redd+2*[1 1 1])/3;
bluee_shade=(bluee+2*[1 1 1])/3;



guttera=[0.14 0.06];
marginsa=[0.14 0.14 0.04 0.04]; %left right bottom top



%% Detection % if type ==1 or 3


close all;
full_figure = 1;

figure
if full_figure
    set(gcf, 'Position', [100 100 500 920])
    guttera2 = guttera;
    marginsa2 = marginsa;
else
    set(gcf, 'Position', [100 100 500 190])
    guttera2 = guttera;
    marginsa2 = [0.14 0.14 0.2 0.2];
end




for cond = [1 3]
    if full_figure
        
        % a)   set size
        tight_subplot(5,2,1,0.5+cond/2, guttera2, marginsa2)
        if model_pred
            h_p=fill([Nvec Nvec(end:-1:1)], [[squeeze(mean(p_pr_tp(:,cond,:),1))]'-[squeeze(std(p_pr_tp(:,cond,:),1))]'/sqrt(Nsubj)...
                [squeeze(mean(p_pr_tp(:,cond,end:-1:1),1))]'+[squeeze(std(p_pr_tp(:,cond,:,end:-1:1),1))]'/sqrt(Nsubj)],bluee_shade,'EdgeColor', 'None'); hold on;
            h_a=fill([Nvec Nvec(end:-1:1)], [[squeeze(mean(p_pr_ta(:,cond,:),1))]'-[squeeze(std(p_pr_ta(:,cond,:),1))]'/sqrt(Nsubj)...
                [squeeze(mean(p_pr_ta(:,cond,end:-1:1),1))]'+[squeeze(std(p_pr_ta(:,cond,:,end:-1:1),1))]'/sqrt(Nsubj)],redd_shade ,'EdgeColor', 'None'); hold on;
      
        end
        
        plot(Nvec, squeeze(mean(p_tp(:, cond, :),1)), 'Color', bluee, 'Linestyle', 'none'); hold on;
        errorbar(Nvec, squeeze(mean(p_tp(:, cond, :),1)),squeeze(std(p_tp(:, cond, :),1))/sqrt(Nsubj), 'o','MarkerSize', sz_dot, 'MarkerFaceColor', bluee, 'MarkerEdgeColor', bluee,'Color',bluee , 'LineWidth',1); hold on;
        plot(Nvec, squeeze(mean(p_ta(:, cond, :),1)), 'Color', redd, 'Linestyle', 'none'); hold on;
        errorbar(Nvec, squeeze(mean(p_ta(:, cond, :),1)),squeeze(std(p_ta(:, cond, :),1))/sqrt(Nsubj), 'o','MarkerSize', sz_dot, 'MarkerFaceColor', redd, 'MarkerEdgeColor',redd,'Color',redd, 'LineWidth',1); hold on;
        chance=plot(Nvec, 0.5*ones(1, length(Nvec)),'--k'); hold on;
        box off
        ylim([0.3 1])
        xlim([min(Nvec)-0.5 max(Nvec)+0.5]);
        box off
        set(gca, 'tickdir', 'out')
        set(gca, 'xtick', Nvec)
        set(gca, 'xticklabels', {'2', '3', '4', '6'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
        set(gca, 'ytick', [0:0.1:1])
        set(gca, 'yticklabels', {})
        set(gca,'ticklength', [0.04 0.04])
        xlabel('Set size','FontName', 'Helvetica', 'FontSize', fontsz)
        if cond == 1
            ylabel('Proportion target present','FontName', 'Helvetica', 'FontSize', fontsz)
            set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
        end
        text(3, 1.16, text_labels{(1+cond)/2},  'FontName', 'Helvetica', 'FontSize', fontsz)
        
        
        %b) T-most similar of distractors
        
        for Nind = 1:4
            bincentersl=[squeeze(bincenterz(cond, Nind,1,:))]';
            
            tight_subplot(5,2,2,0.5+cond/2, guttera2, marginsa2)
            
            if model_pred
                hacc=fill([bincentersl bincentersl(end:-1:1)], [[squeeze(mean(accuracy_pred_all(:,cond,Nind,:),1))]'-[squeeze(std(accuracy_pred_all(:,cond,Nind,:),1))./sqrt(squeeze(sum(~isnan(accuracy_pred_all(:,cond,Nind,:)),1)))]'...
                    [squeeze(mean(accuracy_pred_all(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(accuracy_pred_all(:,cond,Nind,end:-1:1),1))./sqrt(squeeze(sum(~isnan(accuracy_pred_all(:,cond,Nind,:)),1)))]'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
                
            end
        end
        for Nind = 1:4
            bincentersl=[squeeze(bincenterz(cond, Nind,1,:))]';
            
            tight_subplot(5,2,2,0.5+cond/2, guttera2, marginsa2)
            
            plot(bincentersl, squeeze(mean(accuracy_all(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
            leg(Nind)=errorbar(bincentersl,squeeze(mean(accuracy_all(:,cond,Nind,:),1)), squeeze(std(accuracy_all(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none');
            %plot(bincentersl,squeeze(mean(p_data_tp(:, cond, Nind, :),1)), '-o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',3*sz_dot); hold on;
            %plot(bincentersl,squeeze(mean(p_data_ta(:, cond, Nind, :),1)), '--o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',3*sz_dot); hold on;
            chance = plot(0:10:90, 0.5*ones(1, 10),'--k'); hold on;
            box off
            ylim([0.4 1])
            xlim([0 90])
            set(gca, 'xtick', [0:15:90])
            set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', fontsz)
            set(gca, 'ytick', [0:0.1:1])
            set(gca, 'yticklabels', {})
            if cond == 1
                set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
                ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', [0.04 0.04])
            xlabel('T - MSD orientation distance (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
            
            
            if cond == 4
                llpa = legend([leg(1) leg(2) leg(3) leg(4) chance],'N = 2', 'N = 3','N = 4', 'N = 6','chance','FontName','Helvetica','FontSize', fontsz)
                legend boxoff
                set(llpa, 'FontName','Helvetica','FontSize', fontsz*0.9, 'Position', [0.9 0.84 0.04 0.04])
            end
            
            
            
            
            %c) T-mean of distractors
            for Nind = 2:4
                
                binzz_mean_mean = [mean(squeeze(binzz_meanE(:,cond,Nind,:)),1)]';
                bincenterz_mean = [(binzz_mean_mean(1:end-1)+binzz_mean_mean(2:end))/2]';
                
                tight_subplot(5,2,3,0.5+cond/2, guttera2, marginsa2)
                
                if model_pred
                    h_m(Nind)=fill([bincenterz_mean bincenterz_mean(end:-1:1)], [[squeeze(mean(accuracy_pred_mean_all(:,cond,Nind,:),1))]'-[squeeze(std(accuracy_pred_mean_all(:,cond,Nind,:),1))]'./sqrt(squeeze(sum(~isnan(accuracy_pred_mean_all(:,cond,Nind,:)),1)))'...
                        [squeeze(mean(accuracy_pred_mean_all(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(accuracy_pred_mean_all(:,cond,Nind,end:-1:1),1))]'./sqrt(squeeze(sum(~isnan(accuracy_pred_mean_all(:,cond,Nind,:)),1)))'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
                
                end
            end
            
            for Nind = 2:4
                
                binzz_mean_mean = [mean(squeeze(binzz_meanE(:,cond,Nind,:)),1)]';
                bincenterz_mean = [(binzz_mean_mean(1:end-1)+binzz_mean_mean(2:end))/2]';
                
                tight_subplot(5,2,3,0.5+cond/2, guttera2, marginsa2)
                plot(bincenterz_mean, squeeze(mean(accuracy_mean_all(:,cond,Nind,:),1)),'Color', colors_sz(Nind,:),'Linestyle', 'none','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
                %plot(bincenterz_mean,squeeze(mean(p_data_tp_mean(:, cond, Nind, :),1)), '-o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',3*sz_dot); hold on;
                %plot(bincenterz_mean,squeeze(mean(p_data_ta_mean(:, cond, Nind, :),1)), '--o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',3*sz_dot); hold on;
                leg(Nind) = errorbar(bincenterz_mean,squeeze(mean(accuracy_mean_all(:,cond,Nind,:),1)), squeeze(std(accuracy_mean_all(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none');
                chance = plot(0:10:90, 0.5*ones(1, 10),'--k'); hold on;
                box off
                ylim([0.4 1])
                xlim([0 90])
                set(gca, 'xtick', [0:15:90])
                set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', fontsz)
                set(gca, 'ytick', [0:0.1:1])
                set(gca, 'yticklabels', {})
                set(gca, 'tickdir', 'out')
                set(gca,'ticklength', [0.04 0.04])
                xlabel('T-mean of distractors distance (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
                
                if cond == 1
                    set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
                    ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
                end
                
            end
            
            
            
            % d) circ var of distractors
            for Nind = 2:4
                
                binzz_cstd_mean = [mean(squeeze(binzz_cstdE(:,cond,Nind,:)),1)]';
                bincenterz_cstd = [(binzz_cstd_mean(1:end-1)+binzz_cstd_mean(2:end))/2]';
                
                tight_subplot(5,2,4,0.5+cond/2, guttera2, marginsa2)
                
                if model_pred
                    h_v(Nind)=fill([bincenterz_cstd bincenterz_cstd(end:-1:1)], [[squeeze(nanmean(accuracy_pred_cstd_all(:,cond,Nind,:),1))]'-[squeeze(nanstd(accuracy_pred_cstd_all(:,cond,Nind,:),1))]'./sqrt(squeeze(sum(~isnan(accuracy_pred_cstd_all(:,cond,Nind,:)),1)))'...
                        [squeeze(nanmean(accuracy_pred_cstd_all(:,cond,Nind,end:-1:1),1))]'+[squeeze(nanstd(accuracy_pred_cstd_all(:,cond,Nind,end:-1:1),1))]'./sqrt(squeeze(sum(~isnan(accuracy_pred_cstd_all(:,cond,Nind,:)),1)))'],colors_sz_shade(Nind,:),'EdgeColor', 'None', 'Linestyle', 'none'); hold on;
                
                end
                
            end
            for Nind = 2:4
                
                binzz_cstd_mean = [mean(squeeze(binzz_cstdE(:,cond,Nind,:)),1)]';
                bincenterz_cstd = [(binzz_cstd_mean(1:end-1)+binzz_cstd_mean(2:end))/2]';
                
                tight_subplot(5,2,4,0.5+cond/2, guttera2, marginsa2)
                plot(bincenterz_cstd, squeeze(nanmean(accuracy_cstd_all(:,cond,Nind,:),1)),  'o', 'Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
                leg(Nind) = errorbar(bincenterz_cstd,squeeze(nanmean(accuracy_cstd_all(:,cond,Nind,:),1)), squeeze(nanstd(accuracy_cstd_all(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:),'Linestyle', 'none');
                chance = plot(0:0.10:0.90, ones(1,10)*1/Nvec(Nind),'--', 'Color',colors_sz(Nind,:)); hold on;
                box off
                ylim([0.4 1])
                xlim([0 1])
                set(gca, 'xtick',[0:0.2:1])
                set(gca, 'ytick', [0:0.1:1])
                set(gca, 'yticklabels', {})
                set(gca, 'tickdir', 'out')
                set(gca,'ticklength', [0.04 0.04])
                xlabel('circular variance of distractors','FontName', 'Helvetica', 'FontSize', fontsz)
                if cond == 1
                    set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
                    ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
                end
            end
            
            
            
            
        end
    else
        %b) T-most similar of distractors
        
        for Nind = 1:4
            bincentersl=[squeeze(bincenterz(cond, Nind,1,:))]';
            
            tight_subplot(1,2,1,0.5+cond/2, guttera2, marginsa2)
            
            if model_pred
                hacc=fill([bincentersl bincentersl(end:-1:1)], [[squeeze(mean(accuracy_pred_all(:,cond,Nind,:),1))]'-[squeeze(std(accuracy_pred_all(:,cond,Nind,:),1))./sqrt(squeeze(sum(~isnan(accuracy_pred_all(:,cond,Nind,:)),1)))]'...
                    [squeeze(mean(accuracy_pred_all(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(accuracy_pred_all(:,cond,Nind,end:-1:1),1))./sqrt(squeeze(sum(~isnan(accuracy_pred_all(:,cond,Nind,:)),1)))]'],colors_sz_shade(Nind,:),'EdgeColor', 'None', 'Linestyle', 'none'); hold on;
            
            end
        end
        for Nind = 1:4
            bincentersl=[squeeze(bincenterz(cond, Nind,1,:))]';
            
            tight_subplot(1,2,1,0.5+cond/2, guttera2, marginsa2)
            
            plot(bincentersl, squeeze(mean(accuracy_all(:,cond,Nind,:),1)), 'o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
            leg(Nind)=errorbar(bincentersl,squeeze(mean(accuracy_all(:,cond,Nind,:),1)), squeeze(std(accuracy_all(:,cond,Nind,:),1))/sqrt(Nsubj),'Color', colors_sz(Nind,:),'Linestyle', 'none');
            chance = plot(0:10:90, 0.5*ones(1, 10),'--k'); hold on;
            box off
            ylim([0.4 1])
            xlim([0 90])
            set(gca, 'xtick', [0:15:90])
            set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', fontsz)
            set(gca, 'ytick', [0:0.1:1])
            set(gca, 'yticklabels', {})
            if cond == 1
                set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
                ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', [0.04 0.04])
            xlabel('T - MSD orientation distance (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
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
%%
psname = ['DET_full_figure_',num2str(full_figure),'_model_',num2str(mi),'_type_',num2str(type) ,'_exp_',num2str(exp_i), '.pdf']
%print_pdf(psname)



%% Localization  %% if type = 2 or 3
close all;


full_figure = 1;

figure
if full_figure
    set(gcf, 'Position', [100 100 500 920])
    guttera2 = guttera;
    marginsa2 = marginsa;
else
    set(gcf, 'Position', [100 100 500 190])
    guttera2 = guttera;
    marginsa2 = [0.14 0.14 0.2 0.2];
end


if model_predR
    p_resp_ori_diff_rank_simim = squeeze(mean(p_resp_ori_diff_rank_simi,5));
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


for cond = [2 4]
    
    if full_figure
        %a) with set size
        tight_subplot(5,2,1,cond/2, guttera2, marginsa2)
        if model_pred
            h_p=fill([Nvec Nvec(end:-1:1)], [[squeeze(mean(p_pr_c(:,cond,:),1))]'-[squeeze(std(p_pr_c(:,cond,:),1))]'/sqrt(Nsubj)...
                [squeeze(mean(p_pr_c(:,cond,end:-1:1),1))]'+[squeeze(std(p_pr_c(:,cond,:,end:-1:1),1))]'/sqrt(Nsubj)],col_gray,'EdgeColor', 'None'); hold on;
        end
        plot(Nvec, squeeze(mean(accuracy(:, cond, :),1)), 'Color', 'k', 'Linestyle', 'none'); hold on;
        errorbar(Nvec, squeeze(mean(accuracy(:, cond, :),1)),squeeze(std(accuracy(:, cond, :),1))/sqrt(Nsubj), 'o','MarkerSize', sz_dot, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','Color', 'k', 'LineWidth',1, 'Linestyle', 'none'); hold on;
        chance = plot(Nvec, 1./Nvec,'--k'); hold on;
        
        box off
        ylim([0 1])
        xlim([min(Nvec)-0.5 max(Nvec)+0.5]);
        box off
        set(gca, 'tickdir', 'out')
        set(gca, 'xtick', Nvec)
        set(gca, 'xticklabels', {'2','3', '4', '6'},'FontName', 'Helvetica', 'FontSize', fontsz)
        set(gca, 'ytick', [0:0.1:1])
        set(gca, 'yticklabels',{})
        set(gca, 'tickdir', 'out')
        set(gca,'ticklength', [0.04 0.04])
        xlabel('Set size','FontName', 'Helvetica', 'FontSize', fontsz)
        
        if cond == 2
            set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
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
            plot(bincentersl, squeeze(mean(p_data_w(:,cond,Nind,:),1)), 'o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
            leg(Nind) = errorbar(bincentersl,squeeze(mean(p_data_w(:,cond,Nind,:),1)), squeeze(std(p_data_w(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none');
            chance = plot(0:10:90, ones(1,10)*1/Nvec(Nind),'--', 'Color',colors_sz(Nind,:)); hold on;
            box off
            ylim([0 1])
            xlim([0 90])
            set(gca, 'xtick', [0:15:90])
            set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', fontsz)
            set(gca, 'ytick', [0:0.1:1])
            set(gca, 'yticklabels', {})
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', [0.04 0.04])
            if cond == 2
                set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
                ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            xlabel('T - MSD orientation distance (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
            
        end
        
        if cond == 4
            llpa = legend([leg(1) leg(2) leg(3) leg(4) chance],'N = 2', 'N = 3','N = 4', 'N = 6','chance','FontName','Helvetica','FontSize', fontsz)
            legend boxoff
            set(llpa, 'FontName','Helvetica','FontSize', fontsz*0.9, 'Position', [0.9 0.84 0.04 0.04])
        end
        

        
        %c) T-mean of distractors
        for Nind = 2:4
            
            binzz_mean_mean = [mean(squeeze(binzz_meanE(:,cond,Nind,:)),1)]';
            bincenterz_mean = [(binzz_mean_mean(1:end-1)+binzz_mean_mean(2:end))/2]';
            
            tight_subplot(5,2,3,cond/2, guttera2, marginsa2)
            if model_pred
                h_p=fill([bincenterz_mean bincenterz_mean(end:-1:1)], [[squeeze(mean(p_pred_w_mean(:,cond,Nind,:),1))]'-[squeeze(std(p_pred_w_mean(:,cond,Nind,:),1))]'./sqrt(squeeze(sum(~isnan(p_pred_w_mean(:,cond,Nind,:)),1)))'...
                    [squeeze(mean(p_pred_w_mean(:,cond,Nind,end:-1:1),1))]'+[squeeze(std(p_pred_w_mean(:,cond,Nind,end:-1:1),1))]'./sqrt(squeeze(sum(~isnan(p_pred_w_mean(:,cond,Nind,:)),1)))'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
            end
            plot(bincenterz_mean, squeeze(mean(p_data_w_mean(:,cond,Nind,:),1)), 'o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
            leg(Nind) = errorbar(bincenterz_mean,squeeze(mean(p_data_w_mean(:,cond,Nind,:),1)), squeeze(std(p_data_w_mean(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none');
            chance = plot(0:10:90, ones(1,10)*1/Nvec(Nind),'--', 'Color',colors_sz(Nind,:)); hold on;
            box off
            ylim([0 1])
            xlim([0 90])
            set(gca, 'xtick', [0:15:90])
            set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', fontsz)
            set(gca, 'ytick', [0:0.1:1])
            set(gca, 'yticklabels', {})
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', [0.04 0.04])
            
            xlabel('T-mean of distractors distance (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
            
            if cond == 2
                set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
                ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            
        end
        
        
        
        % d) circ var of distractors
        for Nind = 2:4
            
            binzz_cstd_mean = [mean(squeeze(binzz_cstdE(:,cond,Nind,:)),1)]';
           
            bincenterz_cstd = [(binzz_cstd_mean(1:end-1)+binzz_cstd_mean(2:end))/2]';
            
            tight_subplot(5,2,4,cond/2, guttera2, marginsa2)
           
            if model_pred
                h_p=fill([bincenterz_cstd bincenterz_cstd(end:-1:1)], [[squeeze(nanmean(p_pred_w_cstd(:,cond,Nind,:),1))]'-[squeeze(nanstd(p_pred_w_cstd(:,cond,Nind,:),1))]'./sqrt(squeeze(sum(~isnan(p_pred_w_cstd(:,cond,Nind,:)),1)))'...
                    [squeeze(nanmean(p_pred_w_cstd(:,cond,Nind,end:-1:1),1))]'+[squeeze(nanstd(p_pred_w_cstd(:,cond,Nind,end:-1:1),1))]'./sqrt(squeeze(sum(~isnan(p_pred_w_cstd(:,cond,Nind,:)),1)))'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
            end
            plot(bincenterz_cstd, squeeze(nanmean(p_data_w_cstd(:,cond,Nind,:),1)), 'o','Color', colors_sz(Nind,:), 'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Linestyle', 'none'); hold on;
            leg(Nind) = errorbar(bincenterz_cstd,squeeze(nanmean(p_data_w_cstd(:,cond,Nind,:),1)), squeeze(nanstd(p_data_w_cstd(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none');
            chance = plot(0:0.10:0.90, ones(1,10)*1/Nvec(Nind),'--', 'Color',colors_sz(Nind,:)); hold on;
            box off
            ylim([0 1])  
            xlim([0 1])
            set(gca, 'xtick',[0:0.2:1])
            set(gca, 'ytick', [0:0.1:1])
            set(gca, 'yticklabels', {})
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', [0.04 0.04])
            
     
            xlabel('circular variance of distractors','FontName', 'Helvetica', 'FontSize', fontsz)
            
            if cond == 2
                set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
                ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            
        end
        if cond == 4
            llpa = legend([leg(1) leg(2) leg(3) leg(4) chance],'N = 2', 'N = 3','N = 4', 'N = 6','chance','FontName','Helvetica','FontSize', fontsz)
            legend boxoff
            set(llpa, 'FontName','Helvetica','FontSize', fontsz*0.9, 'Position', [0.9 0.84 0.04 0.04])
        end
        
        
        
        %e) prop responses
        guttera3=[0.04 0.09]; %horiz vert
        marginsa3=[0.08 0.08 0.06 0.09]; %left right bottom top
        for Nind = 2:4
            
            tight_subplot(5,6, 5, 3*(cond/2-1) + Nind-1, guttera3, marginsa3)
            
            if model_predR
                h_p=fill([1:1:loc_vec(Nind) loc_vec(Nind):-1:1], [[[squeeze(mean(p_resp_ori_diff_rank_simim(:,cond,Nind,1:loc_vec(Nind)),1))]-[squeeze(std(p_resp_ori_diff_rank_simim(:,cond,Nind,1:loc_vec(Nind)),1))]/sqrt(Nsubj)]'...
                    [[squeeze(mean(p_resp_ori_diff_rank_simim(:,cond,Nind,loc_vec(Nind):-1:1),1))]+[squeeze(std(p_resp_ori_diff_rank_simim(:,cond,Nind,loc_vec(Nind):-1:1),1))]/sqrt(Nsubj)]'],colors_sz_shade(Nind,:),'EdgeColor', 'None'); hold on;
                %bar(1:1:(max(loc_vec)+0), squeeze(mean(p_resp_ori_diff_rank_simim(:,cond,Nind,1:(max(loc_vec)+0)),1)), 'EdgeColor', colors_sz(Nind,:), 'FaceColor', colors_sz_shade(Nind,:), 'BarWidth',0.8); hold on;
            end
            
            
            bar(1:1:(max(loc_vec)+0), squeeze(mean(p_resp_ori_diff_rank(:,cond,Nind,1:(max(loc_vec)+0)),1)), 'EdgeColor', colors_sz(Nind,:),'LineWidth', 1.1,  'FaceColor', 'None', 'BarWidth',0.46); hold on;
            errorbar(1:1:(max(loc_vec)+0),squeeze(mean(p_resp_ori_diff_rank(:,cond,Nind,1:(max(loc_vec)+0)),1)), squeeze(std(p_resp_ori_diff_rank(:,cond,Nind,1:(max(loc_vec)+0)),1))/sqrt(Nsubj),'Color', colors_sz(Nind,:));
            if Nind>1
                plot(1:1:loc_vec(Nind), (1-mean(p_c(:,cond, Nind),1))*chance_ori(Nind,1:loc_vec(Nind)), '--','Color',  'k', 'Linewidth',0.5); hold on;
            end
            
            
            box off
            ylim([0 p_resp_max])
            xlim([0.4 max(loc_vec)+0.6])
            set(gca, 'xtick', (0:1:(max(loc_vec)+1)))
            set(gca, 'ytick', [0:0.1:0.4])
            if Nind>2
                set(gca, 'yticklabels', []);
            end
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', [0.04 0.04])
            
            
            if cond ==2 & Nind == 4
                xlabel('Rank of T - R orientation distance','FontName', 'Helvetica', 'FontSize', fontsz)
            elseif Nind == 2
                ylabel('Proportion response','FontName', 'Helvetica', 'FontSize', fontsz)
                
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
            plot(bincentersl, squeeze(mean(p_data_w(:,cond,Nind,:),1)), 'o','Color', colors_sz(Nind,:),'MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot,'Linestyle', 'none'); hold on;
            leg(Nind) = errorbar(bincentersl,squeeze(mean(p_data_w(:,cond,Nind,:),1)), squeeze(std(p_data_w(:,cond,Nind,:),1))/sqrt(Nsubj), 'Color', colors_sz(Nind,:), 'Linestyle', 'none');
            chance = plot(0:10:90, ones(1,10)*1/Nvec(Nind),'--', 'Color',colors_sz(Nind,:)); hold on;
            box off
            ylim([0 1])
            xlim([0 90])
            set(gca, 'xtick', [0:15:90])
            set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', fontsz)
            set(gca, 'ytick', [0:0.1:1])
            set(gca, 'yticklabels', {})
            set(gca, 'tickdir', 'out')
            set(gca,'ticklength', [0.04 0.04])
            if cond == 2
                set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
                ylabel('Proportion correct','FontName', 'Helvetica', 'FontSize', fontsz)
            end
            xlabel('T - MSD orientation distance (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
            
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
%%
psname = ['LOC_full_figure_',num2str(full_figure),'_model_',num2str(mi),'_type_',num2str(type) ,'_exp_',num2str(exp_i), '.pdf']
%print_pdf(psname)
%% model comparison and parameter values
cd([curr_dir, '/model_fitting/'])
load(['nll_params_best_all.mat']) 
cd ..
cd([curr_dir, '/plotting/'])
%plot_params(exp_i, mi, type)
%plot_model_comparison(exp_i, type)


%% reaction times
%Localization
figure
set(gcf, 'Position', [100 100 500 220])
guttera2 = guttera;
marginsa2 = [0.14 0.14 0.2 0.2];
%b) T-most similar of distractors

for cond = [2 4]%[4 2]
    for Nind = 1:4
        bincentersl = [squeeze(bincenterz(cond, Nind,1,:))]';
        
        tight_subplot(1,2,1,cond/2, guttera2, marginsa2)
        
        leg(Nind) = plot(bincentersl, squeeze(mean(rt_data_w(:,cond,Nind,:),1)), 'o-','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Color', colors_sz(Nind,:)); hold on;
        errorbar(bincentersl,squeeze(mean(rt_data_w(:,cond,Nind,:),1)), squeeze(std(rt_data_w(:,cond,Nind,:),1))/sqrt(Nsubj),'Color', colors_sz(Nind,:),'Linestyle', 'none'); hold on;
        
        box off
        %ylim([0.2 0.7]) % for det
        %ylim([0.4 1]) % for loc, exp 1
        ylim([0.8 1.4]) % for loc, exp 2
        xlim([0 90])
        set(gca, 'xtick', [0:15:90])
        set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', fontsz)
        set(gca, 'ytick', [0:0.1:2])
        set(gca, 'yticklabels', {})
        if cond == 2
            set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1', '','1.2', '', '1.4', '', '1.6', '', '1.8' }, 'FontName', 'Helvetica', 'FontSize', fontsz)
            ylabel('Reaction time (sec)','FontName', 'Helvetica', 'FontSize', fontsz)
        end
        set(gca, 'tickdir', 'out')
        set(gca,'ticklength', [0.04 0.04])
        xlabel('T - MSD orientation distance (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
        if Nind == 1
            text(30, 1.46, text_labels{(cond)/2},  'FontName', 'Helvetica', 'FontSize', fontsz)
        end
        if cond == 2 & Nind == 4
            text(-40, 1.53, 'Localization',  'FontName', 'Helvetica', 'FontSize', fontsz*1.2)
        end
        
        if cond == 4
            llpa = legend([leg(1) leg(2) leg(3) leg(4)],'N = 2', 'N = 3','N = 4', 'N = 6','FontName','Helvetica','FontSize', fontsz)
            legend boxoff
            set(llpa, 'FontName','Helvetica','FontSize', fontsz*0.9, 'Position', [0.9 0.84 0.04 0.04])
        end
        
    end
    
end
%%
psname = ['RT_type_',num2str(type) ,'_exp_',num2str(exp_i), '.pdf']
%print_pdf(psname)


%%
%Detection
figure
set(gcf, 'Position', [100 100 500 220])
guttera2 = guttera;
marginsa2 = [0.14 0.14 0.2 0.2];
%b) T-most similar of distractors

for cond = [1 3]%[4 2]
    for Nind = 1:4
        bincentersl=[squeeze(bincenterz(cond, Nind,1,:))]';
        
        tight_subplot(1,2,1,0.5+cond/2, guttera2, marginsa2)
        
        leg(Nind)=plot(bincentersl, squeeze(mean(rt_data_w(:,cond,Nind,:),1)), '-o','MarkerEdgeColor', colors_sz(Nind,:), 'MarkerFaceColor', colors_sz(Nind,:), 'MarkerSize',sz_dot, 'Color', colors_sz(Nind,:)); hold on;
        errorbar(bincentersl,squeeze(mean(rt_data_w(:,cond,Nind,:),1)), squeeze(std(rt_data_w(:,cond,Nind,:),1))/sqrt(Nsubj),'Color', colors_sz(Nind,:),'Linestyle', 'none');
        chance = plot(0:10:90, 0.5*ones(1, 10),'--k'); hold on;
        box off
        ylim([0.15 0.7])
        xlim([0 90])
        set(gca, 'xtick', [0:15:90])
        set(gca, 'xticklabels', {'0', '', '30', '', '60', '', '90' },'FontName', 'Helvetica', 'FontSize', fontsz)
        set(gca, 'ytick', [0:0.1:1])
        set(gca, 'yticklabels', {})
        if cond == 1
            set(gca, 'yticklabels', {'0', '','0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
            ylabel('Reaction time (sec)','FontName', 'Helvetica', 'FontSize', fontsz)
        end
        set(gca, 'tickdir', 'out')
        set(gca,'ticklength', [0.04 0.04])
        xlabel('T - MSD orientation distance (deg)','FontName', 'Helvetica', 'FontSize', fontsz)
        if Nind ==1
            text(30, 0.76, text_labels{(1+cond)/2},  'FontName', 'Helvetica', 'FontSize', fontsz)
        end
        if cond == 1 & Nind == 4
            text(-40, 0.83, 'Detection',  'FontName', 'Helvetica', 'FontSize', fontsz*1.2)
        end
        
        if cond == 3
            llpa = legend([leg(1) leg(2) leg(3) leg(4)],'N = 2', 'N = 3','N = 4', 'N = 6','FontName','Helvetica','FontSize', fontsz)
            legend boxoff
            set(llpa, 'FontName','Helvetica','FontSize', fontsz*0.9, 'Position', [0.9 0.84 0.04 0.04])
        end
        
    end
    
end

%%
psname = ['RT_type_',num2str(type) ,'_exp_',num2str(exp_i), '.pdf']
%print_pdf(psname)




