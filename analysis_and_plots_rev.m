clear all; close all;

%exp_i = 1; % exp 1: stimulus spacing = 60 dva, exp 2: 30 dva
exp_i = 1;
load(['alldata_exp',num2str(exp_i),'.mat'])
mi = 1; % 1: wo decision noise, 2: decision noise
type = 1; %1: detection, 2: localization, 3: joint

model_pred = 1; % flag if we want to also calculate and plot model predictions
model_predR = 1;  % for model of specific localization responses, beyond prop corr
% added in June ---to not give error when asking for ddr for type 1, wont
% be present only in the model_pred matrices of type 2
if type == 1
    model_predR = 0;
end
Nsamp = 1200;

Nsubj = size(alldata,1);
Ncond = size(alldata,2);
Ntrials = length(alldata(1).data.N);
Nvec = unique(alldata(1).data.N)';

nbinz = 6;
nbinz2 = 6;
nbinz3 = 3;
nbinz3_2 = 2;
if exp_i == 1
    nbinz_jj = 5;%nbinz_jj = 3;
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
p_pr_c = nan(Nsubj, Ncond, length(Nvec));
preds = nan(Nsubj, n_models, Ncond, Ntrials);

p_resp_ori_diff_rank = nan(Nsubj, Ncond, length(Nvec), nbinz_jj+1);
p_resp_ori_diff_rank_simi = nan(Nsubj, Ncond, length(Nvec),nbinz_jj+1,Nsamp);   %sbjid, cond, Nind, jj, simi


binzz_meanE = nan(Nsubj, Ncond, length(Nvec), nbinz2+1);
binzzcvarE = nan(Nsubj, Ncond, length(Nvec), nbinz2+1);


b_fit_all = nan(Nsubj, Ncond, 3,3);  %we decided to do several regressions
% each with set size and either msd, mean and cvar


curr_dir = pwd;

cd([curr_dir, '/model_fitting/'])
load(['nll_params_best_all.mat'])
cd ..

n_types = 3;
n_vars = [6 5 6;  7 6 7]; % for models and types 1,2,3
params = nan(n_types,Nsubj,2,7); % 7 = max of npars
for modi = 1:2
    for ti = 1:n_types
        params(ti,1:Nsubj,1:2,1:n_vars(modi,ti)) = squeeze(params_all(exp_i,modi,ti,1:Nsubj,1:2,1:n_vars(modi,ti))); % 1 & 2 = perception & memory
    end
end

indi_rm = {[1],[1 2], [1 3], [1 4], [1 2 3 4], [1 3 4], [1 2 4], [1 2 3]};
coeff_len = [2 4 4 4 16 8 8 8] - 1; % we substract the intercept
var_indi_rm = {'N','min','mean', 'var','acc'};

regr_models_beta = nan(Nsubj, Ncond,length(indi_rm),max(coeff_len));
regr_models_beta_se = nan(Nsubj, Ncond,length(indi_rm),max(coeff_len));
regr_models_tStat = nan(Nsubj, Ncond,length(indi_rm),max(coeff_len));
regr_models_pValue = nan(Nsubj, Ncond,length(indi_rm), max(coeff_len));

regr_models_LL = nan(Nsubj, Ncond,length(indi_rm));
regr_models_AIC = nan(Nsubj, Ncond,length(indi_rm));
regr_models_BIC = nan(Nsubj, Ncond,length(indi_rm));
regr_models_Rsquared_adjusted = nan(Nsubj, Ncond,length(indi_rm));


load('mean_min_cvar_hist_type_1.mat') % will just vary by set size
binzz_mean_2_E = NaN(4,3);%binzz_mean_3 = [];
binzz_cvar_3_E = NaN(4,4);

for Nind = 2:4
    binzz_mean_2_E(Nind,1) = min(mean_ori_diff_hist(Nind,:));
    binzz_mean_2_E(Nind,2) = median(mean_ori_diff_hist(Nind,:));
    binzz_mean_2_E(Nind,3) = max(mean_ori_diff_hist(Nind,:));
    
    binzz_cvar_3_E(Nind,1) = min(cvar_ori_diff_hist(Nind,:));
    binzz_cvar_3_E(Nind,2) = quantile(cvar_ori_diff_hist(Nind,:),1/3);
    binzz_cvar_3_E(Nind,3) = quantile(cvar_ori_diff_hist(Nind,:),2/3);
    binzz_cvar_3_E(Nind,4) = quantile(cvar_ori_diff_hist(Nind,:),3/3);
    
    binzz_min_2_E(Nind,1) = min(min_ori_diff_hist(Nind,:));
    binzz_min_2_E(Nind,2) = median(min_ori_diff_hist(Nind,:));
    binzz_min_2_E(Nind,3) = max(min_ori_diff_hist(Nind,:));
end


% we doubled all our stimuli to the space (-pi,pi) for modeling ---
%0 has 3 categories, the rest have 2
nbinz_tval = 4;
bincenterz_target_val_deg = [0 45 90 135];%;./180*pi;
bincenterz_target_val = bincenterz_target_val_deg./180*pi;
tval_unit = (bincenterz_target_val_deg(2)-bincenterz_target_val_deg(1))/2;%10;%(bincenterz_target_val_deg(2)-bincenterz_target_val_deg(1))/2; % 22.5;
bin_tval_0 = [bincenterz_target_val_deg(1)-tval_unit bincenterz_target_val_deg(1)+tval_unit...
    180-tval_unit 180 -180 -180+tval_unit]./180*pi; % 0 and 180(=-180) 157.5 202.5
bin_tval_45 = [bincenterz_target_val_deg(2)-tval_unit bincenterz_target_val_deg(2)+tval_unit -180+bincenterz_target_val_deg(2)-tval_unit -180+bincenterz_target_val_deg(2)+tval_unit ]./180*pi; % 45 and 225 (=-135) 202.5 247.5 -157.5 -112.5
bin_tval_90 = [bincenterz_target_val_deg(3)-tval_unit bincenterz_target_val_deg(3)+tval_unit -180+bincenterz_target_val_deg(3)-tval_unit -180+bincenterz_target_val_deg(3)+tval_unit]./180*pi; % 90 and 270 (-90) 247.5 292.5
bin_tval_135 = [bincenterz_target_val_deg(4)-tval_unit bincenterz_target_val_deg(4)+tval_unit -180+bincenterz_target_val_deg(4)-tval_unit -180+bincenterz_target_val_deg(4)+tval_unit]./180*pi; %[112.5 157.5 -22.5 -67.5]./180*pi; % 135 is the same as -45
binz_target_val = {bin_tval_0, bin_tval_45, bin_tval_90, bin_tval_135};


load('mean_min_cvar_hist_type_2.mat')
binzz_mean_2_E_Loc = NaN(4,3);%binzz_mean_3 = [];
binzz_cvar_3_E_Loc = NaN(4,4);

for Nind = 2:4
    binzz_mean_2_E_Loc(Nind,1) = min(mean_ori_diff_hist_Loc(Nind,:));
    binzz_mean_2_E_Loc(Nind,2) = median(mean_ori_diff_hist_Loc(Nind,:));
    binzz_mean_2_E_Loc(Nind,3) = max(mean_ori_diff_hist_Loc(Nind,:));
    
    binzz_cvar_3_E_Loc(Nind,1) = min(cvar_ori_diff_hist_Loc(Nind,:));
    binzz_cvar_3_E_Loc(Nind,2) = quantile(cvar_ori_diff_hist_Loc(Nind,:),1/3);
    binzz_cvar_3_E_Loc(Nind,3) = quantile(cvar_ori_diff_hist_Loc(Nind,:),2/3);
    binzz_cvar_3_E_Loc(Nind,4) = quantile(cvar_ori_diff_hist_Loc(Nind,:),3/3);
end


min_ori_diff_on_lin_sep = NaN(Nsubj, Ncond, length(Nvec));
min_ori_diff_on_non_lin_sep = NaN(Nsubj, Ncond, length(Nvec));

for sbjid = 1:Nsubj
    
    sbjid
    
    % read in the output of modelling
    if model_pred
        
        model_fits_name = ['/model_pred/'];
        
        dirname = [curr_dir,model_fits_name];
        filez = dir([curr_dir,model_fits_name]);
        flst = {filez.name};
        
        cd(dirname)
        
        load(['model_pred_exp_', num2str(exp_i),'_model_',num2str(mi), '_type_', num2str(type),'_sbj_',num2str(sbjid),'.mat'])
        cd ..
        predz(mi,1,:) = pred(1,:);
        predz(mi,2,:) = pred(2,:);
        predz(mi,3,:) = pred(3,:);
        predz(mi,4,:) = pred(4,:);
        
        
        if type == 1 % det
            preds(sbjid,mi,1,:) = pred(1,:);
            preds(sbjid,mi,3,:) = pred(3,:);
        elseif type == 2 % loc
            preds(sbjid,mi,2,:) = pc(2,:); %pred(2,:);
            preds(sbjid,mi,4,:) = pc(4,:);%pred(4,:);
        elseif type == 3 %both
            preds(sbjid,mi,1,:) = pred(1,:);
            preds(sbjid,mi,3,:) = pred(3,:);
            
            preds(sbjid,mi,2,:) = pc(2,:); %pred(2,:);
            preds(sbjid,mi,4,:) = pc(4,:); %pred(4,:);
        end
        
    end
    
    
    for cond = 1:Ncond
        
        data = alldata(sbjid,cond).data;
        lin_sep1 = nan(Ntrials,1);
        
        % characterize prop target present and absent relative to
        % 3 summary statistics related to the distractors distribution
        % 1) ori distance from the target to the most similar distractor
        % 2) ori distance from target to the circular mean of the distractors
        % 3) circular standard deviation of the distractors
        acc_regr = nan(Ntrials,1);
        min_ori_diff = nan(Ntrials,1);
        mean_ori_diff = nan(Ntrials,1);
        cvar_ori_diff = nan(Ntrials,1);
        
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
            % 1) ori distance from the target to the most similar distractor
            % 2) ori distance from target to the circular mean of the distractors
            % 3) circular standard deviation of the distractors
            
            
            
            data_corr(sbjid, cond, 1:Ntrials) = (data.response == 1 & data.target_pres == 1) | ...
                (data.response == 0 & data.target_pres == 0);
            
            % split by set size
            for Nind = 1:length(Nvec)
                N = Nvec(Nind);
                idx_N = data.N == N;
                
                %data_corr(sbjid, cond, Nind,1:Ntrials/Ncond) = (data.response(idx_N) == 1 & data.target_pres(idx_N) == 1) | ...
                %   (data.response(idx_N) == 0 & data.target_pres(idx_N) == 0);
                
                p_tp(sbjid, cond, Nind) = sum(data.response(idx_pres & idx_N)==1)/sum(idx_pres & idx_N);
                p_ta(sbjid, cond, Nind) = sum(data.response(idx_abs & idx_N)==1)/sum(idx_abs & idx_N);
                
                n_tp(sbjid,cond,Nind) = sum(idx_pres & idx_N);
                n_ta(sbjid,cond,Nind) = sum(idx_abs & idx_N);
                
                rt(sbjid,cond, Nind) = median(data.reaction_time(idx_N));
                
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
                    
                    acc_regr(idx) = [data.response(idx) == 0];
                    
                    vall = [abs(bsxfun(@minus,data.stims(idx,:),data.target_val(idx)))]'; % in the interval (0,2*pi)
                    vall(vall>pi) = 2*pi-vall(vall>pi);
                    min_ori_diff(idx) = min(vall);% in the interval (0,pi)
                    
                    vall_mean = circ_mean(setdiff(data.stims(idx,find(~isnan(data.stims(idx,:)))), data.target_val(idx))');
                    vall_mean = abs(bsxfun(@minus,vall_mean, data.target_val(idx)));
                    vall_mean(vall_mean>pi) = 2*pi- vall_mean(vall_mean>pi);
                    mean_ori_diff(idx) = vall_mean;% in the interval (0,pi)
                    
                    vall_cvar_vec = (setdiff(data.stims(idx,find(~isnan(data.stims(idx,:)))), data.target_val(idx)));
                    vall_cvar_vec(vall_cvar_vec>pi) = 2*pi-vall_cvar_vec(vall_cvar_vec>pi);
                    vall_cvar = circ_var(vall_cvar_vec,[],[],[]);
                    cvar_ori_diff(idx) = vall_cvar; % in the interval (0,1)
                    
                    
                    vall_lin_sep = [(bsxfun(@minus,data.stims(idx,:),data.target_val(idx)))];%' % in the interval (0,2*pi)
                    vall_lin_sep(vall_lin_sep>pi) = vall_lin_sep(vall_lin_sep>pi) - 2* pi;
                    vall_lin_sep(vall_lin_sep<-pi) = vall_lin_sep(vall_lin_sep<-pi) + 2* pi;
                    if (sum(vall_lin_sep>=0)/ length(data.stims(idx,~isnan(data.stims(idx,:)))) == 1) |...
                            (sum(vall_lin_sep<=0)/ length(data.stims(idx,~isnan(data.stims(idx,:)))) == 1)
                        lin_sep1(idx) = 1;
                    end
                    
                end
                
                
                idx_target_pres = find(data.N ==N & data.target_pres==1);
                for j = 1:length(idx_target_pres)
                    idx = idx_target_pres(j);
                    acc_regr(idx) = [data.response(idx) == 1];
                    
                    vall = [abs(bsxfun(@minus,setdiff(data.stims(idx,:),data.target_val(idx)),data.target_val(idx)))]'; %btwn (0,2*pi)
                    vall(vall>pi) = 2*pi-vall(vall>pi);
                    min_ori_diff(idx) = min(vall);
                    
                    vall_mean = circ_mean(setdiff(data.stims(idx,find(~isnan(data.stims(idx,:)))), data.target_val(idx))');
                    vall_mean = abs(bsxfun(@minus,vall_mean, data.target_val(idx)));
                    vall_mean(vall_mean>pi) = 2*pi- vall_mean(vall_mean>pi);
                    mean_ori_diff(idx) = vall_mean;
                    
                    
                    vall_cvar_vec = (setdiff(data.stims(idx,find(~isnan(data.stims(idx,:)))), data.target_val(idx)));
                    vall_cvar_vec(vall_cvar_vec>pi)= 2*pi-vall_cvar_vec(vall_cvar_vec>pi);
                    vall_cvar = circ_var(vall_cvar_vec,[],[],[]);
                    cvar_ori_diff(idx) = vall_cvar;
                    
                    vall_lin_sep = [(bsxfun(@minus,data.stims(idx,:),data.target_val(idx)))];%' % in the interval (0,2*pi)
                    vall_lin_sep(vall_lin_sep>pi) = vall_lin_sep(vall_lin_sep>pi) - 2* pi;
                    vall_lin_sep(vall_lin_sep<-pi) = vall_lin_sep(vall_lin_sep<-pi) + 2* pi;
                    if (sum(vall_lin_sep>=0)/ length(data.stims(idx,~isnan(data.stims(idx,:)))) == 1) |...
                            (sum(vall_lin_sep<=0)/ length(data.stims(idx,~isnan(data.stims(idx,:)))) == 1)
                        lin_sep1(idx) = 1;
                    end
                end
                
                acc_lin_sep(sbjid, cond, Nind) = sum(acc_regr == 1 & lin_sep1 == 1 & data.N == N)/ sum(lin_sep1 == 1 & data.N == N);
                acc_non_lin_sep(sbjid, cond, Nind) = sum(acc_regr == 1 & lin_sep1~=1 & data.N == N)/ sum(lin_sep1 ~= 1 & data.N == N);
                n_trials_lin_sep(sbjid, cond, Nind) = sum(lin_sep1 == 1 & data.N == N); % out of 200 total
                
                min_ori_diff_on_lin_sep(sbjid, cond, Nind) =  mean(min_ori_diff(lin_sep1 == 1 & idx_N)); %median(min_ori_diff(lin_sep1 == 1 & idx_N));
                min_ori_diff_on_non_lin_sep(sbjid, cond, Nind) = mean(min_ori_diff(lin_sep1 ~= 1 & idx_N)); %median(min_ori_diff(lin_sep1 ~= 1 & idx_N));
                
                if model_pred  
                    
                    p_pr_tp_lin_sep(sbjid, cond, Nind) = mean(ptp(lin_sep1 == 1 & idx_pres & idx_N));
                    p_pr_ta_lin_sep(sbjid, cond, Nind) = mean(ptp(lin_sep1 == 1 & idx_abs & idx_N));
                    acc_lin_sep_pred(sbjid, cond, Nind) = 1-((1-p_pr_tp_lin_sep(sbjid, cond, Nind))*sum(lin_sep1 == 1 & idx_pres & idx_N) + p_pr_ta_lin_sep(sbjid, cond, Nind)* sum(lin_sep1 == 1 & idx_abs & idx_N))/ sum(lin_sep1 == 1 & idx_N);
                
                    p_pr_tp_non_lin_sep(sbjid, cond, Nind) = mean(ptp(lin_sep1 ~= 1 & idx_pres & idx_N));
                    p_pr_ta_non_lin_sep(sbjid, cond, Nind) = mean(ptp(lin_sep1 ~= 1 & idx_abs & idx_N));
                    acc_non_lin_sep_pred(sbjid, cond, Nind) = 1-((1-p_pr_tp_non_lin_sep(sbjid, cond, Nind))*sum(lin_sep1 ~= 1 & idx_pres & idx_N) + p_pr_ta_non_lin_sep(sbjid, cond, Nind)* sum(lin_sep1 ~= 1 & idx_abs & idx_N))/ sum(lin_sep1 ~= 1 & idx_N);
                
                end
                
                % split by target_val
                for j = 1:nbinz_tval
                    bin_target_val = binz_target_val{j};
                    if j == 1
                        idx_bin_tval  = data.target_val <bin_target_val(2) & data.target_val>bin_target_val(1) |...
                            data.target_val <bin_target_val(4) & data.target_val>bin_target_val(3) |...
                            data.target_val <bin_target_val(6) & data.target_val>bin_target_val(5) ;
                    elseif j> 1
                        idx_bin_tval  = data.target_val <bin_target_val(2) & data.target_val>bin_target_val(1) |...
                            data.target_val <bin_target_val(4) & data.target_val>bin_target_val(3) ;
                    end
                    
                    idx_t = idx_N & idx_bin_tval; % & idx_pres;
                    
                    p_data_tval_pc(sbjid, cond, Nind, j) = sum(data_corr(sbjid, cond, idx_t)==1)/sum(idx_t);
                    %p_data_tp(sbjid, cond, Nind, j) = sum(data.response(idx_p)==1)/sum(idx_p);
                    %n_data_tp(sbjid, cond, Nind, j) = sum(idx_p);
                    
                    %idx_a = idx_N & idx_bin_a & idx_abs;
                    %p_data_ta(sbjid, cond, Nind, j) = sum(data.response(idx_a)==1)/sum(idx_a);
                    %n_data_ta(sbjid, cond, Nind, j) = sum(idx_a);
                    
                    %error_rate_all(sbjid, cond, Nind,j)=((1-p_data_tp(sbjid, cond, Nind,j))*sum(idx_p) + p_data_ta(sbjid, cond, Nind,j)* sum(idx_a))/ (sum(idx_p)+sum(idx_a));
                    
                    %rt_data_w_msd(sbjid, cond, Nind, j) = median(data.reaction_time(idx_N & idx_bin_p));
                    
                    rt_data_tval_pc(sbjid, cond, Nind, j) = median(data.reaction_time(idx_t));
                    
                    if model_pred
                        
                        p_pred_tval_tp(sbjid, cond, Nind, j) = mean(ptp(idx_N & idx_bin_tval & idx_pres));
                        p_pred_tval_ta(sbjid, cond, Nind, j) = mean(ptp(idx_N & idx_bin_tval & idx_abs));
                        error_rate_pred_tval_all(sbjid, cond, Nind,j) = ((1-p_pred_tval_tp(sbjid, cond, Nind,j))*sum(idx_N & idx_bin_tval & idx_pres) + p_pred_tval_ta(sbjid, cond, Nind,j)* sum(idx_N & idx_bin_tval & idx_abs))/ (sum(idx_N & idx_bin_tval & idx_pres)+sum(idx_N & idx_bin_tval & idx_abs));
                        p_pred_tval_pc(sbjid, cond, Nind, j) = 1 - error_rate_pred_tval_all(sbjid, cond, Nind,j);
                    end   
                end
                
                
                % split by target_val 
                if Nind >= 2
                    binzz_mean_2 = squeeze(binzz_mean_2_E(Nind,:));
                    binzz_min_2= squeeze(binzz_min_2_E(Nind,:));
                    for jj  = 1: nbinz3_2
                        idx_bin_meanE3  = min_ori_diff<binzz_min_2(jj+1) & min_ori_diff>binzz_min_2(jj);
                        
                        
                        for j = 1:nbinz_tval
                            bin_target_val = binz_target_val{j};
                            if j == 1
                                idx_bin_tval  = data.target_val <bin_target_val(2) & data.target_val>bin_target_val(1) |...
                                    data.target_val <bin_target_val(4) & data.target_val>bin_target_val(3) |...
                                    data.target_val <bin_target_val(6) & data.target_val>bin_target_val(5) ;
                            elseif j> 1
                                idx_bin_tval  = data.target_val <bin_target_val(2) & data.target_val>bin_target_val(1) |...
                                    data.target_val <bin_target_val(4) & data.target_val>bin_target_val(3) ;
                            end
                            
                            idx_t = idx_N & idx_bin_tval & idx_bin_meanE3; % & idx_pres;
                            
                            p_data_tval_EE_pc(sbjid, cond, Nind, jj,j) = sum(data_corr(sbjid, cond, idx_t)==1)/sum(idx_t);
                            rt_data_tval_EE_pc(sbjid, cond, Nind, jj, j) = median(data.reaction_time(idx_t));
                            
                            
                        end
                    end
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
                    
                    rt_data_w_msd(sbjid, cond, Nind, j) = median(data.reaction_time(idx_N & idx_bin_p));
                    
                    if model_pred
                        p_pred_tp(sbjid, cond, Nind, j) = mean(ptp(idx_p));
                        p_pred_ta(sbjid, cond, Nind, j) = mean(ptp(idx_a));
                        error_rate_pred_all(sbjid, cond, Nind,j) = ((1-p_pred_tp(sbjid, cond, Nind,j))*sum(idx_p) + p_pred_ta(sbjid, cond, Nind,j)* sum(idx_a))/ (sum(idx_p)+sum(idx_a));
                    end

                end
                
                binzz_mean = [];
                binzz_cvar = [];
                for j = 1: nbinz2
                    binzz_mean(j) = quantile(mean_ori_diff(data.N==N),j/(nbinz2));
                    binzz_cvar(j) = quantile(cvar_ori_diff(data.N==N),j/(nbinz2));
                end
                binzz_mean = [min(mean_ori_diff(data.N==N)) binzz_mean];
                binzz_cvar = [min(cvar_ori_diff(data.N==N)) binzz_cvar];
                
                binzz_meanE(sbjid,cond,Nind,:) = binzz_mean;
                binzz_cvarE(sbjid,cond,Nind,:) = binzz_cvar;
                
                %split by mean of distractors and cvar of distractors
                for j = 1:nbinz2
                    
                    idx_bin_mean  = mean_ori_diff<binzz_mean(j+1) & mean_ori_diff>binzz_mean(j);
                    idx_bin_cvar  = cvar_ori_diff<binzz_cvar(j+1) & cvar_ori_diff>binzz_cvar(j);
                    
                    idxb_mean_p = idx_pres & idx_N & idx_bin_mean;
                    idxb_mean_a = idx_abs & idx_N & idx_bin_mean;
                    
                    idxb_cvar_p = idx_pres & idx_N & idx_bin_cvar;
                    idxb_cvar_a = idx_abs & idx_N & idx_bin_cvar;
                    
                    p_data_tp_mean(sbjid, cond, Nind, j) = sum(data.response(idxb_mean_p)==1)/sum(idxb_mean_p);
                    p_data_ta_mean(sbjid, cond, Nind, j) = sum(data.response(idxb_mean_a)==1)/sum(idxb_mean_a);
                    
                    p_data_tp_cvar(sbjid, cond, Nind, j) = sum(data.response(idxb_cvar_p)==1)/sum(idxb_cvar_p);
                    p_data_ta_cvar(sbjid, cond, Nind, j) = sum(data.response(idxb_cvar_a)==1)/sum(idxb_cvar_a);
                    
                    n_data_w_mean(sbjid, cond, Nind, j) = sum(idx_N & idx_bin_mean);
                    n_data_w_cvar(sbjid, cond, Nind, j) = sum(idx_N & idx_bin_cvar);
                    
                    rt_data_w_mean(sbjid,cond, Nind, j) =  median(data.reaction_time(idx_N & idx_bin_mean));
                    rt_data_w_cvar(sbjid,cond, Nind, j) =  median(data.reaction_time(idx_N & idx_bin_cvar));
                    
                    error_rate_mean_all(sbjid, cond, Nind,j) = ((1-p_data_tp_mean(sbjid, cond, Nind,j))*sum(idxb_mean_p) + p_data_ta_mean(sbjid, cond, Nind,j)* sum(idxb_mean_a))/ (sum(idxb_mean_p)+sum(idxb_mean_a));
                    error_rate_cvar_all(sbjid, cond, Nind,j) = ((1-p_data_tp_cvar(sbjid, cond, Nind,j))*sum(idxb_cvar_p) + p_data_ta_cvar(sbjid, cond, Nind,j)* sum(idxb_cvar_a))/ (sum(idxb_cvar_p)+sum(idxb_cvar_a));
                    
                    
                    if model_pred
                        p_pred_tp_mean(sbjid, cond, Nind, j) = mean(ptp(idxb_mean_p));
                        p_pred_ta_mean(sbjid, cond, Nind, j) = mean(ptp(idxb_mean_a));
                        
                        error_rate_mean_pred_all(sbjid, cond, Nind,j)=((1-p_pred_tp_mean(sbjid, cond, Nind,j))*sum(idxb_mean_p) + p_pred_ta_mean(sbjid, cond, Nind,j)* sum(idxb_mean_a))/ (sum(idxb_mean_p)+sum(idxb_mean_a));
                        
                        p_pred_tp_cvar(sbjid, cond, Nind, j) = mean(ptp(idxb_cvar_p));
                        p_pred_ta_cvar(sbjid, cond, Nind, j) = mean(ptp(idxb_cvar_a));
                        
                        error_rate_cvar_pred_all(sbjid, cond, Nind,j) = ((1-p_pred_tp_cvar(sbjid, cond, Nind,j))*sum(idxb_cvar_p) + p_pred_ta_cvar(sbjid, cond, Nind,j)* sum(idxb_cvar_a))/ (sum(idxb_cvar_p)+sum(idxb_cvar_a));
                    end
                    
                end
                
                if Nind >= 2
                    binzz_mean_2 = squeeze(binzz_mean_2_E(Nind,:));
                    binzz_cvar_3 = squeeze(binzz_cvar_3_E(Nind,:));
                    
                    for jj  = 1: nbinz3_2
                        idx_bin_meanE3  = mean_ori_diff<binzz_mean_2(jj+1) & mean_ori_diff>binzz_mean_2(jj);
                        
                        for ll = 1: nbinz3
                            idx_bin_cvarE3  = cvar_ori_diff<binzz_cvar_3(ll+1) & cvar_ori_diff>binzz_cvar_3(ll);
                            
                            idxb_pres_mean_cvar = idx_pres & idx_N & idx_bin_meanE3 & idx_bin_cvarE3;
                            idxb_abs_mean_cvar = idx_abs & idx_N & idx_bin_meanE3 & idx_bin_cvarE3;
                            
                            p_data_tp_mean_cvar(sbjid, cond, Nind, jj,ll) = sum(data.response(idxb_pres_mean_cvar)==1)/sum(idxb_pres_mean_cvar);
                            p_data_ta_mean_cvar(sbjid, cond, Nind, jj,ll) = sum(data.response(idxb_abs_mean_cvar)==1)/sum(idxb_abs_mean_cvar);
                            
                            error_rate_mean_cvar_all(sbjid, cond, Nind,jj,ll) = ((1-p_data_tp_mean_cvar(sbjid, cond, Nind,jj,ll))*sum(idxb_pres_mean_cvar) + p_data_ta_mean_cvar(sbjid, cond, Nind,jj,ll)* sum(idxb_abs_mean_cvar))/ (sum(idxb_pres_mean_cvar)+sum(idxb_abs_mean_cvar));
                            
                            % now overall, not broken down also by T-D mean
                            %idx_bin_cvarE3  = cvar_ori_diff<binzz_cvar_3(ll+1) & cvar_ori_diff>binzz_cvar_3(ll);
                            
                            idxb_pres_mean_cvar_ONLY = idx_pres & idx_N &  idx_bin_cvarE3;
                            idxb_abs_mean_cvar_ONLY = idx_abs & idx_N & idx_bin_cvarE3;
                            
                            p_data_tp_mean_cvar_ONLY(sbjid, cond, Nind, jj,ll) = sum(data.response(idxb_pres_mean_cvar_ONLY)==1)/sum(idxb_pres_mean_cvar_ONLY);
                            p_data_ta_mean_cvar_ONLY(sbjid, cond, Nind, jj,ll) = sum(data.response(idxb_abs_mean_cvar_ONLY)==1)/sum(idxb_abs_mean_cvar_ONLY);
                            
                            error_rate_mean_cvar_all_ONLY(sbjid, cond, Nind,jj,ll) = ((1-p_data_tp_mean_cvar_ONLY(sbjid, cond, Nind,jj,ll))*sum(idxb_pres_mean_cvar_ONLY) + p_data_ta_mean_cvar_ONLY(sbjid, cond, Nind,jj,ll)* sum(idxb_abs_mean_cvar_ONLY))/ (sum(idxb_pres_mean_cvar_ONLY)+sum(idxb_abs_mean_cvar_ONLY));
                     
                        end
                    end
                    
                end
                
            end
            
            for szi = 1:3
                %Nvec_szi = Nvec(szi+1);
                [r_corr4(sbjid,cond,1,szi), p_corr4(sbjid,cond,1,szi)] = circ_corrcc(min_ori_diff(data.N==Nvec(szi+1)),mean_ori_diff(data.N==Nvec(szi+1)));%,'type', 'Spearman');
                [r_corr4(sbjid,cond,2,szi), p_corr4(sbjid,cond,2,szi)]= circ_corrcc(min_ori_diff(data.N==Nvec(szi+1)),cvar_ori_diff(data.N==Nvec(szi+1)));%,'type', 'Spearman');
                [r_corr4(sbjid,cond,3,szi), p_corr4(sbjid,cond,3,szi)] = circ_corrcc(mean_ori_diff(data.N==Nvec(szi+1)),cvar_ori_diff(data.N==Nvec(szi+1))); %,'type', 'Spearman');
            end
            
            
            x_regrr = [(data.N(data.N>2) - mean(data.N(data.N>2)))/std(data.N(data.N>2)) (min_ori_diff(data.N>2) - mean(min_ori_diff(data.N>2)))/std(min_ori_diff(data.N>2))...
                (mean_ori_diff(data.N>2) - mean(mean_ori_diff(data.N>2)))/std(mean_ori_diff(data.N>2))...
                (cvar_ori_diff(data.N>2) - mean(cvar_ori_diff(data.N>2)))/std(cvar_ori_diff(data.N>2))];
            
            for rmi = 1:length(indi_rm) % regression model index
                dsa = [x_regrr(:, indi_rm{rmi}) acc_regr(data.N>2)];
                varNames = cell(1+length(indi_rm{rmi}),1);
                
                modelspec = [ var_indi_rm{end},'~'];
                indi_rm_rmi = indi_rm{rmi};
                varNames{1} = var_indi_rm{indi_rm_rmi(1)};
                modelspec = [modelspec var_indi_rm{indi_rm_rmi(1)}];
                if length(indi_rm{rmi})>1
                    for rmij = 2: length(indi_rm{rmi})
                        varNames{rmij} =  var_indi_rm{indi_rm{rmi}(rmij)};
                        modelspec = [modelspec '*', var_indi_rm{indi_rm{rmi}(rmij)}];  %'acc~ N*min';
                    end
                else
                    rmij = 1;
                end
                
                varNames{rmij+1} = var_indi_rm{end};
                dsaT = array2table(dsa, 'VariableNames',varNames);
                mdlALL{sbjid, cond, rmi} = fitglm(dsaT,modelspec,'Distribution','binomial', 'link', 'logit');
                
                
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
            
            
            all_min_ori_diff = [all_min_ori_diff min_ori_diff];
        else %cond even, localization data
            
            idx_corr  = (data.target_loc == data.response) ;
            idx_incorr  = (data.target_loc ~= data.response) ;
            
            
            loc_diff = nan(Ntrials,1);
            resp_ori_diff_rank = nan(Ntrials,1);
            resp_loc_diff_rank = nan(Ntrials,1);
            
            for Nind = 1:length(Nvec)
                
                N = Nvec(Nind);
                idx_N = data.N == N;
                
                if model_pred
                    pc = squeeze(preds(sbjid,mi,cond,:));
                    
                end
                
                p_c(sbjid, cond, Nind) = sum(idx_corr & idx_N)/sum(idx_N);
                rt(sbjid, cond, Nind) = median(data.reaction_time(idx_N));
                %rt_c_corr(sbjid, cond, Nind) = median(data.reaction_time(idx_corr & idx_N));
                %rt_c_incorr(sbjid, cond, Nind) = median(data.reaction_time(idx_incorr & idx_N));
                
                n_c(sbjid,cond,Nind) = sum(idx_N);
                
                
                error_rate(sbjid,cond,Nind) = 1-p_c(sbjid, cond, Nind);
                if model_pred
                    p_pr_c(sbjid, cond, Nind) = mean(pc(idx_corr & idx_N));
                end
                
                idx_N_val = find(data.N==N);
                
                for j = 1:length(idx_N_val)
                    idx = idx_N_val(j);
                    
                    acc_regr(idx) = [data.response(idx) == data.target_loc(idx)];
                    
                    vall = [abs(bsxfun(@minus,data.stims(idx,:),data.target_val(idx)))]'; %btwn (0,2*pi)
                    vall(vall>pi) = 2*pi-vall(vall>pi);
                    min_ori_diff(idx) = min(vall(vall>0)); %min(vall); %btwn (0,pi)
                    
                    vall2 = abs(data.stims(idx,data.response(idx)) - data.target_val(idx));
                    if vall2 > pi
                        vall2 = 2*pi-vall2;
                    end
                    vall_sorted = sort(vall);
                    
                    vall_mean = circ_mean(setdiff(data.stims(idx,find(~isnan(data.stims(idx,:)))), data.target_val(idx))');
                    vall_mean = abs(bsxfun(@minus,vall_mean, data.target_val(idx)));
                    vall_mean(vall_mean>pi) = 2*pi- vall_mean(vall_mean>pi);
                    mean_ori_diff(idx) = vall_mean;
                    
                    vall_cvar_vec = (setdiff(data.stims(idx,find(~isnan(data.stims(idx,:)))), data.target_val(idx)));
                    vall_cvar_vec(vall_cvar_vec>pi) = 2*pi-vall_cvar_vec(vall_cvar_vec>pi);
                    vall_cvar = circ_var(vall_cvar_vec,[],[],[]);
                    cvar_ori_diff(idx) = vall_cvar;

                    vall_lin_sep =   bsxfun(@minus,data.stims(idx,:),data.target_val(idx));
                    vall_lin_sep(vall_lin_sep>pi) = vall_lin_sep(vall_lin_sep>pi) - 2* pi;
                    vall_lin_sep(vall_lin_sep<-pi) = vall_lin_sep(vall_lin_sep<-pi) + 2* pi;
                    if (sum(vall_lin_sep>=0)/ length(data.stims(idx,~isnan(data.stims(idx,:)))) == 1) |...
                            (sum(vall_lin_sep<=0)/ length(data.stims(idx,~isnan(data.stims(idx,:)))) == 1)
                        lin_sep1(idx) = 1;
                    end
                    
                    if model_predR
                        for simi = 1: Nsamp
                            data_stims_idx = data.stims(idx,~isnan(data.stims(idx,:)));
                            vall3 = abs(data_stims_idx(ddr(cond,idx,simi)) - data.target_val(idx));
                            if vall3 > pi
                                vall3 = 2*pi-vall3;
                            end
                            resp_ori_diff_rank_simi(idx,simi) = find(vall3==vall_sorted)-1;
                            resp_loc_diff_rank_simi(idx,simi) = abs(find(vall==0) - find(vall==vall3));
                            
                        end
                    end
                    
                    resp_ori_diff_rank(idx) = find(vall2==vall_sorted)-1;
                    
                end
                
                acc_lin_sep(sbjid, cond, Nind) = sum(acc_regr == 1 & lin_sep1 == 1 & data.N == N)/ sum(lin_sep1 == 1 & data.N == N);
                acc_non_lin_sep(sbjid, cond, Nind) = sum(acc_regr == 1 & lin_sep1 ~= 1 & data.N == N)/ sum(lin_sep1 ~= 1 & data.N == N);
                n_trials_lin_sep(sbjid, cond, Nind) = sum(lin_sep1 == 1 & data.N == N); % out of 200 total
                
                min_ori_diff_on_lin_sep(sbjid, cond, Nind) = mean(min_ori_diff(lin_sep1 == 1 & idx_N)); %median(min_ori_diff(lin_sep1 == 1 & idx_N));
                min_ori_diff_on_non_lin_sep(sbjid, cond, Nind) = mean(min_ori_diff(lin_sep1 ~= 1 & idx_N)); %median(min_ori_diff(lin_sep1 ~= 1 & idx_N));
                
                if model_pred
                    acc_lin_sep_pred(sbjid, cond, Nind) = mean(pc(lin_sep1 == 1 & data.N == N));
                    acc_non_lin_sep_pred(sbjid, cond, Nind) = mean(pc(lin_sep1 ~= 1 & data.N == N));
                end
                
                %split by order, prepare to analyze response choices
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
                    
                    if jj <= Nvec(Nind)-1
                        idx_loco = [(resp_ori_diff_rank==jj)];
                        idx_lo = idx_loco & idx_N;
                        p_resp_ori_diff_rank(sbjid, cond, Nind, jj) = sum(idx_lo)/length(idx_N_val);
                    end
                    if jj <= loc_vec(Nind)
                        idx_locl = [(resp_loc_diff_rank==jj)];
                        idx_ll = idx_locl & idx_N;
                        p_resp_loc_diff_rank(sbjid, cond, Nind,jj) = sum(idx_ll)/length(idx_N_val);
                    end
                    
                    if model_predR
                        clear idx_loco; clear idx_lo;  clear idx_ll;
                        
                        for simi = 1:Nsamp
                            if jj <= Nvec(Nind)-1
                                idx_loco = [(resp_ori_diff_rank_simi(:,simi)==jj)];
                                idx_lo = idx_loco & idx_N;
                                p_resp_ori_diff_rank_simi(sbjid, cond, Nind, jj, simi) = sum(idx_lo)/length(idx_N_val);    
                            end   
                        end
                    end
                    
                end
                
                for j = 1:nbinz
                    idx_bin  = min_ori_diff<binz(cond,Nind,1,j+1) & min_ori_diff>binz(cond,Nind,1,j);
                    idxb = idx_corr & idx_N & idx_bin;
                    idx_bin_rt = idx_bin;
                    
                    p_data_w(sbjid, cond, Nind, j) = sum(idxb)/sum(idx_N & idx_bin);
                    n_data_w(sbjid, cond, Nind, j) = sum(idx_N & idx_bin);
                    
                    rt_data_w_msd(sbjid, cond, Nind, j) = median(data.reaction_time(idx_N & idx_bin));
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
                binzz_cvar = [];
                for j = 1: nbinz2
                    binzz_mean(j) = quantile(mean_ori_diff(data.N==N),j/(nbinz2));
                    binzz_cvar(j) = quantile(cvar_ori_diff(data.N==N),j/(nbinz2));
                end
                binzz_mean = [min(mean_ori_diff(data.N==N)) binzz_mean];
                binzz_cvar = [min(cvar_ori_diff(data.N==N)) binzz_cvar];
                
                binzz_meanE(sbjid,cond,Nind,:) = binzz_mean;
                binzz_cvarE(sbjid,cond,Nind,:) = binzz_cvar;
                
                for j = 1:nbinz2
                    idx_bin_mean  = mean_ori_diff<binzz_mean(j+1) & mean_ori_diff>binzz_mean(j);
                    idx_bin_cvar  = cvar_ori_diff<binzz_cvar(j+1) & cvar_ori_diff>binzz_cvar(j);
                    
                    idxb_mean = idx_corr & idx_N & idx_bin_mean;
                    idxb_cvar = idx_corr & idx_N & idx_bin_cvar;
                    
                    p_data_w_mean(sbjid, cond, Nind, j) = sum(idxb_mean)/sum(idx_N & idx_bin_mean);
                    p_data_w_cvar(sbjid, cond, Nind, j) = sum(idxb_cvar)/sum(idx_N & idx_bin_cvar);
                    
                    n_data_w_mean(sbjid, cond, Nind, j) = sum(idx_N & idx_bin_mean);
                    n_data_w_cvar(sbjid, cond, Nind, j) = sum(idx_N & idx_bin_cvar);
                    
                    rt_data_w_mean(sbjid,cond, Nind, j) =  median(data.reaction_time(idx_N & idx_bin_mean));
                    rt_data_w_cvar(sbjid,cond, Nind, j) =  median(data.reaction_time(idx_N & idx_bin_cvar));
                    
                    if model_pred
                        
                        p_pred_w_mean(sbjid, cond, Nind, j) = mean(pc(idxb_mean));
                        p_pred_w_cvar(sbjid, cond, Nind, j) = mean(pc(idxb_cvar));
                    end
                end
                
            end
            
            
            for szi = 1:3
                [r_corr4(sbjid,cond,1,szi), p_corr4(sbjid,cond,1,szi)] = circ_corrcc(min_ori_diff(data.N==Nvec(szi+1)),mean_ori_diff(data.N==Nvec(szi+1))); %,'type', 'Spearman');
                [r_corr4(sbjid,cond,2,szi), p_corr4(sbjid,cond,2,szi)]= circ_corrcc(min_ori_diff(data.N==Nvec(szi+1)),cvar_ori_diff(data.N==Nvec(szi+1))); %,'type', 'Spearman');
                [r_corr4(sbjid,cond,3,szi), p_corr4(sbjid,cond,3,szi)] = circ_corrcc(mean_ori_diff(data.N==Nvec(szi+1)),cvar_ori_diff(data.N==Nvec(szi+1))); %,'type', 'Spearman');
            end
            
            x_regrr = [(data.N(data.N>2) - mean(data.N(data.N>2)))/std(data.N(data.N>2)) (min_ori_diff(data.N>2) - mean(min_ori_diff(data.N>2)))/std(min_ori_diff(data.N>2))...
                (mean_ori_diff(data.N>2) - mean(mean_ori_diff(data.N>2)))/std(mean_ori_diff(data.N>2))...
                (cvar_ori_diff(data.N>2) - mean(cvar_ori_diff(data.N>2)))/std(cvar_ori_diff(data.N>2))];
            
            
            for rmi = 1:length(indi_rm) % regression model index
                dsa = [x_regrr(:, indi_rm{rmi}) acc_regr(data.N>2)];
                varNames = cell(1+length(indi_rm{rmi}),1);
                
                modelspec = [ var_indi_rm{end},'~'];
                indi_rm_rmi = indi_rm{rmi};
                varNames{1} = var_indi_rm{indi_rm_rmi(1)};
                modelspec = [modelspec var_indi_rm{indi_rm_rmi(1)}];
                if length(indi_rm{rmi})>1
                    for rmij = 2: length(indi_rm{rmi})
                        varNames{rmij} =  var_indi_rm{indi_rm{rmi}(rmij)};
                        modelspec = [modelspec '*', var_indi_rm{indi_rm{rmi}(rmij)}];  %'acc~ N*min';
                    end
                else
                    rmij = 1;
                    
                end
                
                varNames{rmij+1} = var_indi_rm{end};
                dsaT = array2table(dsa, 'VariableNames',varNames);
                mdlALL{sbjid, cond, rmi} = fitglm(dsaT,modelspec,'Distribution','binomial', 'link', 'logit');
                
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
        
        all_min_ori_diff = [all_min_ori_diff min_ori_diff];
        all_mean_ori_diff = [all_mean_ori_diff mean_ori_diff];
        
        
    end
end


if model_pred
    cd ..
    close all;
end



%%
binzz_rank = [0 2 3 5];
for si = 1: Nsubj
    for pmi = 1:2
        for ci = 2:4
            x = squeeze(p_resp_ori_diff_rank(si,2*pmi,ci,:));
            x = x(~isnan(x));
            nbinss = binzz_rank(ci); % number of bin
            edges = linspace(0,nansum(squeeze(p_resp_ori_diff_rank(si,2*pmi,ci,:))'),nbinss+1); % edges of the bins
            E = nansum(squeeze(p_resp_ori_diff_rank(si,2*pmi,ci,:)))'/nbinss*ones(nbinss,1); % expected value (equal for uniform dist)
            
            [h,p_chi(si,pmi,ci),stats] = chi2gof(x','Expected',E','Edges',edges)
        end
    end
end

%% define bins


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

%% accuracy with metrics
accuracy=bsxfun(@minus,1,error_rate);
accuracy_msd_all = bsxfun(@minus,1,error_rate_all); % really accuracy_msd_all
accuracy_mean_all = bsxfun(@minus,1,error_rate_mean_all);
accuracy_cvar_all = bsxfun(@minus,1,error_rate_cvar_all);

if model_pred
    accuracy_pred = bsxfun(@minus,1,error_rate_pred);
    accuracy_pred_msd_all = bsxfun(@minus,1,error_rate_pred_all);
    accuracy_pred_mean_all = bsxfun(@minus,1, error_rate_mean_pred_all);
    accuracy_pred_cvar_all = bsxfun(@minus,1, error_rate_cvar_pred_all);
end

%%
cd(curr_dir)

savefilename_all_vars = ['all_vars_exp_',num2str(exp_i),'_type_',num2str(type) ,'_mi_',num2str(mi), '.mat']
%save(savefilename_pars_nll,'params','nll_sbj','-mat')
save(savefilename_all_vars,'exp_i', 'mi','type','Nvec','Nsubj','nbinz_jj','loc_vec','p_tp', 'p_ta','p_pr_tp', 'p_pr_ta',...
    'bincenterz','binzz_meanE', 'binzz_cvarE', 'binzz_cvar_3_E','binzz_mean_2_E','bincenterz_target_val',...
    'accuracy', 'acc_lin_sep', 'acc_non_lin_sep', 'n_trials_lin_sep', 'acc_lin_sep_pred', 'acc_non_lin_sep_pred', 'p_data_tval_pc', 'p_data_tval_EE_pc', 'accuracy_msd_all', 'accuracy_pred_msd_all',...
    'accuracy_mean_all','accuracy_pred_mean_all','accuracy_cvar_all', 'accuracy_pred_cvar_all',...
    'error_rate_mean_cvar_all','error_rate_mean_cvar_all_ONLY',...
    'p_pr_tp', 'p_pr_ta','p_data_tp_cvar', 'p_data_ta_cvar', 'p_data_tp_mean_cvar', 'p_data_ta_mean_cvar',...
    'p_data_tp_mean_cvar_ONLY', 'p_data_ta_mean_cvar_ONLY', 'p_pred_tp_cvar', 'p_pred_ta_cvar',...
    'p_pred_tval_pc', ...
    'p_c','p_pr_c','p_data_w','p_data_w_mean','p_data_w_cvar','p_pred_w','p_pred_w_mean','p_pred_w_cvar',...
    'p_data_tp','p_data_ta','p_pred_tp','p_pred_ta','p_data_tp_mean','p_data_ta_mean', 'p_data_tp_cvar','p_data_ta_cvar',...
    'p_pred_tp_mean', 'p_pred_ta_mean','p_pred_tp_cvar', 'p_pred_ta_cvar',....
    'p_resp_ori_diff_rank','p_resp_ori_diff_rank_simi',...
    'rt','rt_data_w_msd','rt_data_w_mean','rt_data_w_cvar','-mat')
%


%%
%{
exp_i = 1;
type = 2;
mi = 1;
load(['all_vars_exp_',num2str(exp_i),'_type_',num2str(type) ,'_mi_',num2str(mi), '.mat'])
%}

%% Detection % if type ==1 or 3
model_pred = 1;
if type == 1 || 3
    full_figure_flag = 1;%0; % alternative is reduced figure, just with T-MSD
    save_flag = 1;
    plot_det(exp_i,type,mi, model_pred, full_figure_flag, save_flag) % will load all_vars inside fcn
end

%% Localization  %% if type = 2 or 3
model_pred = 1;
model_predR = 1;
save_flag = 1;
if type == 2 || 3
    full_figure_flag =  1; % alternative is reduced figure, just with T-MSD
    save_flag = 1;
    plot_loc(exp_i,type,mi,model_pred, model_predR,full_figure_flag, save_flag)
end



%%
close all;
plot_det_div(exp_i,type,mi,model_pred,full_figure_flag, save_flag)



