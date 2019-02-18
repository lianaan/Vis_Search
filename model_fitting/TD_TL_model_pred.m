function TD_TL_model_pred(mpr_index)
%params_all(ei,mi,ti,sbjid,1,1:length(squeeze(params_runs(min_r_ind(1),1,:))))

%exp_i, type
Nexp = 2;
Nmodels = 2; %1 wo decision noise, 2: w decision noise
Ntypes = 3;
% type = 1  % fit det data only
% type = 2  % fit loc data only
% type = 3  % fit jointly det and loc data

if mpr_index < 66 % exp 1
    exp_i = 1;
    Nsubjects = 11;
else % exp 2
    exp_i = 2;
    Nsubjects = 7;
    mpr_index = mpr_index - 66; % 2 models *3 types * 11 subjects for exp1
    % 2 models * 3 types * 7 subjects for exp 2  = 42
end

load(['alldata_exp',num2str(exp_i),'.mat'])
%load('nll_params_best_all.mat')
load('nll_params_best_all_miE_2_type_3.mat')

mi = mod(floor(([mpr_index]-0)/(Ntypes*Nsubjects)), Nmodels)+1 % model index
type = mod(floor(([mpr_index]-0)/(Nsubjects)), Ntypes)+1
par = mod(([mpr_index]-0), Nsubjects)+1 % subject index


datac(1)=alldata(par,1).data; % attn, det
datac(2)=alldata(par,2).data; % attn, loc
datac(3)=alldata(par,3).data; % vstm, det
datac(4)=alldata(par,4).data; % vstm, loc


Nvec = unique(datac(1).N);
Ntrials=length(datac(1).N);

params_fit = nan(4,8); % changed 12/19/2018
pred = nan(4, Ntrials);
for cond=[1 3]
    
    
    data_d=datac(cond);
    data_l=datac(cond+1);
    
    
    
    if type == 1
        params_fit(cond,:) = squeeze(params_all(exp_i,mi,type,par,(cond+1)/2,:));
        pred(cond,:)=Predict_all(params_fit(cond,:), {data_d},mi, type)';
    elseif type == 2
        params_fit(cond+1,:) = squeeze(params_all(exp_i,mi,type,par,(cond+1)/2,:));
        pred(cond+1,:)=Predict_all(params_fit(cond+1,:), {data_l},mi, type)'; 
    elseif type == 3
        params_fit(cond,:) = squeeze(params_all(exp_i,mi,type,par,(cond+1)/2,:));
        pred(cond:cond+1,:)=Predict_all(params_fit(cond,:), {data_d, data_l},mi, type)';  
    end
    
end


savefilename=['model_pred_exp_',num2str(exp_i),'_model_', num2str(mi), '_type_', num2str(type), '_sbj_',num2str(par),'.mat'];
save(savefilename,'pred','-mat')
