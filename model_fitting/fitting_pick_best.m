clear all; close all;
%exp: 1: stimulus spacing = 60 dva, exp 2: 30 dva
%models: 1: wo decision noise, 2: decision noise
%types: %1: detection, 2: localization, 3: joint

curr_dir = pwd;
model_fits_name = '/model_fits/'

dirname = [curr_dir,model_fits_name];
cd(dirname)
filez = dir([dirname]);
flst = {filez.name};


n_exp = 2;
n_models = 2;
n_types = [3 3];
n_subj = [11 7];

nll_all=nan(n_exp,n_models,max(n_types), max(n_subj),2); % 2 for Perception and Memory
params_all=nan(n_exp,n_models,max(n_types),max(n_subj),2,7);

for ei = 1:2
    load(['alldata_exp',num2str(ei),'.mat']);
    
    
    for mi =  1:2
        for ti = 1:n_types(mi)
            ri_len = nan(n_subj(ei),1);
            for sbjid = 1:n_subj(ei)
                strtofind = ['model_fits_exp', num2str(ei),'_modelE_',num2str(mi), '_type_', num2str(ti),'_sbj_',num2str(sbjid), '_run'];
                ix = regexp(flst,strtofind);
                ix =~ cellfun('isempty',ix);
                ri_len(sbjid) = sum(ix);
            end
            
            
            nll_runs = nan(ri_len(sbjid),2);
            params_runs = nan(ri_len(sbjid),7);
            
            for sbjid = 1: n_subj(ei)
                
                for ri = 1:ri_len(sbjid)
                    
                    load(['model_fits_exp', num2str(ei),'_modelE_',num2str(mi), '_type_', num2str(ti),'_sbj_',num2str(sbjid), '_run_',num2str(ri),'.mat'])
                    
                    if ismember(ti, [1 3])
                        nll_runs(ri,:) = [nll(1) nll(3)];
                        params_runs(ri,1,1:length(params_fit(1,:)))=params_fit(1,:);
                        params_runs(ri,2,1:length(params_fit(3,:)))=params_fit(3,:);
                    elseif ti == 2
                        nll_runs(ri,:) = [nll(2) nll(4)];
                        params_runs(ri,1,1:length(params_fit(2,:)))=params_fit(2,:);
                        params_runs(ri,2,1:length(params_fit(4,:)))=params_fit(4,:);
                    end
                    
                end
                val1=find(squeeze(nll_runs(:,1))==min(nll_runs(:,1))); %Perception condition
                min_r_ind(1)=val1(1);
                val2=find(squeeze(nll_runs(:,2))==min(nll_runs(:,2))); %Memory condition
                min_r_ind(2)=val2(1);
                
                nll_all(ei,mi,ti,sbjid,1)=nll_runs(min_r_ind(1),1);
                nll_all(ei,mi,ti,sbjid,2)=nll_runs(min_r_ind(2),2);
                
                params_all(ei,mi,ti,sbjid,1,1:length(squeeze(params_runs(min_r_ind(1),1,:))))=squeeze(params_runs(min_r_ind(1),1,:));
                params_all(ei,mi,ti,sbjid,2,1:length(squeeze(params_runs(min_r_ind(2),2,:))))=squeeze(params_runs(min_r_ind(2),2,:));
                
            end
            
            
        end
    end
end

%%
savefilename = 'nll_params_best_all.mat'
%save(savefilename, 'nll_all','params_all','-mat')




