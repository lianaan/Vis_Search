function TD_TL_model_fitting(mpr_index)

Nexp = 2;
Nmodels = 2; %1 wo decision noise, 2: w decision noise
Ntypes = 3;
% type = 1  % fit det data only
% type = 2  % fit loc data only
% type = 3  % fit jointly det and loc data

if mpr_index < 1320 % exp 1
    exp_i = 1
    Nsubjects = 11;
else % exp 2
    exp_i = 2
    Nsubjects = 6; 
    
end

if exp_i == 1
    load('alldata_exp1.mat')
elseif exp_i == 2
    load('alldata_exp2.mat')
    mpr_index = mpr_index - 1320; % 2 models *3 types * 11 subjects for exp1 * 20 runs
    % 2 models * 3 types * 6 subjects for exp 2 * 20 runs = 720 % 1320+720 = 2040
end
Nruns = 1200;


% unpack mpr_index
mi = mod(floor(([mpr_index]-0)/(Ntypes*Nsubjects*Nruns)), Nmodels)+1 % model index
type = mod(floor(([mpr_index]-0)/(Nsubjects*Nruns)), Ntypes)+1
par = mod(floor(([mpr_index]-0)/Nruns), Nsubjects)+1 % subject index
runi = mod(([mpr_index]-0),Nruns)+1  % model fitting run index, bads needs several starting points


rng(runi);



addpath('/bads-dev-master/')
options = bads('defaults');              % Default options
options.Ninit = 2;                      % Only 2 points for initial mesh
options.UncertaintyHandling = 1;        % Activate noise handling
options.NoiseSize = 1;                  % Estimated noise magnitude


datac(1)=alldata(par,1).data; % attn, det
datac(2)=alldata(par,2).data; % attn, loc
datac(3)=alldata(par,3).data; % vstm, det
datac(4)=alldata(par,4).data; % vstm, loc


Nvec = unique(datac(1).N);
Ntrials=length(datac(1).N);


if ismember(type, [1 3])
    nvars = 6;
elseif type == 2
    nvars = 5;
end
if mi == 2
    nvars = nvars+1;
end


Jbar_min=log(1.1); %added log to ensure matched parametrization across dimensions
Jbar_max=log(200);
tau_min=log(1.1); 
tau_max=log(200);
pp_min=0.48;
pp_add=0.15; 


PLB = [Jbar_min Jbar_min Jbar_min Jbar_min tau_min ];  % Plausible lower bound
PUB = [Jbar_max Jbar_max Jbar_max Jbar_max tau_max ];  % Plausible upper bound
LB = PLB;    % Lower bound
UB = PUB;    % Upper bound
if ismember(type, [1, 3])
    PLB = [PLB pp_min];  % Plausible lower bound
    PUB = [PUB pp_min+pp_add];  % Plausible upper bound
    LB = [LB 0];      % Lower bound
    UB = [UB 1];
end
if mi == 2
    PLB = [PLB log(0.001)];
    PUB = [PUB log(100)]; %max power to raise posterior
    LB = PLB;   % Lower bound
    UB = PUB;   % Upper bound
end


params_fit=nan(4,nvars); %Only joint model
pred=nan(4, Ntrials);


for cond=[1 3]
    
    
    data_d=datac(cond);
    data_l=datac(cond+1);
    
    start_pars = (PUB-PLB).*rand(1,nvars) + PLB;    % Initial point
    
    
    if type == 1
        fun= @(pars) -sum(sum(Loglike_all(pars, {data_d},mi, type)))
        [params_fit(cond,:), nll(cond), exitflag, output, funValues,gpstruct]=bads(fun,start_pars,LB,UB,PLB,PUB,options);
    
    elseif type == 2
        fun= @(pars) -sum(sum(Loglike_all(pars,  {data_l},mi, type)))
        [params_fit(cond+1,:), nll(cond+1), exitflag, output, funValues,gpstruct]=bads(fun,start_pars,LB,UB,PLB,PUB,options);
                
    elseif type == 3
        fun= @(pars) -sum(sum(Loglike_all(pars, {data_d, data_l},mi, type)))
        [params_fit(cond,:), nll(cond), exitflag, output, funValues,gpstruct]=bads(fun,start_pars,LB,UB,PLB,PUB,options);
                
    end
    
end


savefilename=['model_fits_exp',num2str(exp_i),'_model_', num2str(mi), '_type_', num2str(type), '_sbj_',num2str(par),'_run_',num2str(runi),'.mat'];
save(savefilename,'params_fit','nll','-mat')
