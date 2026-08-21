
clear all;
Ntrials = 1000000;
kappa_prior = 0; 
mi = 1;

type = 1; % Det
%type = 2; % Loc
%type=3 -> Det & Loc joint

pars = [4.1887 2.6760 1.8 1.8 1.5 0.56];

Nsamp = 1; % Number of simulated measurements for each trial


J1bar = exp(pars(1:4)); % mean of the gamma distribution of J1
tau   = exp(pars(5)); % scale of the gamma distribution of J1
if ismember(type, [1 3])
    pp    = pars(6);
    nvars = 6;
elseif type == 2
    nvars = 5;
end
% model 2 with decision noise; for type 3 - 2 diff decision noise values
% for det and loc
if mi == 2 & type <3 
    alpha = exp(pars(nvars+1));
    nvars = nvars+1;
elseif mi == 2 & type == 3
    alpha1 = exp(pars(nvars+1));
    alpha2 = exp(pars(nvars+2));
    nvars = nvars+2;   
end


lapse_yes = 0;


% Create part of the data except the responses
% set sizes are only 4
Nvec = [2 3 4 6];
Nind_sel = 3;
Nind = Nind_sel;
for Nti = 1:Ntrials
        stimvec(Nti,1: Nvec(Nind)) = circ_vmrnd(zeros(1,1),kappa_prior,Nvec(Nind)); %* 2;%/pi*90 ;
        targetidx(Nti) = randi(Nvec(Nind));
        
        if type == 1
            if rand <0.5
                targetval(Nti) = stimvec(Nti,targetidx(Nti));
                target_pres(Nti) = 1;
            else
                targetval(Nti) = circ_vmrnd(zeros(1,1),kappa_prior,1);
                target_pres(Nti) = 0;
            end
        elseif type == 2
            targetval(Nti) = stimvec(Nti,targetidx(Nti)); %for loc
        end
        
        vall = [abs(bsxfun(@minus,stimvec(Nti,:)',targetval(Nti)))]'; %btwn (0,2*pi)
        vall(vall>pi) = 2*pi-vall(vall>pi);
        min_ori_diff(Nind,Nti) = min(vall(vall>0)); % btwn 0 and pi
        
        
        vall_mean = circ_mean(setdiff(stimvec(Nti,find(~isnan(stimvec(Nti,:)))), targetval(Nti))');
        vall_mean = abs(bsxfun(@minus,vall_mean, targetval(Nti)));
        vall_mean(vall_mean>pi) = 2*pi- vall_mean(vall_mean>pi);
        mean_ori_diff(Nind,Nti) = vall_mean;

        
        vall_cvar_vec = (setdiff(stimvec(Nti,find(~isnan(stimvec(Nti,:)))), targetval(Nti)));
        vall_cvar = circ_var(vall_cvar_vec,[],[],[]);
        cvar_ori_diff(Nind,Nti) = vall_cvar; % in the interval (0,1)
end

if type==1 | type==3
    
    Ndata_d   = repmat(Nvec, Ntrials, 1);%data_d.N;
    Ntrials_d = Ntrials;
    sdata_d  = stimvec;    
    %sdata_d   = data_d.stims;
    Cdata_d   = target_pres;%randsample([zeros(1,Ntrials/2) ones(1,Ntrials/2)], Ntrials);
    %data_d.target_pres;
    %Ldata_d   = data_d.target_loc; 
    %Resp_d    = data_d.response;
    Nvec_d    = unique(Ndata_d); % only 4
    %target_d  = data_d.target_val;
    target_d = targetval;
    
end 

if type==2 | type==3
    
    Ndata_l   = repmat(Nvec, Ntrials, 1);%data_l.N;
    Ntrials_l = Ntrials; %length(data_l.N);
    sdata_l   = stimvec; %data_l.stims;
    Cdata_l   = ones(1,Ntrials);%data_l.target_pres; %always 1
    Ldata_l   = targetidx';%data_l.target_loc;
    %Resp_l    = data_l.response;
    Nvec_l    = unique(Ndata_l);
    target_l  = targetval';%data_l.target_val;
end

output=[];

if (type==1) | (type==3)
    % Computing probability of responding "target present" on each trial
    ptp = NaN(Ntrials_d,1);
    
    for Nind = Nind_sel%3%4%3%Nvec%1:length(Nvec_d)
        N         = Nvec(Nind); %Nvec_d(Nind);
        idx       = 1:1:Ntrials;%find(Ndata_d==N);
        %J         = gamrnd(J1bar(1)/tau, tau, length(idx), N, Nsamp); %MODIFIED % Precisions: Nidx x N x Nsamp
        J         = gamrnd(J1bar(Nind)/tau, tau, [length(idx), N]);
        kappa     = fisher2kappa(J); % Concentration parameter of von Mises
        
        sdata_de=sdata_d(idx,:)';
        sdata_de=sdata_de(:);
        sdata_de=sdata_de(~isnan(sdata_de));
        sdata_de=reshape(sdata_de,N,length(idx))';
        
        %stims     = repmat(sdata_de(:,1:N),[1,1, Nsamp]); % Read in data and replicate
        stims     = repmat(sdata_de(:,1:N),[1,1]);
        %x         = qrandvm(stims,kappa,[length(idx) N Nsamp]); %circ_vmrnd(stims,kappa);
        x         = qrandvm(stims,kappa,[length(idx) N ]); 
        dec_var = -log(besseli(0,kappa,1)) - kappa + kappa.*cos(bsxfun(@minus,x,target_d(idx)'));
        dec_var_global = squeeze(log(mean(exp(dec_var),2)) + log(pp/(1-pp)));
        
        Resp_d = dec_var_global > 0;
        
        
        if mi == 1
            ptp(idx) = mean(dec_var_global,2); % before it was ptp(idx) = mean(dec_var_global>0,2); NOT OKAY
        elseif mi == 2 & type == 1  
            ptp(idx) = mean(1./(1+ exp(-alpha * dec_var_global)),2);
        elseif mi == 2 & type == 3
            ptp(idx) = mean(1./(1+ exp(-alpha1 * dec_var_global)),2);
        end
        %}
        
        
    end
    
    %{
    ptp(ptp==0) = 1/Nsamp;
    ptp(ptp==1) = 1-1/Nsamp;
    
    output=[output, ptp];
    %}
    output = Resp_d;
end



if (type==2) | (type==3)
    % Computing probability of responding "Lhat" on each trial
    pc = NaN(Ntrials_l,1);
    Ldata_le=NaN(Ntrials_l,1);
    matr=cumsum(isnan(sdata_l),2);
    Ldata_le=Ldata_l-matr(sub2ind(size(matr),[1:1:Ntrials_l]',Ldata_l));
    %Resp_le =Resp_l - matr(sub2ind(size(matr),[1:1:Ntrials_l]',Resp_l));
    
    ddr = NaN(Ntrials_l, Nsamp);
    prob_l = NaN(Ntrials_l,6);
    plhat = NaN(Ntrials_l,1);
    
    for Nind = Nind_sel%3%1:length(Nvec_l)
        N         = Nvec_l(Nind);
        idx       = 1:1:Ntrials;%find(Ndata_l==N);
        
        J         = gamrnd(J1bar(Nind)/tau, tau, length(idx), N, Nsamp); % Precisions: Nidx x max(Nvec_l) x Nsamp
        
        kappa     = fisher2kappa(J); % Concentration parameter of von Mises
        
        sdata_le = sdata_l(idx,:)';
        sdata_le = sdata_le(:);
        sdata_le = sdata_le(~isnan(sdata_le));
        sdata_le = reshape(sdata_le,N,length(idx))';
        
        stims     = repmat(sdata_le(:,1:N),[1,1, Nsamp]); % Read in data -keep all 6 vals, also nan and replicate
        %x         = qrandvm(stims,kappa,[length(idx) N Nsamp]); %circ_vmrnd(stims,kappa);
        x         = qrandvm(stims,kappa,[length(idx) N ]);
        dec_var = -log(besseli(0,kappa,1))-kappa+kappa.*cos(bsxfun(@minus,x,target_l(idx)));
        post_L_given_x = bsxfun(@rdivide, exp(dec_var), sum(exp(dec_var),2));
        if (mi == 1) 
            ddp = NaN(Ntrials_l/4,Nsamp);
            [mx,dd]   = max(post_L_given_x,[],2);%  [mx,dd]   = max(dec_var,[],2); % same
            ddp = squeeze(dd);
            
            
            
            for ti = 1: length(idx)
                for kki = 1: N
                    prob_l(ti,kki) = sum(ddp(ti,:)==kki)/Nsamp;
                end
            end
            
            ddr(idx,:) = ddp;
            %plhat(idx) = prob_l(sub2ind(size(prob_l), [1:1:length(idx)]',Resp_le(idx))); % prob of the chosen response
            pc(idx) =  prob_l(sub2ind(size(prob_l), [1:1:length(idx)]',Ldata_le(idx))); % prob associated with the correct response
            
            Resp_l(idx) = ddr(idx,:);
            
        elseif mi == 2 
            post_L_given_x = bsxfun(@rdivide, exp(dec_var), sum(exp(dec_var),2));
            if type == 2
                post_L_given_x_noisy = post_L_given_x.^alpha;
            elseif type == 3
                post_L_given_x_noisy = post_L_given_x.^alpha2;
            end
            
            post_L_given_x_noisy = bsxfun(@rdivide, post_L_given_x_noisy, sum(post_L_given_x_noisy,2)); %renormalize
            
            for ti = 1: length(idx)
                for si = 1: Nsamp
                    ddr(idx(ti),si)  = find(mnrnd(1,post_L_given_x_noisy(ti,:,si))==1); %generalization of binornd
                end
                
                for kki = 1: N
                    prob_l(ti,kki) = sum(ddr(idx(ti),:)==kki)/Nsamp;
                end
            end
            
            %plhat(idx) = prob_l(sub2ind(size(prob_l), [1:1:length(idx)]',Resp_le(idx))); % prob of the chosen response
            %pc(idx) =  prob_l(sub2ind(size(prob_l), [1:1:length(idx)]',Ldata_le(idx))); % prob associated with the correct response
        end
        
    end
    plhat(plhat==0) = 1/Nsamp;
    plhat(plhat==1) = 1-1/Nsamp;
    
    LL_l = NaN(Ntrials_l,1);
    LL_l = log(plhat);
    
    output=[output, plhat];
end

