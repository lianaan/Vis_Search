function [output ddr pc] = Predict_all(pars, data,mi, type)

%type=1 -> Det
%type=2 -> Loc
%type=3 -> Det & Loc joint


Nsamp = 1200;%1200;%2000; % Number of simulated measurements for each trial



if type==1
    data_d=data{1};
elseif type==2
    data_l=data{1};
elseif type==3
    data_d=data{1};
    data_l=data{2};
end


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



% Unpack data

if type==1 | type==3
    
    Ndata_d   = data_d.N;
    Ntrials_d = length(data_d.N);
    sdata_d   = data_d.stims;
    Cdata_d   = data_d.target_pres;
    Ldata_d   = data_d.target_loc; 
    Resp_d    = data_d.response;
    Nvec_d    = unique(Ndata_d);
    target_d  = data_d.target_val;
end
if type==2 |type==3
    
    Ndata_l   = data_l.N;
    Ntrials_l = length(data_l.N);
    sdata_l   = data_l.stims;
    Cdata_l   = data_l.target_pres; %always 1
    Ldata_l   = data_l.target_loc;
    Resp_l    = data_l.response;
    Nvec_l    = unique(Ndata_l);
    target_l  = data_l.target_val;
end


output=[];

if (type==1) | (type==3)
    % Computing probability of responding "target present" on each trial
    ptp = NaN(Ntrials_d,1);
    
    for Nind = 1:length(Nvec_d)
        N         = Nvec_d(Nind);
        idx       = find(Ndata_d==N);
        J         = gamrnd(J1bar(Nind)/tau, tau, length(idx), N, Nsamp); % Precisions: Nidx x N x Nsamp
        kappa     = fisher2kappa(J); % Concentration parameter of von Mises
        
        sdata_de=sdata_d(idx,:)';
        sdata_de=sdata_de(:);
        sdata_de=sdata_de(~isnan(sdata_de));
        sdata_de=reshape(sdata_de,N,length(idx))';
        
        stims     = repmat(sdata_de(:,1:N),[1,1, Nsamp]); % Read in data and replicate
        x         = qrandvm(stims,kappa,[length(idx) N Nsamp]); %circ_vmrnd(stims,kappa);
        dec_var = -log(besseli(0,kappa,1)) - kappa + kappa.*cos(bsxfun(@minus,x,target_d(idx)));
        dec_var_global = squeeze(log(mean(exp(dec_var),2)) + log(pp/(1-pp)));
        if mi == 1
            ptp(idx) = mean(dec_var_global>0,2);
        elseif mi == 2 & type == 1  
            ptp(idx) = mean(1./(1+ exp(-alpha * dec_var_global)),2);
        elseif mi == 2 & type == 3
            ptp(idx) = mean(1./(1+ exp(-alpha1 * dec_var_global)),2);
        end
        
        
    end
    
    
    ptp(ptp==0) = 1/Nsamp;
    ptp(ptp==1) = 1-1/Nsamp;
    
    output=[output, ptp];
end





if (type==2) | (type==3)
    % Computing probability of responding "Lhat" on each trial
    pc = NaN(Ntrials_l,1);
    Ldata_le=NaN(Ntrials_l,1);
    matr=cumsum(isnan(sdata_l),2);
    Ldata_le=Ldata_l-matr(sub2ind(size(matr),[1:1:Ntrials_l]',Ldata_l));
    Resp_le =Resp_l - matr(sub2ind(size(matr),[1:1:Ntrials_l]',Resp_l));
    
    ddr = NaN(Ntrials_l, Nsamp);
    prob_l = NaN(Ntrials_l,6);
    plhat = NaN(Ntrials_l,1);
    
    for Nind = 1:length(Nvec_l)
        N         = Nvec_l(Nind);
        idx       = find(Ndata_l==N);
        
        J         = gamrnd(J1bar(Nind)/tau, tau, length(idx), N, Nsamp); % Precisions: Nidx x max(Nvec_l) x Nsamp
        
        kappa     = fisher2kappa(J); % Concentration parameter of von Mises
        
        sdata_le = sdata_l(idx,:)';
        sdata_le = sdata_le(:);
        sdata_le = sdata_le(~isnan(sdata_le));
        sdata_le = reshape(sdata_le,N,length(idx))';
        
        stims     = repmat(sdata_le(:,1:N),[1,1, Nsamp]); % Read in data -keep all 6 vals, also nan and replicate
        x         = qrandvm(stims,kappa,[length(idx) N Nsamp]); %circ_vmrnd(stims,kappa);
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
            plhat(idx) = prob_l(sub2ind(size(prob_l), [1:1:length(idx)]',Resp_le(idx))); % prob of the chosen response
            pc(idx) =  prob_l(sub2ind(size(prob_l), [1:1:length(idx)]',Ldata_le(idx))); % prob associated with the correct response
            
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
            
            plhat(idx) = prob_l(sub2ind(size(prob_l), [1:1:length(idx)]',Resp_le(idx))); % prob of the chosen response
            pc(idx) =  prob_l(sub2ind(size(prob_l), [1:1:length(idx)]',Ldata_le(idx))); % prob associated with the correct response
        end
        
    end
    plhat(plhat==0) = 1/Nsamp;
    plhat(plhat==1) = 1-1/Nsamp;
    
    LL_l = NaN(Ntrials_l,1);
    LL_l = log(plhat);
    
    output=[output, plhat];
end


end
