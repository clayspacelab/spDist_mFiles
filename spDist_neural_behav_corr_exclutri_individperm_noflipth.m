function spDist_neural_behav_corr_exclutri_individperm_noflipth(subj,sess,ROIs)

root = spDist_loadRoot;


task_dir = 'spDist';

if nargin < 1 || isempty(subj)
   subj = {'CC','KD','AY','MR','XL','SF','EK'};
    % subj = {'AY','CC','EK','KD','MR','SF','XL'}; %alph


end

if nargin < 2 || isempty(sess)

    sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
  

end

if nargin < 3 || isempty(ROIs)
  ROIs = {'V1','V2','V3','V3AB','hV4','LO1','IPS0','IPS1','IPS2','IPS3','sPCS'};
   %ROIs = {'V1V2V3','V3AB','IPS0IPS1','IPS2IPS3','sPCS'};
end


func_suffix = 'surf';


nchan = 8;
which_vox = 0.1; % if > 1 , look for WHICH_VOXvox string; otherwise, look for VE<100*WHICH_VOX>

smooth_by = 1; % if this is 1, use regular files, otherwise, load smooth_by files

myTR = 0.75;

t_range_to_plot = [-inf 12]; % plot b/w these (s)

trn_tpts = 7:15; % if blank, load files w/ no _trn%ito%i, otherwise,


% set up file loading strings for below
if smooth_by == 1
    smooth_str = '';
else
    smooth_str = sprintf('_smooth%i',smooth_by);
end

if isempty(trn_tpts)
    trn_str = '';
else
    trn_str = sprintf('_trn%ito%i',trn_tpts(1),trn_tpts(end));
end

if which_vox < 1
    vox_str = sprintf('_VE%03.f',100*which_vox);
else
    vox_str = sprintf('_%ivox',which_vox);
end

% seed random number generator:     
rng(spDist_randSeed);

%% load neural data
startidx = 1;
bidx =1;
WHICH_EXCL = [13 20 21 22]; 
all_data_beh = [];
all_subj_beh=[];

for ss = 1:length(subj)
    
    for vv = 1:length(ROIs)
        
            % just one file to load
            fn = sprintf('%s/%s_reconstructions/%s_%s_%s_%s_%ichan%s%s%s_recon_thruTime1.mat',root,task_dir,subj{ss},horzcat(sess{ss}{:}),ROIs{vv},func_suffix,nchan,vox_str,smooth_str,trn_str);
            
            fprintf('loading %s...\n',fn);
            data = load(fn);
            
              
            if vv == 1 && ss == 1
                % initialize variables...
                
                
                nblankt = length(ROIs)*size(data.recons{1},1);
                all_recons = cell(size(data.recons));
                for aa = 1:length(data.recons)
                    all_recons{aa} = nan(nblankt,size(data.recons{aa},2),size(data.recons{aa},3));
                end
                
                all_recons_nodist = nan(nblankt,size(data.recons_nodist,2),size(data.recons_nodist,3));
                
                all_conds = nan(nblankt,size(data.c_all,2));
                all_angs = nan(nblankt,size(data.a_all,2));
                
                all_fidelity = nan(nblankt,size(data.recons{1},3),length(data.recons)); % timecourse of fidelity for each alignment condition
                all_fidelity_nodist = nan(nblankt,size(data.recons_nodist,3));
                
                all_subj = nan(nblankt,1);
                all_ROIs = nan(nblankt,1);
                all_sess = nan(nblankt,1);
              
                
                angs = data.angs;
                tpts = data.delay_tpts;
                
            end

            thisidx = startidx:(startidx+size(data.c_all,1)-1);
   
            for aa = 1:length(all_recons)
                all_recons{aa}(thisidx,:,:) = data.recons{aa};
                all_fidelity(thisidx,:,aa) = squeeze(mean(cosd(angs) .* data.recons{aa},2));
            end
            
            all_recons_nodist(thisidx,:,:) = data.recons_nodist;
            all_fidelity_nodist(thisidx,:) =  squeeze(mean(cosd(angs) .* data.recons_nodist,2));
            
            all_conds(thisidx,:) = data.c_all;
            all_angs(thisidx,:) = data.a_all;
            
            
            all_subj(thisidx) = ss;
            
            
            all_ROIs(thisidx) = vv;
            
            all_sess(thisidx) = data.sess_all;
            all_r(thisidx,1) = (data.sess_all*100+data.r_all)';
            
            %%%inserting 
            for sessidx = 1:length(sess{ss})
              
                fn = sprintf('%s/spDist_behav_92220/%s_%s_scored.mat',root,subj{ss},sess{ss}{sessidx});
                %fprintf('Loading scored eye data from %s\n',fn);
                this_scored = load(fn);
                this_data.s_all = this_scored.ii_sess;
                
                thisbidx = bidx:(bidx+size(this_scored.ii_sess.trialinfo,1)-1);
                  
                this_data.sess_all = sessidx;
                all_data_beh= cat_struct(all_data_beh,this_data);
                all_subj_beh(thisbidx,1) = ss;
                
                bidx =thisbidx(end)+1;
            end
  
            startidx = thisidx(end)+1;
            
            clear data;
 
    end
    
end
all_data_beh.use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, WHICH_EXCL), all_data_beh.s_all.excl_trial, 'UniformOutput',false));
fprintf('neural and behave loaded...\n')

% visualize all_data_beh.use_trial, make sure the pattern repeats exactly
% across all ROIs
% for ss=1:length(subj)
% figure
% for vv = 1:length(ROIs)
%     hold on;
%     thisidx = all_subj==ss & all_ROIs ==vv;
%     myd = all_data_beh.use_trial(thisidx);
%     subplot(1,length(ROIs),vv)
%     imagesc(myd)
%     clear myd thisidx
% end 
% end
% 


%% load recons, no flipping

all_recons_noflip = all_recons{1}; 

t_range_to_plot_subset = [10.5 12];
tpts_to_plot = (tpts*myTR) >= t_range_to_plot_subset(1) & (tpts*myTR) < t_range_to_plot_subset(2);
fprintf('using delay tpts = %i to %i\n', t_range_to_plot_subset(1), t_range_to_plot_subset(2))

cond_to_plot = [0]; % just one, zero bin
store_b = nan(length(cond_to_plot),length(ROIs),length(subj),36);


figure;

for ff = 1:length(cond_to_plot)
    thisd = nan(length(ROIs),length(subj),36, length(angs));
    
    thisb = nan(1,length(subj));
   for vv = 1:length(ROIs) 
       
        for ss = 1:length(subj)
  
            nruns = 0;
            
            % get one example set of trials
            ru = unique(all_r(all_subj==ss & all_ROIs==vv));
            
            % figure out how many runs we have for that subj, add to total
            nruns = nruns+length(ru);
            for nr =1:nruns
                thisidx = all_data_beh.use_trial==1 &  all_r ==ru(nr) & all_subj==ss & all_ROIs==vv & all_conds(:,1)==2 & all_conds(:,6)==0 ;
                
                %thisidx =  all_r ==ru(nr) & all_subj==ss & all_ROIs==vv & all_conds(:,1)==2 & all_conds(:,6)==0 ;
                
                thisd(vv,ss,nr,:) = squeeze(mean(mean(all_recons_noflip(thisidx,:,tpts_to_plot),1),3));
                
                thisb(ss) = atan2d(sum(mean(mean(all_recons_noflip(thisidx,:,tpts_to_plot),1),3).*sind(angs)),...
                    sum(mean(mean(all_recons_noflip(thisidx,:,tpts_to_plot),1),3).*cosd(angs)));
                
                store_b(ff,vv,ss,nr) = thisb(ss);

            end
        end
        
        
    end
    clear thisidx
    clear thisb
    clear thisd
end

fprintf('neural bias computed ...\n')
%% organize beh data

% which param should we collect?

params_of_interest = {'f_sacc'};
param_str = {'f sacc'};

mu = nan(2,length(subj),36); % no subj performed > 36 runs 
err = nan(2,length(subj),36);
dtheta = nan(length(subj),36);
drad = nan(length(subj),36);
dtheta_poldeg = nan(length(subj),36);


for ss = 1:length(subj)
    
                nruns = 0;
           
                % get one example set of trials
                ru = unique(all_r(all_subj==ss));
                
                % figure out how many runs we have for that subj, add to total
                nruns = nruns+length(ru);
                
                for nr = 1:nruns
 
                       % thisorigidx =  all_ROIs ==1 & all_r==ru(nr) & all_subj_beh==ss & all_data_beh.s_all.trialinfo(:,1)==2  & all_data_beh.s_all.trialinfo(:,6)==0; %now, collect all jitters.
                        
                       thisorigidx = all_data_beh.use_trial==1 &  all_ROIs ==1 & all_r==ru(nr) & all_subj==ss & all_data_beh.s_all.trialinfo(:,1)==2  & all_data_beh.s_all.trialinfo(:,6)==0;
                        
                        err(1,ss,nr) = nanstd(all_data_beh.s_all.(params_of_interest{1})(thisorigidx,2));
                        mu(:,ss,nr) = nanmean( all_data_beh.s_all.(params_of_interest{1})(thisorigidx,:),  1 ); %mean of a single number = the number .. leaving this
                        tmp(1,:) = nanmean( all_data_beh.s_all.(params_of_interest{1})(thisorigidx,:),  1 );
                        [dtheta(ss,nr), drad(ss,nr)]= cart2pol(tmp(1),tmp(2));
                        dtheta_poldeg(ss,nr) = dtheta(ss,nr).*360/(2*pi);
    
                end
                
                
 
end 


% squeeze over subj, collect each subj, all trials into one container
rho_group = nan(length(ROIs),1);
rho_group_theta = nan(length(ROIs),1);
pval_group = nan(length(ROIs),1);
pval_group_theta = nan(length(ROIs),1);

bn_corr =  figure('Name','behneurcorr');

bn_corr_theta =  figure('Name','behneurcorrtheta');

bn_corr_theta2sd =  figure('Name','behneurcorrtheta2sd');

for vv=1:length(ROIs)
    
    thisb_neural_tmp =[];
    thisb_behav_tmp =[];
    thisb_behav = [];
    thisb_neural = [];
    thisb_behav_poldegtmp =[];
    thisb_behav_polradtmp =[];
    thisb_behav_poldeg =[];
   % thisb_behav_polrad =[];
    
    for ss = 1:length(subj)   
        thisb_neural_tmp(ss,:) =  squeeze(store_b(1,vv,ss,:));
        thisb_behav_tmp(ss,:)=  squeeze(mu(2,ss,:))';
        thisb_behav_poldegtmp(ss,:)=  dtheta_poldeg(ss,:);
       % thisb_behav_polradtmp(ss,:)=  dtheta(ss,:);
    end
    %thisidx = ~isnan(thisb_behav_tmp) & thisb_behav_tmp > -8 &  thisb_behav_tmp < 8;
%     
%     thisb_behavt = thisb_behav_tmp(thisidx); %exclude outliers
%     thisb_neuralt = thisb_neural_tmp(thisidx); %use if no exclu
%     thisb_behav_poldegt = thisb_behav_poldegtmp(thisidx);
%     thisb_behav_polradt = thisb_behav_polradtmp(thisidx);
 
% %     thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp) & thisb_behav_tmp > -8 &  thisb_behav_tmp < 8); %exclude outliers
% %     thisb_neural = thisb_neural_tmp(~isnan(thisb_behav_tmp) & thisb_behav_tmp > -8 &  thisb_behav_tmp < 8); %use if no exclu
% %     thisb_behav_poldeg = thisb_behav_poldegtmp(~isnan(thisb_behav_tmp) & thisb_behav_tmp > -8 &  thisb_behav_tmp < 8);
%    % thisb_behav_polrad = thisb_behav_polradtmp(~isnan(thisb_behav_tmp) & thisb_behav_tmp > -8 &  thisb_behav_tmp < 8);
   
    thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp)); %exclude outliers
    thisb_neural = thisb_neural_tmp(~isnan(thisb_behav_tmp)); %use if no exclu
    thisb_behav_poldeg = thisb_behav_poldegtmp(~isnan(thisb_behav_tmp));
   
   
   mysd = std(thisb_behav_poldeg);
   thisb_behav_poldeg2sd = thisb_behav_poldegtmp( ~isnan(thisb_behav_tmp) & abs(thisb_behav_poldegtmp) <= 2 * mysd);
   thisb_neural_poldeg2sd = thisb_neural_tmp( ~isnan(thisb_behav_tmp) & abs(thisb_behav_poldegtmp) <= 2 * mysd);
   
   [rho_group_theta2sd(vv), pval_group_theta2sd(vv)] =corr(thisb_neural_poldeg2sd, thisb_behav_poldeg2sd);
   fish_z_theta2sd(vv)  = atanh(rho_group_theta2sd(vv));
   
   
    set(0, 'CurrentFigure', bn_corr_theta2sd)
    subplot(1,length(ROIs),vv);
    scatter(thisb_neural_poldeg2sd,  thisb_behav_poldeg2sd ,30,'filled','MarkerFaceAlpha',.3)
    h = lsline;
    h.LineWidth =2;
    [rho_group(vv), pval_group(vv)] = corr(thisb_neural, thisb_behav);
    set(gcf,'Position',[15        1052        2542         286])
    if vv==1 
    ylabel('DVA, sacc endpt')
    xlabel('Neural bias, Polar angle \circ')
    else
    end 
    
   
    
    set(0, 'CurrentFigure', bn_corr)
    subplot(1,length(ROIs),vv);
    scatter(thisb_neural,thisb_behav,30,'filled','MarkerFaceAlpha',.3)
    h = lsline;
    h.LineWidth =2;
    [rho_group(vv), pval_group(vv)] = corr(thisb_neural, thisb_behav);
    set(gcf,'Position',[15        1052        2542         286])
    if vv==1 
    ylabel('DVA, sacc endpt')
    xlabel('Neural bias, Polar angle \circ')
    else
    end 
    
    set(0, 'CurrentFigure', bn_corr_theta)
    subplot(1,length(ROIs),vv);
    scatter(thisb_neural,thisb_behav_poldeg,30,'filled','MarkerFaceAlpha',.3)
    h = lsline;
    h.LineWidth =2;
    title(ROIs{vv})
    set(gcf,'Position',[15        1052        2542         286])
    if vv==1
    ylabel('Polar angle \circ delta, target, sacc endpt')
    xlabel('Neural bias, Polar angle \circ')
    else 
    end
    [rho_group_theta(vv), pval_group_theta(vv)] = corr(thisb_neural, thisb_behav_poldeg);

   
end
grouprho = figure('Name','grouprho');
plot(rho_group,'-')
hold on;
plot(rho_group_theta,'--')
legend('rho,linear', 'rho,polar')

%% corr on individual subj, convert to z, t-test

fish_fig =  figure('Name','fish');

fish_theta_fig =  figure('Name','fish_theta');

subjrho = figure('Name','subjrho');
    
pval_sg = nan(length(subj),length(ROIs));
rho_sg = nan(length(subj),length(ROIs));
fish_z = nan(length(subj),length(ROIs));
pval_sg_theta = nan(length(subj),length(ROIs));
rho_sg_theta = nan(length(subj),length(ROIs));
fish_z_theta = nan(length(subj),length(ROIs));

for vv=1:length(ROIs)
    
    for ss = 1:length(subj)
        thisb_neural_tmp = [];
        thisb_behav_tmp = [];
        thisb_neural = [];
        thisb_behav = [];
        thisb_behav_poldegtmp =[];
        thisb_behav_polradtmp =[];
        thisb_behav_poldeg =[];
        thisb_behav_polrad =[];
        
        thisb_neural_tmp(1,:) =  squeeze(store_b(1,vv,ss,:));
        thisb_behav_tmp(1,:)=  squeeze(mu(2,ss,:))';
        thisb_behav_poldegtmp(1,:)=  dtheta_poldeg(ss,:);
%     
%         thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp) & thisb_behav_tmp > -8 &  thisb_behav_tmp < 8); %exclude outliers
%         thisb_neural = thisb_neural_tmp(~isnan(thisb_behav_tmp) & thisb_behav_tmp > -8 &  thisb_behav_tmp < 8); %use if no exclu
%        
%         thisb_behav_poldeg = thisb_behav_poldegtmp(~isnan(thisb_behav_tmp) & thisb_behav_tmp > -8 &  thisb_behav_tmp < 8);  
%         
%         
        
        thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp)); %exclude outliers
        thisb_neural = thisb_neural_tmp(~isnan(thisb_behav_tmp)); %use if no exclu
       
        thisb_behav_poldeg = thisb_behav_poldegtmp(~isnan(thisb_behav_tmp));  
        
        
        mysd = std(thisb_behav_poldeg);
        thisb_behav_poldeg2sd = thisb_behav_poldegtmp( ~isnan(thisb_behav_tmp) & abs(thisb_behav_poldegtmp) <= 2 * mysd); 
        thisb_neural_poldeg2sd = thisb_neural_tmp( ~isnan(thisb_behav_tmp) & abs(thisb_behav_poldegtmp) <= 2 * mysd);

        [rho_sg_theta2sd(ss,vv), pval_sg_theta2sd(ss,vv)] =corr(thisb_neural_poldeg2sd', thisb_behav_poldeg2sd');
        fish_z_theta2sd(ss,vv)  = atanh(rho_sg_theta2sd(ss,vv));
        
        
        [rho_sg(ss,vv), pval_sg(ss,vv)] = corr(thisb_neural', thisb_behav');

        fish_z(ss,vv)  = atanh(rho_sg(ss,vv));

        [rho_sg_theta(ss,vv), pval_sg_theta(ss,vv)] =corr(thisb_neural', thisb_behav_poldeg');
        
        fish_z_theta(ss,vv)  = atanh(rho_sg_theta(ss,vv)); % identical
   
        clear mysd
    end
    
    set(0, 'CurrentFigure', fish_fig)
    subplot(1,length(ROIs),vv)
    plot(1,fish_z(:,vv),'bo','markerfacecolor','b','markersize',5)
    hold on;
    my_sem = std(fish_z(:,vv))/sqrt(length(subj));
    ylim([-.4 .7])
    errorbar(1,mean(fish_z(:,vv)),my_sem,'ko','markerfacecolor','k','markersize',10)
    title(ROIs{vv})
    
    [h,p(vv),~,stats] = ttest(fish_z(:,vv));
    realT(vv) = stats.tstat;
    ylabel('Fisher-Z')
    
    set(0, 'CurrentFigure', fish_theta_fig)
    subplot(1,length(ROIs),vv)
    plot(1,fish_z_theta(:,vv),'bo','markerfacecolor','b','markersize',5)
    hold on;
    my_sem_th = std(fish_z_theta(:,vv))/sqrt(length(subj));
    errorbar(1,mean(fish_z_theta(:,vv)),my_sem_th,'ko','markerfacecolor','k','markersize',10)
    title(ROIs{vv})
    ylim([-.4 .7])
    [h_th,p_th(vv),~,stats_th] = ttest(fish_z_theta(:,vv));
    realT_th(vv) = stats_th.tstat;
    

    set(0, 'CurrentFigure', subjrho)
    subplot(1,length(ROIs),vv)
    plot(1:7,rho_sg(:,vv),'-')
    hold on;
    plot(1:7,rho_sg_theta(:,vv),'--')
    if vv== length(ROIs)
    legend('rho,linear', 'rho,polar')
    else
    end 
  
    

end


%%%%%%%%%%%%%%%%%%%%%%%%% PERM TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we're permuting two measures here, linear bias and ciruclar bias, open
% all variables up front 
iter = 10000; 
fish_zperm = nan(length(subj),length(ROIs),iter);
fish_zthperm = nan(length(subj),length(ROIs),iter);
shuff_T =nan(length(ROIs),iter);
shuff_Tth =nan(length(ROIs),iter);
rho_perm = nan(length(subj),length(ROIs),iter);
pval_perm = nan(length(subj),length(ROIs),iter);
rho_permth = nan(length(subj),length(ROIs),iter);
pval_permth = nan(length(subj),length(ROIs),iter);


for xx = 1:iter
    
    if xx  == 1
        fprintf('doing shuffle ...\n')
        tic
    elseif xx == iter/2
        fprintf('halfway ...\n')
        toc
    else
    end
    
    for vv=1:length(ROIs)
        
        for ss = 1:length(subj)
            thisb_neural_tmp_perm = [];
            thisb_behav_tmp_perm  = [];
            thisb_behav_perm  = [];
            thisb_neural_perm  = [];
            shuff_bdata = [];
            shuff_ndata = [];
            thisb_behav_poldegtmp =[];
            thisb_behav_poldegperm =[];
            thisb_neural_poldegperm =[];
            shuff_btheta =[];
            shuff_ntheta =[];
            
            
            thisb_neural_tmp_perm (1,:) =  squeeze(store_b(1,vv,ss,:));
            thisb_behav_tmp_perm (1,:)=  squeeze(mu(2,ss,:))';
            
            thisb_behav_poldegtmp(1,:)=  dtheta_poldeg(ss,:);
            %thisb_behav_polradtmp(1,:)=  dtheta(ss,:);
            
%             thisb_behav_perm  = thisb_behav_tmp_perm (~isnan(thisb_behav_tmp_perm) & thisb_behav_tmp_perm  > -8 &  thisb_behav_tmp_perm  < 8); %exclude outliers
%             thisb_neural_perm  = thisb_neural_tmp_perm (~isnan(thisb_behav_tmp_perm) & thisb_behav_tmp_perm  > -8 &  thisb_behav_tmp_perm  < 8); %use if no exclu
%             
%             thisb_behav_poldegperm = thisb_behav_poldegtmp(~isnan(thisb_behav_tmp_perm) & thisb_behav_tmp_perm > -8 &  thisb_behav_tmp_perm < 8);
%             thisb_neural_poldegperm  = thisb_neural_tmp_perm (~isnan(thisb_behav_tmp_perm) & thisb_behav_tmp_perm  > -8 &  thisb_behav_tmp_perm  < 8);


            thisb_behav_perm  = thisb_behav_tmp_perm (~isnan(thisb_behav_tmp_perm) ); %exclude outliers
            thisb_neural_perm  = thisb_neural_tmp_perm (~isnan(thisb_behav_tmp_perm)); %use if no exclu
            
            thisb_behav_poldegperm = thisb_behav_poldegtmp(~isnan(thisb_behav_tmp_perm));



            %thisb_behav_poldegperm = thisb_behav_poldegtmp(~isnan(thisb_behav_poldegtmp));
           % thisb_neural_poldegperm =   thisb_neural_tmp_perm(~isnan(thisb_behav_poldegtmp));
            thisb_neural_poldegperm  = thisb_neural_tmp_perm (~isnan(thisb_behav_tmp_perm));
            
            % do the shuffle
            
            shuff_bdata = thisb_behav_perm; %reload these each time
            shuff_ndata = thisb_neural_perm;
            subjidx = find(shuff_bdata);
            shuff_bidx = randperm(length(subjidx))';
            shuff_bdata(subjidx) = shuff_bdata(shuff_bidx); %insert shuffled subj indices into bdata (note : in this case, bdata is ONLY this subj. just jumbling idx)
            
            bdata  = shuff_bdata(subjidx);
            ndata  = shuff_ndata(subjidx);
            
            [rho_perm(ss,vv,xx), pval_perm(ss,vv,xx)] = corr(ndata',bdata'); %1st output is RHO, 2nd is pval
            fish_zperm(ss,vv,xx) = atanh(rho_perm(ss,vv,xx));
            
            shuff_btheta = thisb_behav_poldegperm; %reload these each time
            shuff_ntheta = thisb_neural_poldegperm;
            subjidxth = find(shuff_btheta);
            shuff_bidxth = randperm(length(subjidxth))';
            shuff_btheta(subjidxth) = shuff_btheta(shuff_bidxth); %insert shuffled subj indices into bdata (note : in this case, bdata is ONLY this subj. just jumbling idx)
            
            bthdata  = shuff_btheta(subjidxth);
            nthdata  = shuff_ntheta(subjidxth);
            
            [rho_permth(ss,vv,xx), pval_permth(ss,vv,xx)] = corr(nthdata',bthdata');
            fish_zthperm(ss,vv,xx) = atanh(rho_permth(ss,vv,xx));
            
            clear thisidx subjidx shuff_bidx shuff_nidx shuff_bdata shuff_ndata bdata ndata subjidxth shuff_bidxth  shuff_btheta shuff_ntheta bthdata nthdata
        end
        
        [h,p,~,statsz] = ttest(fish_zperm(:,vv,xx));
        shuff_T(vv,xx) = statsz.tstat;
        
        [h,p,~,statsth] = ttest(fish_zthperm(:,vv,xx));
        shuff_Tth(vv,xx) = statsth.tstat;
        
    end
end

p_twotail =nan(length(ROIs),1);
p_onetail =nan(length(ROIs),1);
p_th_twotail =nan(length(ROIs),1);
p_th_onetail =nan(length(ROIs),1);
pmask_twotail = [];
pfdr_twotail = [];
pfdr_onetail = [];
pmask_onetail = [];
pfdr_th_twotail =[];
pmask_th_twotail =[];
pmask_th_onetail =[];
pfdr_th_onetail =[];

tdistlinbias = figure('Name','tdistlinbias');
for vr = 1:length(ROIs)
    hold on;
    subplot(1,length(ROIs),vr)
    histogram(sort(shuff_T))
    line([realT(vr)  realT(vr)], [0 max(ylim)],'LineWidth',0.75,'color','r')
    p_twotail(vr) = 2 * min ( mean( shuff_T(vr,:) <= realT(vr) ), mean ( shuff_T(vr,:) >= realT(vr)));
    p_onetail(vr) = mean( shuff_T(vr,:) >= realT(vr) );
    title(ROIs{vr})
end
  
 
[pfdr_twotail, pmask_twotail] =  fdr(p_twotail(:),0.05);
[pfdr_onetail, pmask_onetail] =  fdr(p_onetail(:),0.05); 


tdistpolbias = figure('Name','tdistpolbias');
for vr = 1:length(ROIs)
    
    hold on;
    subplot(1,length(ROIs),vr)
    histogram(sort(shuff_Tth))
    line([realT_th(vr)  realT_th(vr)], [0 max(ylim)],'LineWidth',0.75,'color','r')
    p_th_twotail(vr) = 2 * min ( mean( shuff_Tth(vr,:) <= realT_th(vr) ), mean ( shuff_Tth(vr,:) >= realT_th(vr)));
    p_th_onetail(vr) = mean( shuff_Tth(vr,:) >= realT_th(vr) );
    title(ROIs{vr})
end
  
 
[pfdr_th_twotail, pmask_th_twotail] =  fdr(p_th_twotail(:),0.05);
[pfdr_th_onetail, pmask_th_onetail] =  fdr(p_th_onetail(:),0.05);

% one last plot, t-dists overlaid

tdistlinpolbias = figure('Name','tdistlinpolbias');

for vr = 1:length(ROIs)
    hold on;
    subplot(1,length(ROIs),vr)
    histogram(sort(shuff_T))
    line([realT(vr)  realT(vr)], [0 max(ylim)],'LineWidth',0.75,'color','r')
    hold on; 
    histogram(sort(shuff_Tth))
    line([realT_th(vr)  realT_th(vr)], [0 max(ylim)],'LineWidth',0.75,'color','b')

    title(ROIs{vr})
    if vr ==length(ROIs)
      legend('lin tval', 'pol tval')
    else
    end
    
end
  

fprintf('done')
end

