function spDist_neural_behav_corr_exclutri_individperm_noflipth_concat(subj,sess,ROIs)

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
 % ROIs = {'V1','V2','V3','V3AB','hV4','LO1','IPS0','IPS1','IPS2','IPS3','sPCS'};
   ROIs = {'V1V2V3','V3AB','hV4','LO1','IPS0IPS1','IPS2IPS3','sPCS'};
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

all_recons_noflip = all_recons{1}; % simply reiterate, we are NOT flipping a thing 

%fprintf('chaotic tpt range grl')
%delay_tpt_range = [3.75 5.25; 7.5 9; 10.5 12];
delay_tpt_range = [0 .75; .75 1.5; 1.5 2.25; 2.25 3; 3 3.75; 3.75 4.5; 4.5 5.25; 5.25 6; 6 6.75; 6.75 7.5; 7.5 8.25; 8.25 9; 9 9.75; 9.75 10.5; 10.5 11.25; 11.25 12;]
myTR = 0.75;

delay_tpts = cell(size(delay_tpt_range,1),1);
for dd = 1:size(delay_tpt_range,1)
    delay_tpts{dd} = (tpts*myTR) >= delay_tpt_range(dd,1) & (tpts*myTR) < delay_tpt_range(dd,2);
end



cond = [1 2]; % do no distractor, then distractor 
store_b = nan(length(cond),length(ROIs),length(subj),36);


figure;

for cc= 1:length(cond)
    
    thisd = nan(length(cond),length(ROIs),length(delay_tpts),length(subj),36, length(angs));
    thisb = nan(1,length(subj));
    
   for vv = 1:length(ROIs) 
       
       for dd =1:length(delay_tpts)
           
        for ss = 1:length(subj)
  
            nruns = 0;
            
            % get one example set of trials
            ru = unique(all_r(all_subj==ss & all_ROIs==vv));
            
            % figure out how many runs we have for that subj, add to total
            nruns = nruns+length(ru);
            for nr =1:nruns
                if cc ==1
                    thisidx = all_data_beh.use_trial==1 &  all_r ==ru(nr) & all_subj==ss & all_ROIs==vv & all_conds(:,1)==1;
                else
                    thisidx = all_data_beh.use_trial==1 &  all_r ==ru(nr) & all_subj==ss & all_ROIs==vv & all_conds(:,1)==2 & all_conds(:,6)==0 ;
                end
                
                thisd(cc,vv,dd,ss,nr,:) = squeeze(mean(mean(all_recons_noflip(thisidx,:,delay_tpts{dd}),1),3));
                
                thisb(ss) = atan2d(sum(mean(mean(all_recons_noflip(thisidx,:, delay_tpts{dd}),1),3).*sind(angs)),...
                    sum(mean(mean(all_recons_noflip(thisidx,:,delay_tpts{dd}),1),3).*cosd(angs)));
                
                store_b(cc,vv,dd,ss,nr) = thisb(ss);

            end
        end
        end
           clear thisidx
           clear thisb
        
    end
 
end

fprintf('neural bias computed ...\n')
%% organize beh data

% which param should we collect?
dist_colors = [ 0.7100 0.2128 0.4772; 0 0 1]; % red - no dist,  blue - dist

params_of_interest = {'f_sacc'};

param_str = {'fsacc'};

mu = nan(2,length(subj),length(cond),36); % no subj performed > 36 runs 
err = nan(2,length(subj),length(cond),36);
dtheta = nan(length(subj),length(cond),36);
drad = nan(length(subj),length(cond),36);
dtheta_poldeg = nan(length(subj),length(cond),36);

for cc = 1:length(cond)
    for ss = 1:length(subj)
        
        nruns = 0;
        
        % get one example set of trials
        ru = unique(all_r(all_subj==ss));
        
        % figure out how many runs we have for that subj, add to total
        nruns = nruns+length(ru);
        
        for nr = 1:nruns
            
            if cc ==1
                thisorigidx = all_data_beh.use_trial==1 &  all_ROIs ==1 & all_r==ru(nr) & all_subj==ss & all_data_beh.s_all.trialinfo(:,1)==1;
            else
                thisorigidx = all_data_beh.use_trial==1 &  all_ROIs ==1 & all_r==ru(nr) & all_subj==ss & all_data_beh.s_all.trialinfo(:,1)==2  & all_data_beh.s_all.trialinfo(:,6)==0;
                
            end
            err(1,ss,cc,nr) = nanstd(all_data_beh.s_all.(params_of_interest{1})(thisorigidx,2));
            mu(:,ss,cc,nr) = nanmean( all_data_beh.s_all.(params_of_interest{1})(thisorigidx,:),  1 ); %mean of a single number = the number .. leaving this
            tmp(1,:) = nanmean( all_data_beh.s_all.(params_of_interest{1})(thisorigidx,:),  1 );
            [dtheta(ss,cc,nr), drad(ss,cc,nr)]= cart2pol(tmp(1),tmp(2));
            dtheta_poldeg(ss,cc,nr) = dtheta(ss,cc,nr).*360/(2*pi);
            
        end
        
        
        
    end
end

% squeeze over subj, collect each subj, all trials into one container
rho_group = nan(length(cond),length(ROIs),1);
rho_group_theta = nan(length(cond),length(ROIs),1);
pval_group = nan(length(cond),length(ROIs),1);
pval_group_theta = nan(length(cond),length(ROIs),1);

% bn_corr =  figure('Name','behneurcorr');
% 
bn_corr_theta_nodist =  figure('Name','behneurcorrthetanodist');
bn_corr_theta_dist =  figure('Name','behneurcorrthetadist');
bn_corr_theta_both =  figure('Name','behneurcorrthetaboth');

% grouprho = figure('Name','grouprho');

for cc = 1:length(cond)
    
for vv=1:length(ROIs)

 for dd =1:length(delay_tpts)   
    thisb_neural_tmp =[];
    thisb_behav_tmp =[];
    thisb_behav = [];
    thisb_neural = [];
    thisb_behav_poldegtmp =[];
    thisb_behav_polradtmp =[];
    thisb_behav_poldeg =[];
    
    for ss = 1:length(subj)   
        thisb_neural_tmp(ss,:) =  squeeze(store_b(cc,vv,dd,ss,:));
        thisb_behav_tmp(ss,:)=  squeeze(mu(2,ss,cc,:))';
        thisb_behav_poldegtmp(ss,:)=  dtheta_poldeg(ss,cc,:);
    end

    thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp)); 
    thisb_neural = thisb_neural_tmp(~isnan(thisb_behav_tmp)); 
    thisb_behav_poldeg = thisb_behav_poldegtmp(~isnan(thisb_behav_tmp));
  

%     set(0, 'CurrentFigure', bn_corr)
%     subplot(length(delay_tpts),length(ROIs),(dd-1)*length(ROIs)+vv); hold on;
%     scatter(thisb_neural,thisb_behav,30,dist_colors(cc,:),'filled','MarkerFaceAlpha',.3)
%     h = lsline;
%     h.LineWidth =2;
%     [rho_group(cc,vv,dd), pval_group(cc,vv,dd)] = corr(thisb_neural, thisb_behav);
%     set(gcf,'Position',[15        1052        2542         286])
%     if vv==1 && cc ==1 && dd ==1
%     ylabel('DVA, sacc endpt')
%     xlabel('Neural bias, Polar angle \circ')
%     else
%     end 
%     clear h 
%     
    if cc ==1 
    set(0, 'CurrentFigure', bn_corr_theta_nodist)
    subplot(length(delay_tpts),length(ROIs),(dd-1)*length(ROIs)+vv); hold on;
    scatter(thisb_neural,thisb_behav_poldeg,30,dist_colors(cc,:),'filled','MarkerFaceAlpha',.3)
    elseif cc ==2
    set(0, 'CurrentFigure', bn_corr_theta_dist)
    subplot(length(delay_tpts),length(ROIs),(dd-1)*length(ROIs)+vv); hold on;
    scatter(thisb_neural,thisb_behav_poldeg,30,dist_colors(cc,:),'filled','MarkerFaceAlpha',.3) 
    end
    q = lsline;
    q.LineWidth =2;
    clear q
    title(ROIs{vv})
    set(gcf,'Position',[  -110         706        1882         624])
    if vv==1 && dd==1
    ylabel('Epoch 1')
    elseif vv==1 && dd==2
    ylabel('Epoch 2')
    elseif vv==1 && dd ==3
    ylabel({'Polar angle \circ delta, target, sacc endpt', 'Epoch 3'})
    xlabel('Neural bias, Polar angle \circ')
    else 
    end
    [rho_group_theta(cc,vv,dd), pval_group_theta(cc,vv,dd)] = corr(thisb_neural, thisb_behav_poldeg);
    
    
%  
    set(0, 'CurrentFigure', bn_corr_theta_both)
    subplot(length(delay_tpts),length(ROIs),(dd-1)*length(ROIs)+vv); hold on;
    scatter(thisb_neural,thisb_behav_poldeg,30,dist_colors(cc,:),'filled','MarkerFaceAlpha',.1)
    title(ROIs{vv})
    lsline
    set(gcf,'Position',[  -110         706        1882         624])
    if vv==1 && dd==1
    ylabel('Epoch 1')
    elseif vv==1 && dd==2
    ylabel('Epoch 1')
    elseif vv==1 && dd ==3
    ylabel({'Polar angle \circ delta, target, sacc endpt', 'Epoch 3'})
    xlabel('Neural bias, Polar angle \circ')
    else 
    end
    [rho_group_theta(cc,vv,dd), pval_group_theta(cc,vv,dd)] = corr(thisb_neural, thisb_behav_poldeg);
 
 end
   
   
% set(0, 'CurrentFigure', grouprho)
% subplot(1,length(ROIs),vv)
% plot(squeeze(rho_group(cc,vv,:)),'-','color',dist_colors(cc,:))
% hold on;
% plot(squeeze(rho_group_theta(cc,vv,:)),'--','color',dist_colors(cc,:))
% legend('rho,linear', 'rho,polar') 

   
   
   
   
end



end

% grouprho = figure('Name','grouprho');
% plot(rho_group,'-')
% hold on;
% plot(rho_group_theta,'--')
% legend('rho,linear', 'rho,polar')

%% corr on individual subj, convert to z, t-test

% fish_fig =  figure('Name','fish');
% 
fish_theta_fig =  figure('Name','fish_theta');
% 
% subjrho = figure('Name','subjrho');
p=[];

pval_sg = nan(length(cond),length(subj),length(ROIs),length(delay_tpts));
rho_sg = nan(length(cond),length(subj),length(ROIs),length(delay_tpts));
fish_z = nan(length(cond),length(subj),length(ROIs),length(delay_tpts));
pval_sg_theta = nan(length(cond),length(subj),length(ROIs),length(delay_tpts));
rho_sg_theta = nan(length(cond),length(subj),length(ROIs),length(delay_tpts));
fish_z_theta = nan(length(cond),length(subj),length(ROIs),length(delay_tpts));

for cc = 1:length(cond)
    
for vv=1:length(ROIs)
    
    for dd = 1:length(delay_tpts)
        
    for ss = 1:length(subj)
        
        thisb_neural_tmp = [];
        thisb_behav_tmp = [];
        thisb_neural = [];
        thisb_behav = [];
        thisb_behav_poldegtmp =[];
        thisb_behav_polradtmp =[];
        thisb_behav_poldeg =[];
        thisb_behav_polrad =[];
        
        thisb_neural_tmp(1,:) =  squeeze(store_b(cc,vv,dd,ss,:));
        thisb_behav_tmp(1,:)=  squeeze(mu(2,ss,cc,:))';
        thisb_behav_poldegtmp(1,:)=  dtheta_poldeg(ss,cc,:);

        
        thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp)); %exclude outliers
        thisb_neural = thisb_neural_tmp(~isnan(thisb_behav_tmp)); %use if no exclu
       
        thisb_behav_poldeg = thisb_behav_poldegtmp(~isnan(thisb_behav_tmp));  
        

        [rho_sg(cc,ss,vv,dd), pval_sg(cc,ss,vv,dd)] = corr(thisb_neural', thisb_behav');

        fish_z(cc,ss,vv,dd)  = atanh(rho_sg(cc,ss,vv,dd));

        [rho_sg_theta(cc,ss,vv,dd), pval_sg_theta(cc,ss,vv,dd)] =corr(thisb_neural', thisb_behav_poldeg');
        
        fish_z_theta(cc,ss,vv,dd)  = atanh(rho_sg_theta(cc,ss,vv,dd)); % identical
   
        clear mysd
    end
%   
%     set(0, 'CurrentFigure', fish_fig)
%     subplot(1,length(ROIs),vv); hold on;
%     plot(dd,fish_z(cc,:,vv,dd),'o','color',dist_colors(cc,:),'markerfacecolor',dist_colors(cc,:),'markersize',1)
%     hold on;
%     my_sem = std(fish_z(cc,:,vv,dd))/sqrt(length(subj));
%     ylim([-.4 .7])
%     errorbar(dd,mean(fish_z(cc,:,vv,dd)),my_sem,'o','color',dist_colors(cc,:),'markerfacecolor',dist_colors(cc,:),'markersize',7)
%     title(ROIs{vv})
%     xlim([0.05 3.5])
    [h,p(cc,vv,dd),~,stats] = ttest(fish_z(cc,:,vv,dd));
    realT(cc,vv,dd) = stats.tstat;
    %ylabel('Fisher-Z')
    
    set(0, 'CurrentFigure', fish_theta_fig)
    subplot(1,length(ROIs),vv); hold on;
    plot(dd,fish_z_theta(cc,:,vv,dd),'o','color',dist_colors(cc,:),'markerfacecolor',dist_colors(cc,:),'markersize',4)
    hold on;
    my_sem_th = std(fish_z_theta(cc,:,vv,dd))/sqrt(length(subj));
    %errorbar(dd,mean(fish_z_theta(cc,:,vv,dd)),my_sem_th,'o','color',dist_colors(cc,:),'markerfacecolor',dist_colors(cc,:),'markersize',7)
    plot(dd,mean(fish_z_theta(cc,:,vv,dd)), 'o', 'color',dist_colors(cc,:),'markerfacecolor',dist_colors(cc,:),'markersize',8)
    plot(dd*[1 1],[mean(fish_z_theta(cc,:,vv,dd))+1.*my_sem_th, mean(fish_z_theta(cc,:,vv,dd))-1.*my_sem_th], '-', 'color',dist_colors(cc,:),'linewidth',1)
    title(ROIs{vv})
    ylim([-1 1])
    xlim([0.05 3.5])
    if vv==1
    set(gca,'Xtick',[0 1 2 3 4],'Xticklabel',{'','Epoch 1','Epoch 2','Epoch 3',''},'XTickLabelRotation',45,'TickDir','out');
    ylabel('Fisher-Z Corr')
    else
    set(gca,'Xtick',[0 1 2 3 4],'Xticklabel',{'','','','',''},'XTickLabelRotation',45,'TickDir','out');
    end 
    
    [h_th,p_th(cc,vv,dd),~,stats_th] = ttest(fish_z_theta(cc,:,vv,dd));
    realT_th(cc,vv,dd) = stats_th.tstat;
    
% 
%     set(0, 'CurrentFigure', subjrho)
%     subplot(length(delay_tpts),length(ROIs),(dd-1)*length(ROIs)+vv); hold on;
%     pl(1)  = plot(1:7,rho_sg(cc,:,vv,dd),'-','color',dist_colors(cc,:))
%     hold on;
%     pl(2) = plot(1:7,rho_sg_theta(cc,:,vv,dd),'--','color',dist_colors(cc,:))
%     if vv==1 && dd==1 && cc ==1
%         title(ROIs{vv})
%     elseif vv==1 && dd ==3 && cc==1
%      ylabel('untransformed Rho')
%       set(gca,'Xtick',[0 1 2 3 4 5 6 7],'Xticklabel',{'','1','2','3','4','5','6','7',''},'XTickLabelRotation',45,'TickDir','out');
%     elseif vv==length(ROIs) && dd==3 && cc==2 
%     legend(pl,'rho,linear', 'rho,polar')
%     end 
%   
     end

end
end 

%%%%%%%%%%%%%%%%%%%%%%%%% PERM TEST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we're permuting two measures here, linear bias and ciruclar bias, open
iter = 10000; 
fish_zperm = nan(length(cond),length(subj),length(ROIs),length(delay_tpts),iter);
fish_zthperm = nan(length(cond),length(subj),length(ROIs),length(delay_tpts),iter);
shuff_T =nan(length(cond),length(ROIs),length(delay_tpts),iter);
shuff_Tth =nan(length(cond),length(ROIs),length(delay_tpts),iter);
rho_perm = nan(length(cond),length(subj),length(ROIs),length(delay_tpts),iter);
pval_perm = nan(length(cond),length(subj),length(ROIs),length(delay_tpts),iter);
rho_permth = nan(length(cond),length(subj),length(ROIs),length(delay_tpts),iter);
pval_permth = nan(length(cond),length(subj),length(ROIs),length(delay_tpts),iter);


for xx = 1:iter
    
    if xx  == 1
        fprintf('doing shuffle ...\n')
        tic
    elseif xx == iter/2
        fprintf('halfway ...\n')
        toc
    else
    end
    for cc = 1:length(cond)
        for vv=1:length(ROIs)
            
            for dd =1:length(delay_tpts)
                
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
                    
                    
                    thisb_neural_tmp_perm (1,:) =  squeeze(store_b(cc,vv,dd,ss,:));
                    thisb_behav_tmp_perm (1,:)=  squeeze(mu(2,ss,cc,:))';
                    
                    thisb_behav_poldegtmp(1,:)=  dtheta_poldeg(ss,cc,:); %this obv cant change across ROI, delay, lose two dimens (vv,dd) here 
                    
                    thisb_behav_perm  = thisb_behav_tmp_perm (~isnan(thisb_behav_tmp_perm) ); %exclude outliers
                    thisb_neural_perm  = thisb_neural_tmp_perm (~isnan(thisb_behav_tmp_perm)); %use if no exclu
                    
                    thisb_behav_poldegperm = thisb_behav_poldegtmp(~isnan(thisb_behav_tmp_perm));
                    
                    thisb_neural_poldegperm  = thisb_neural_tmp_perm (~isnan(thisb_behav_tmp_perm));
                    
                    % do the shuffle
                    
                    shuff_bdata = thisb_behav_perm; %reload these each time
                    shuff_ndata = thisb_neural_perm;
                    subjidx = find(shuff_bdata);
                    shuff_bidx = randperm(length(subjidx))';
                    shuff_bdata(subjidx) = shuff_bdata(shuff_bidx); %insert shuffled subj indices into bdata (note : in this case, bdata is ONLY this subj. just jumbling idx)
                    
                    bdata  = shuff_bdata(subjidx);
                    ndata  = shuff_ndata(subjidx);
                    
                    [rho_perm(cc,ss,vv,dd,xx), pval_perm(cc,ss,vv,dd,xx)] = corr(ndata',bdata'); %1st output is RHO, 2nd is pval
                    fish_zperm(cc,ss,vv,dd,xx) = atanh(rho_perm(cc,ss,vv,dd,xx));
                    
                    shuff_btheta = thisb_behav_poldegperm; %reload these each time
                    shuff_ntheta = thisb_neural_poldegperm;
                    subjidxth = find(shuff_btheta);
                    shuff_bidxth = randperm(length(subjidxth))';
                    shuff_btheta(subjidxth) = shuff_btheta(shuff_bidxth); %insert shuffled subj indices into bdata (note : in this case, bdata is ONLY this subj. just jumbling idx)
                    
                    bthdata  = shuff_btheta(subjidxth);
                    nthdata  = shuff_ntheta(subjidxth);
                    
                    [rho_permth(cc,ss,vv,dd,xx), pval_permth(cc,ss,vv,dd,xx)] = corr(nthdata',bthdata');
                    fish_zthperm(cc,ss,vv,dd,xx) = atanh(rho_permth(cc,ss,vv,dd,xx));
                    
                    clear thisidx subjidx shuff_bidx shuff_nidx shuff_bdata shuff_ndata bdata ndata subjidxth shuff_bidxth  shuff_btheta shuff_ntheta bthdata nthdata
                end
                
                [h,p,~,statsz] = ttest(fish_zperm(cc,:,vv,dd,xx));
                shuff_T(cc,vv,dd,xx) = statsz.tstat;
                
                [h,p,~,statsth] = ttest(fish_zthperm(cc,:,vv,dd,xx));
                shuff_Tth(cc,vv,dd,xx) = statsth.tstat;
            end
        end
    end
end
p_twotail =nan(length(cond),length(ROIs),length(delay_tpts));
p_onetail =nan(length(cond),length(ROIs),length(delay_tpts));
p_th_twotail =nan(length(cond),length(ROIs),length(delay_tpts));
p_th_onetail =nan(length(cond),length(ROIs),length(delay_tpts));
pmask_twotail = [];
pfdr_twotail = [];
pfdr_onetail = [];
pmask_onetail = [];
pfdr_th_twotail =[];
pmask_th_twotail =[];
pmask_th_onetail =[];
pfdr_th_onetail =[];

tdistlinbias = figure('Name','tdistlinbias');
for cc = 1:length(cond)
for dd =1:length(delay_tpts)
for vr = 1:length(ROIs)
    hold on;
    subplot(1,length(ROIs),vr)
    histogram(sort(shuff_T))
    line([realT(cc,vr,dd)  realT(cc,vr,dd)], [0 max(ylim)],'LineWidth',0.75,'color','r')
    p_twotail(cc,vr,dd) = 2 * min ( mean( shuff_T(cc,vr,dd,:) <= realT(vr) ), mean ( shuff_T(cc,vr,dd,:) >= realT(cc,vr,dd)));
    p_onetail(cc,vr,dd) = mean( shuff_T(cc,vr,dd,:) >= realT(cc,vr,dd) );
    title(ROIs{vr})
end

[pfdr_twotail(cc,:,dd), pmask_twotail(cc,:,dd)] =  fdr(p_twotail(cc,:,dd),0.05);
[pfdr_onetail(cc,:,dd), pmask_onetail(cc,:,dd)] =  fdr(p_onetail(cc,:,dd),0.05); 

end
end 



tdistpolbias = figure('Name','tdistpolbias');
for cc =1:length(cond)
for dd = 1:length(delay_tpts)
for vr = 1:length(ROIs)
    
    hold on;
    subplot(1,length(ROIs),vr)
    histogram(sort(shuff_Tth(cc,vr,dd,:)))
    line([realT_th(cc,vr,dd)  realT_th(cc,vr,dd)], [0 max(ylim)],'LineWidth',0.75,'color','r')
    p_th_twotail(cc,vr,dd) = 2 * min ( mean( shuff_Tth(cc,vr,dd,:) <= realT_th(cc,vr,dd) ), mean ( shuff_Tth(cc,vr,dd,:) >= realT_th(cc,vr,dd)));
    p_th_onetail(cc,vr,dd) = mean( shuff_Tth(cc,vr,dd,:) >= realT_th(cc,vr,dd) );
    title(ROIs{vr})
end
  
 
[pfdr_th_twotail(cc,:,dd), pmask_th_twotail(cc,:,dd)] =  fdr(p_th_twotail(cc,:,dd),0.05);
[pfdr_th_onetail(cc,:,dd), pmask_th_onetail(cc,:,dd)] =  fdr(p_th_onetail(cc,:,dd),0.05);
end 
end

fprintf('done')
 
%fn2s = sprintf('%s/%s_corr/neurbehcorr_fsacc_%iter.mat',root,task_dir,iter);
       
%save(fn2s);
        
   
       % fprintf('saving to %s...\n',fn2s);
end

