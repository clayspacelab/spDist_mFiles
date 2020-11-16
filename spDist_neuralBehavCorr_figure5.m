function spDist_neuralBehavCorr_figure5(subj,sess,ROIs)

root = spDist_loadRoot;


task_dir = 'spDist';

if nargin < 1 || isempty(subj)
   subj = {'AY','CC','EK','KD','MR','SF','XL'}; %alph


end

if nargin < 2 || isempty(sess)

    sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
  

end

if nargin < 3 || isempty(ROIs)
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

% seed random number generator 

rng(spDist_randSeed);

%% load neural data
startidx = 1;
bidx =1;
WHICH_EXCL = [13 20 21 22]; 
all_data_beh = [];
all_subj_beh= [];

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


%% load recons, no flipping

all_recons_noflip = all_recons{1}; % simply reiterate, we are NOT flipping a thing 


delay_tpt_range = [3.75 5.25; 8.25 9.75; 10.5 12]; %updated on nov 092020
myTR = 0.75;

delay_tpts = cell(size(delay_tpt_range,1),1);
for dd = 1:size(delay_tpt_range,1)
    delay_tpts{dd} = (tpts*myTR) >= delay_tpt_range(dd,1) & (tpts*myTR) < delay_tpt_range(dd,2);
end



cond = [1 2]; % do no distractor, then distractor 
store_b = nan(length(cond),length(ROIs),length(delay_tpts),length(subj),36);


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
                if cc ==1 %%%%%% NOTE THIS COLLECTS THREE TRIALS IN MGS CONDITION %%%% NEED TO FIX!!!!
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
dist_colors = spDist_condColors; 

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
            err(1,ss,cc,nr) = std(all_data_beh.s_all.(params_of_interest{1})(thisorigidx,2));
            mu(:,ss,cc,nr) = mean( all_data_beh.s_all.(params_of_interest{1})(thisorigidx,:),  1 ); %mean of a single number = the number .. leaving this
            tmp(1,:) = mean( all_data_beh.s_all.(params_of_interest{1})(thisorigidx,:),  1 );
            [dtheta(ss,cc,nr), drad(ss,cc,nr)]= cart2pol(tmp(1),tmp(2));
            dtheta_poldeg(ss,cc,nr) = dtheta(ss,cc,nr).*360/(2*pi);
            
        end
        
        
        
    end
end

%% bias  %% Figure 1F : Near Distractor trials: Bias, polar degrees 
 
params_of_interest = {'f_sacc'};
param_str = {'final sacc'};

mu_bias = nan(1,1,2,length(subj)); % no subj performed > 36 runs 

drad_bias = [];
dtheta_bias = [];

for pp =1:length(params_of_interest)
    for bb =1
        for ss = 1:length(subj)
            
            tmpd =  all_data_beh.s_all.(params_of_interest{pp});
            tmpidx = all_subj==ss & all_data_beh.s_all.trialinfo(:,1)==2 & all_data_beh.use_trial==1 & all_data_beh.s_all.trialinfo(:,6)==0 & all_data_beh.s_all.trialinfo(:,10) < 0; %negative CW jitter
            orig_y =  tmpd(tmpidx,2);
            orig_y_flip = orig_y*-1;
            tmpd(tmpidx,2) = orig_y_flip;  % here is where the actual y flip is inserted for CW bins, into tmpd ONLY
            
            thisidx = all_subj==ss & all_data_beh.s_all.trialinfo(:,1)==2 & all_data_beh.use_trial==1 & all_data_beh.s_all.trialinfo(:,6)==0; %now, collect all jitters.
            
            % distractor bin x param x [radial; tangential] x subj
            tmp_mu(bb,pp,:,ss) = mean( tmpd(thisidx,:),  1 ); % use idx with all zero bin from flipped data
          
            [dtheta_bias, drad_bias] = cart2pol(tmpd(thisidx,1), tmpd(thisidx,2));
            mu_bias(bb,pp,:,ss) = mean([drad_bias rad2deg(dtheta_bias)]); %put rad as x FIRST, theta as Y second, in keeping with later indexing
            
            clear dtheta_err drad_err drad_mu dtheta_mu tmpd thisidx orig_y orig_y_flip
            
            
            
        end
        
    end
end

figure;

for pp = 1:length(params_of_interest)
    
    subplot(1,length(params_of_interest),pp); hold on;
    
    % nbins x nsubj
    thism_near = squeeze(mu_bias(:,pp,2,:))'; %this is the only condition that is information for this analysis.
    plot(1+0.15,thism_near,'o','MarkerSize',5,'Markerfacecolor',[0.3 0.3 0.3], 'Color',[0.3 0.3 0.3]);
    plot(1,mean(thism_near,2),'o','MarkerSize',15,'Color',dist_colors(2,:),'MarkerFaceColor',dist_colors(2,:));
    hold on;
    
    tmpe = std(thism_near)/sqrt(length(subj));
    plot([1 1], mean(thism_near)+tmpe*[-1 1],'-','LineWidth',1.5,'Color',dist_colors(2,:));

    %errorbar(1, mean(thism_near), [std(thism_near)/sqrt(length(subj))],'linewidth',2,'Color',cond_colors(2,:))
    title(param_str{pp});
    
   
    ylabel('Bias, Polar Angle (\circ), toward distractor');
end
xlabel('Near Distractor')
xlim([0 2])

[h_bias_near p_bias_near,CI,STATS] = ttest(thism_near(1,:)') %t-test to determine if bias is different than zero
text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',p_bias_near),'color','k','fontsize',15) %filled

match_ylim(get(gcf,'Children'));
set(gcf,'position',[ 549   724   499   571])
%% Fgiure 5 B & C

rho_group_theta = nan(length(cond),length(ROIs),1);
pval_group_theta = nan(length(cond),length(ROIs),1);
rho_theta = nan(length(subj),length(ROIs),1);
pval_theta = nan(length(subj),length(ROIs),1);
subj_col = lines(7);

figure 

for ss = 1:length(subj)
    hold on;
    for cc = 2 % 1:length(cond)
        
        for vv = 1:length(ROIs)
            
            for dd = 3 % 1:length(delay_tpts)
                
                thisb_neural_tmp =[];
                thisb_behav_tmp =[];
                thisb_behav = [];
                thisb_neural = [];
                thisb_behav_poldegtmp =[];
                thisb_behav_poldeg =[];
                
                
                thisb_neural_tmp(1,:) =  squeeze(store_b(cc,vv,dd,ss,:));
                thisb_behav_tmp(1,:)=  squeeze(mu(2,ss,cc,:))';
                thisb_behav_poldegtmp(1,:)=  dtheta_poldeg(ss,cc,:);
                
                
                thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp));
                thisb_neural = thisb_neural_tmp(~isnan(thisb_behav_tmp));
                thisb_behav_poldeg = thisb_behav_poldegtmp(~isnan(thisb_behav_tmp));
              
                subplot(1,length(ROIs),vv); hold on;
             
                scatter(thisb_neural,thisb_behav_poldeg,30,subj_col(ss,:),'filled','MarkerFaceAlpha',.4)
                
                title(ROIs{vv})
                lsline
                set(gcf,'Position',[  -110         706        1882         624])
         
                if vv==1 && ss ==7 
                    ylabel({'Polar angle \circ delta, target, sacc endpt', 'Epoch 3'})
                    xlabel('Neural bias, Polar angle \circ')
                else
                end
                [rho_theta(cc,ss,vv,dd), pval_theta(cc,ss,vv,dd)] = corr(thisb_neural', thisb_behav_poldeg');
                clear this neural thisb_behav_poldeg
            end
            
            
        end
  
    end
   
    
    
end

 set(get(gcf,'Children'),'XLim',[-180 180],'YLim',[-25 25],'xtick',-180:90:180) %,'xtick',{'180','0','180'},

 %% make group subplot - supp fig 5 A
% 
% % squeeze over subj, collect each subj, all trials into one container
% rho_group_theta = nan(length(cond),length(ROIs),length(delay_tpts));
% pval_group_theta = nan(length(cond),length(ROIs),length(delay_tpts));
% 
% subj_col = lines(7);
% 
% figure 
% 
% 
%     hold on;
%         
%     for vv = 1:length(ROIs)
%         for cc = 2 % 1:length(cond)
%             
%             for dd = 1:length(delay_tpts)
%                 
%                 thisb_neural_tmp =[];
%                 thisb_behav_tmp =[];
%                 thisb_behav = [];
%                 thisb_neural = [];
%                 thisb_behav_poldegtmp =[];
%                 thisb_behav_poldeg =[];
%                 
%                 for ss = 1:length(subj)
%                 thisb_neural_tmp(ss,:) =  squeeze(store_b(cc,vv,dd,ss,:));
%                 thisb_behav_tmp(ss,:)=  squeeze(mu(2,ss,cc,:))';
%                 thisb_behav_poldegtmp(ss,:)=  dtheta_poldeg(ss,cc,:);
%                 
%                 end 
%                 thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp));
%                 thisb_neural = thisb_neural_tmp(~isnan(thisb_behav_tmp));
%                 thisb_behav_poldeg = thisb_behav_poldegtmp(~isnan(thisb_behav_tmp));
%     
%                 subplot(length(delay_tpts),length(ROIs),(dd-1)*length(ROIs)+vv); hold on;
%                 scatter(thisb_neural,thisb_behav_poldeg,30,dist_colors(cc,:),'filled','MarkerFaceAlpha',.4)
%                 title(ROIs{vv})
%                 lsline
%                 set(gcf,'Position',[  -110         706        1882         624])
%                
%                
%                 if vv==1 && ss ==1 
%                     ylabel({'Polar angle \circ delta, target, sacc endpt', 'Epoch 3'})
%                     xlabel('Neural bias, Polar angle \circ')
%                 else
%                 end
%                 [rho_group_theta(cc,vv,dd), pval_group_theta(cc,vv,dd)] = corr(thisb_neural, thisb_behav_poldeg);
%                 clear this neural thisb_behav_poldeg
%             end
%             
%             
%         end
%         
%         
%         
%         
%     end
%    
%     
%  set(get(gcf,'Children'),'XLim',[-180 180],'YLim',[-25 25],'xtick',-180:90:180)

 %% supp fig 5 c
% 
% store_fish_ztheta= [];
% pval_sg_theta = nan(length(cond),length(subj),length(ROIs),length(delay_tpts));
% rho_sg_theta = nan(length(cond),length(subj),length(ROIs),length(delay_tpts));
% fish_z_theta = nan(length(cond),length(subj),length(ROIs),length(delay_tpts));
% 
% figure
% for cc = 2 %1:length(cond)
%     
% for vv=1:length(ROIs)
%     
%     for dd = 3 %1:length(delay_tpts)
%         
%     for ss = 1:length(subj)
%         
%         thisb_neural_tmp = [];
%         thisb_behav_tmp = [];
%         thisb_neural = [];
%         thisb_behav = [];
%         thisb_behav_poldegtmp =[];
%         thisb_behav_polradtmp =[];
%         thisb_behav_poldeg =[];
%         thisb_behav_polrad =[];
%         
%         thisb_neural_tmp(1,:) =  squeeze(store_b(cc,vv,dd,ss,:));
%         thisb_behav_tmp(1,:)=  squeeze(mu(2,ss,cc,:))';
%         thisb_behav_poldegtmp(1,:)=  dtheta_poldeg(ss,cc,:);
%         
%         
%         thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp)); 
%         thisb_neural = thisb_neural_tmp(~isnan(thisb_behav_tmp));
%        
%         thisb_behav_poldeg = thisb_behav_poldegtmp(~isnan(thisb_behav_tmp));  
%         
% % 
% %         [rho_sg(cc,ss,vv,dd), pval_sg(cc,ss,vv,dd)] = corr(thisb_neural', thisb_behav');
% % 
% %         fish_z(cc,ss,vv,dd)  = atanh(rho_sg(cc,ss,vv,dd));
% 
%         [rho_sg_theta(cc,ss,vv,dd), pval_sg_theta(cc,ss,vv,dd)] =corr(thisb_neural', thisb_behav_poldeg');
%         
%         fish_z_theta(cc,ss,vv,dd)  = atanh(rho_sg_theta(cc,ss,vv,dd)); % identical
%         store_fish_ztheta = [store_fish_ztheta; atanh(rho_sg_theta(cc,ss,vv,dd)) cc ss vv dd ];
%         clear mysd
%     end
% 
% %     
%    
%     hold on;
%  
%     plot(vv+0.15,fish_z_theta(cc,:,vv,dd),'o','MarkerSize',5,'Markerfacecolor',[0.3 0.3 0.3], 'Color',[0.3 0.3 0.3]);
%     hold on;
%     my_sem_th = std(fish_z_theta(cc,:,vv,dd))/sqrt(length(subj));
%     plot(vv,mean(fish_z_theta(cc,:,vv,dd)), 'o', 'color',dist_colors(cc,:),'markerfacecolor',dist_colors(cc,:),'markersize',10)
%     plot(vv*[1 1],[mean(fish_z_theta(cc,:,vv,dd))+1.*my_sem_th, mean(fish_z_theta(cc,:,vv,dd))-1.*my_sem_th], '-', 'color',dist_colors(cc,:),'linewidth',2)
%     %title(ROIs{vv})
%     ylim([-1 1])
%     xlim([0.05 12.5])
%     line([0 12.5], [0 0], 'color',[0 0 0],'linewidth',0.5,'linestyle','-') % geh add y =0 baseline 
%     set(gca,'Xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12],'Xticklabel',{'','V1','V2','V3','V3AB','hV4','LO1','IPS0','IPS1','IPS2','IPS3','sPCS',''},'XTickLabelRotation',45,'TickDir','out');
%     ylabel('Fisher Z-transformed Rho')
%  
%     set(gcf,'position',[ 549   724   499   571])
% 
%     
%     [h_th,p_th(cc,vv,dd),~,stats_th] = ttest(fish_z_theta(cc,:,vv,dd));
%     realT_th(cc,vv,dd) = stats_th.tstat;
% 
%      end
% 
% end
% end 

%% figure  5c 

fish_theta_fig =  figure('Name','fish_theta');
store_fish_ztheta= [];
pval_sg_theta = nan(length(cond),length(subj),length(ROIs),length(delay_tpts));
rho_sg_theta = nan(length(cond),length(subj),length(ROIs),length(delay_tpts));
fish_z_theta = nan(length(cond),length(subj),length(ROIs),length(delay_tpts));

figure
for cc = 2 %1:length(cond)
    
for vv=1:length(ROIs)
    
    for dd = 3 %1:length(delay_tpts)
        
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
        
        
        thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp)); 
        thisb_neural = thisb_neural_tmp(~isnan(thisb_behav_tmp));
       
        thisb_behav_poldeg = thisb_behav_poldegtmp(~isnan(thisb_behav_tmp));  
        

        [rho_sg(cc,ss,vv,dd), pval_sg(cc,ss,vv,dd)] = corr(thisb_neural', thisb_behav');

        fish_z(cc,ss,vv,dd)  = atanh(rho_sg(cc,ss,vv,dd));

        [rho_sg_theta(cc,ss,vv,dd), pval_sg_theta(cc,ss,vv,dd)] =corr(thisb_neural', thisb_behav_poldeg');
        
        fish_z_theta(cc,ss,vv,dd)  = atanh(rho_sg_theta(cc,ss,vv,dd)); % identical
        store_fish_ztheta = [store_fish_ztheta; atanh(rho_sg_theta(cc,ss,vv,dd)) cc ss vv dd ];
        clear mysd
    end

%     
   
    hold on;
    plot(vv+0.15,fish_z_theta(cc,:,vv,dd),'o','MarkerSize',5,'Markerfacecolor',[0.3 0.3 0.3], 'Color',[0.3 0.3 0.3]);
    hold on;
    my_sem_th = std(fish_z_theta(cc,:,vv,dd))/sqrt(length(subj));
    plot(vv,mean(fish_z_theta(cc,:,vv,dd)), 'o', 'color',dist_colors(cc,:),'markerfacecolor',dist_colors(cc,:),'markersize',10)
    plot(vv*[1 1],[mean(fish_z_theta(cc,:,vv,dd))+1.*my_sem_th, mean(fish_z_theta(cc,:,vv,dd))-1.*my_sem_th], '-', 'color',dist_colors(cc,:),'linewidth',2)
    %title(ROIs{vv})
    ylim([-1 1])
    xlim([0.05 7.5])
    line([0 8], [0 0], 'color',[0 0 0],'linewidth',0.5,'linestyle','-') % geh add y =0 baseline 

   % if vv==1
    set(gca,'Xtick',[0 1 2 3 4 5 6 7 8 ],'Xticklabel',{'','V1-V3','V3AB','hV4','LO1','IPS0-IPS1','IPS2-IPS3','sPCS',''},'XTickLabelRotation',45,'TickDir','out');
    ylabel('Fisher Z-transformed Rho')
    %else
    set(gcf,'position',[ 549   724   499   571])
    %end 
    
    [h_th,p_th(cc,vv,dd),~,stats_th] = ttest(fish_z_theta(cc,:,vv,dd));
    realT_th(cc,vv,dd) = stats_th.tstat;

     end

end
end 
%%
%%
% % squeeze over subj, collect each subj, all trials into one container
% rho_group_theta = nan(length(cond),length(ROIs),1);
% pval_group_theta = nan(length(cond),length(ROIs),1);
% rho_theta = nan(length(subj),length(ROIs),1);
% pval_theta = nan(length(subj),length(ROIs),1);
% 
% 
% figure 
% 
% for ss = 1:length(subj)
%     hold on;
%     for cc = 2 % 1:length(cond)
%         
%         for vv = 1:length(ROIs)
%             
%             for dd = 3 % 1:length(delay_tpts)
%                 
%                 thisb_neural_tmp =[];
%                 thisb_behav_tmp =[];
%                 thisb_behav = [];
%                 thisb_neural = [];
%                 thisb_behav_poldegtmp =[];
%                 thisb_behav_poldeg =[];
%                 
%                 
%                 thisb_neural_tmp(1,:) =  squeeze(store_b(cc,vv,dd,ss,:));
%                 thisb_behav_tmp(1,:)=  squeeze(mu(2,ss,cc,:))';
%                 thisb_behav_poldegtmp(1,:)=  dtheta_poldeg(ss,cc,:);
%                 
%                 
%                 thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp));
%                 thisb_neural = thisb_neural_tmp(~isnan(thisb_behav_tmp));
%                 thisb_behav_poldeg = thisb_behav_poldegtmp(~isnan(thisb_behav_tmp));
%               
%                 %subplot(length(delay_tpts),length(ROIs),(dd-1)*length(ROIs)+vv); hold on;
%                 subplot(length(subj),length(ROIs),(ss-1)*length(ROIs)+vv); hold on;
%                 %subplot(1,length(ROIs),vv); hold on;
%                 scatter(thisb_neural,thisb_behav_poldeg,30,subj_colors(ss,:),'filled','MarkerFaceAlpha',.4)
%                 
%                 title(ROIs{vv})
%                 lsline
%                 set(gcf,'Position',[  -110         706        1882         624])
%               % xlim([-180 180])
%                %ylim([-25 25])
%                
%                
%                
%                 if vv==1 && ss ==7 
%                     ylabel({'Polar angle \circ delta, target, sacc endpt', 'Epoch 3'})
%                     xlabel('Neural bias, Polar angle \circ')
%                 else
%                 end
%                 [rho_theta(cc,ss,vv,dd), pval_theta(cc,ss,vv,dd)] = corr(thisb_neural', thisb_behav_poldeg');
%                 clear this neural thisb_behav_poldeg
%             end
%             
%             
%         end
%         
%         
%         
%         
%     end
%    
%     
%     
% end
%  set(get(gcf,'Children'),'XLim',[-180 180],'YLim',[-25 25])
%  



%% supplmental figure  5


% fish_theta_fig =  figure('Name','fish_theta');
% store_fish_ztheta= [];
% pval_sg_theta = nan(length(cond),length(subj),length(ROIs),length(delay_tpts));
% rho_sg_theta = nan(length(cond),length(subj),length(ROIs),length(delay_tpts));
% fish_z_theta = nan(length(cond),length(subj),length(ROIs),length(delay_tpts));
% 
% figure
% for cc = 2 %1:length(cond)
%     
% for vv=1:length(ROIs)
%     
%     for dd = 3 %1:length(delay_tpts)
%         
%     for ss = 1:length(subj)
%         
%         thisb_neural_tmp = [];
%         thisb_behav_tmp = [];
%         thisb_neural = [];
%         thisb_behav = [];
%         thisb_behav_poldegtmp =[];
%         thisb_behav_polradtmp =[];
%         thisb_behav_poldeg =[];
%         thisb_behav_polrad =[];
%         
%         thisb_neural_tmp(1,:) =  squeeze(store_b(cc,vv,dd,ss,:));
%         thisb_behav_tmp(1,:)=  squeeze(mu(2,ss,cc,:))';
%         thisb_behav_poldegtmp(1,:)=  dtheta_poldeg(ss,cc,:);
%         
%         
%         thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp)); 
%         thisb_neural = thisb_neural_tmp(~isnan(thisb_behav_tmp));
%        
%         thisb_behav_poldeg = thisb_behav_poldegtmp(~isnan(thisb_behav_tmp));  
%         
% 
%         [rho_sg(cc,ss,vv,dd), pval_sg(cc,ss,vv,dd)] = corr(thisb_neural', thisb_behav');
% 
%         fish_z(cc,ss,vv,dd)  = atanh(rho_sg(cc,ss,vv,dd));
% 
%         [rho_sg_theta(cc,ss,vv,dd), pval_sg_theta(cc,ss,vv,dd)] =corr(thisb_neural', thisb_behav_poldeg');
%         
%         fish_z_theta(cc,ss,vv,dd)  = atanh(rho_sg_theta(cc,ss,vv,dd)); % identical
%         store_fish_ztheta = [store_fish_ztheta; atanh(rho_sg_theta(cc,ss,vv,dd)) cc ss vv dd ];
%         clear mysd
%     end
% 
% %     
%    
%     hold on;
%     plot(vv+0.15,fish_z_theta(cc,:,vv,dd),'o','MarkerSize',5,'Markerfacecolor',[0.3 0.3 0.3], 'Color',[0.3 0.3 0.3]);
%     hold on;
%     my_sem_th = std(fish_z_theta(cc,:,vv,dd))/sqrt(length(subj));
%     plot(vv,mean(fish_z_theta(cc,:,vv,dd)), 'o', 'color',dist_colors(cc,:),'markerfacecolor',dist_colors(cc,:),'markersize',10)
%     plot(vv*[1 1],[mean(fish_z_theta(cc,:,vv,dd))+1.*my_sem_th, mean(fish_z_theta(cc,:,vv,dd))-1.*my_sem_th], '-', 'color',dist_colors(cc,:),'linewidth',2)
%     %title(ROIs{vv})
%     ylim([-1 1])
%     xlim([0.05 7.5])
%     line([0 8], [0 0], 'color',[0 0 0],'linewidth',0.5,'linestyle','-') % geh add y =0 baseline 
% 
%    % if vv==1
%     set(gca,'Xtick',[0 1 2 3 4 5 6 7 8 ],'Xticklabel',{'','V1-V3','V3AB','hV4','LO1','IPS0-IPS1','IPS2-IPS3','sPCS',''},'XTickLabelRotation',45,'TickDir','out');
%     ylabel('Fisher Z-transformed Rho')
%     %else
%     set(gcf,'position',[ 549   724   499   571])
%     %end 
%     
%     [h_th,p_th(cc,vv,dd),~,stats_th] = ttest(fish_z_theta(cc,:,vv,dd));
%     realT_th(cc,vv,dd) = stats_th.tstat;
% 
%      end
% 
% end
% end 
%%
% 
y_var = store_fish_ztheta(:,1); 
cond_var = store_fish_ztheta(:,2); %condition
subj_var = store_fish_ztheta(:,3); %subject
roi_var = store_fish_ztheta(:,4); %roi 
epoch_var = store_fish_ztheta(:,5); %epoch 

% perform real T test

    for vv =1:length(ROIs)
        for cc =1:length(cond)
            for ee = 1:3
                thisidx = roi_var ==vv & cond_var ==cc & epoch_var ==ee;
                y_shuf =y_var(thisidx); %these are already shuffled, see above
                [~,~,~,stats] = ttest(y_shuf);
                real_t(vv,cc,ee) = stats.tstat;
                
                clear thisidx y_shuf
            end   
        end
    end
    
% perform true 1-way anova 


f_real_one = [];
p_real_one = [];


    y = y_var;
    thissubj = subj_var;
    thisroi = roi_var;
    [f_real_one, pval_real_one] = RMAOV1_gh([y,thisroi,thissubj],0.05);
    clear thisroiidx y thiscond thissubj thisepoch 
  

% 
% 
% clear cond_var roi_var epoch_var subj_var y_var
% 
% y_var = store_fish_ztheta(:,1); %fish z 
% cond_var = store_fish_ztheta(:,2); %condition
% subj_var = store_fish_ztheta(:,3); %subject
% roi_var = store_fish_ztheta(:,4); %roi 
% epoch_var = store_fish_ztheta(:,5); %epoch 
% 
% % 1-way 
% 
% f_real_one = [];
% pval_real_one =[];
% 
% for cc =1:length(cond)
% for vv =1:length(ROIs)
% 
%     thisroiidx = roi_var==vv & cond_var==cc;
%     y = y_var(thisroiidx);
%   
%     thissubj = subj_var(thisroiidx);
%     thisepoch = epoch_var(thisroiidx); 
%   
%     [f_real_one(cc,vv), pval_real_one(cc,vv)] = RMAOV1_gh([y,thisepoch,thissubj],0.05);
%     clear thisroiidx y thissubj thisepoch 
%   
% end 
% end 
% 
% 
% clear cond_var roi_var epoch_var subj_var y_var
% check the usage of 'random' per the grouping var 
%%%%%%%%%%%%%%%%%%%%%%%%% PERM TEST # 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we're permuting two measures here, linear bias and ciruclar bias, open
iter = 1000; 
t_iter =cell(length(ROIs),iter);
p_iter =cell(length(ROIs),iter);
f_store_iter_2 = nan(length(ROIs),iter,3);
fish_zthstore_perm =[];
fish_zthperm = nan(length(cond),length(subj),length(ROIs),length(delay_tpts),iter);
shuff_Tth =nan(length(cond),length(ROIs),length(delay_tpts),iter);
rho_permth = nan(length(cond),length(subj),length(ROIs),length(delay_tpts),iter);
pval_permth = nan(length(cond),length(subj),length(ROIs),length(delay_tpts),iter);


for xx = 1:iter
    
    if xx  == 1
        fprintf('starting @ %s \n',datestr(now,'HH:MM:SS'))
        tic
    elseif xx == 2
        toc
    elseif xx == 3
        toc
    elseif xx == iter/4
        fprintf('quarter-way ...\n')
        toc
    elseif xx == iter/2
        fprintf('halfway ...\n')
        toc
    else
    end
    
    for cc = 2 %1:length(cond)
        for vv=1:length(ROIs)
            
            for dd = 3 %1:length(delay_tpts)
                
                for ss = 1:length(subj)
                    thisb_neural_tmp_perm = [];
                    thisb_behav_tmp_perm  = [];
                    thisb_behav_perm  = [];
                    thisb_neural_perm  = [];
                    thisb_behav_poldegtmp = [];
                    thisb_behav_poldegperm = [];
                    thisb_neural_poldegperm = [];
                    shuff_btheta =[];
                    shuff_ntheta =[];
                    
                    
                    thisb_neural_tmp_perm (1,:) =  squeeze(store_b(cc,vv,dd,ss,:));
                    thisb_behav_tmp_perm (1,:)=  squeeze(mu(2,ss,cc,:))';
                    
                    thisb_behav_poldegtmp(1,:)=  dtheta_poldeg(ss,cc,:); %this obv cant change across ROI, delay, lose two dimens (vv,dd) here 

                    
                    thisb_behav_poldegperm = thisb_behav_poldegtmp(~isnan(thisb_behav_tmp_perm));
                    
                    thisb_neural_poldegperm  = thisb_neural_tmp_perm (~isnan(thisb_behav_tmp_perm));
                    
                    % do the shuffle
                    shuff_btheta = thisb_behav_poldegperm; %reload these each time
                    shuff_ntheta = thisb_neural_poldegperm;
                    subjidxth = find(shuff_btheta);
                    shuff_bidxth = randperm(length(subjidxth))';
                    shuff_btheta(subjidxth) = shuff_btheta(shuff_bidxth); %insert shuffled subj indices into bdata (note : in this case, bdata is ONLY this subj. just jumbling idx)
                    
                    bthdata  = shuff_btheta(subjidxth);
                    nthdata  = shuff_ntheta(subjidxth);
                    
                    [rho_permth(cc,ss,vv,dd,xx), pval_permth(cc,ss,vv,dd,xx)] = corr(nthdata',bthdata');
                    fish_zthperm(cc,ss,vv,dd,xx) = atanh(rho_permth(cc,ss,vv,dd,xx));
                    
                    
                    fish_zthstore_perm =[fish_zthstore_perm; fish_zthperm(cc,ss,vv,dd,xx) cc ss vv dd xx]; 
 
                   
                   clear thisidx subjidxth shuff_bidxth  shuff_btheta shuff_ntheta bthdata nthdata
                end
    
                [h,p,~,statsth] = ttest(fish_zthperm(cc,:,vv,dd,xx));
                shuff_Tth(cc,vv,dd,xx) = statsth.tstat;
                
            end
        end
    end
end



% collect original perm pvalues 

p_th_twotail =nan(length(cond),length(ROIs),length(delay_tpts));
p_th_onetail =nan(length(cond),length(ROIs),length(delay_tpts));

pfdr_th_twotail =[];
pmask_th_twotail =[];
pmask_th_onetail =[];
pfdr_th_onetail =[];

p_th_twotaile = [];
p_th_onetaile = [];
pmask_th_twotaile = [];
pmask_th_onetaile = [];
pfdr_th_twotaile = [];
pfdr_th_onetaile = [];


tdistpolbias = figure('Name','tdistpolbias');

realT_col = lines(3);
for cc =1:length(cond)
    for dd = 1:length(delay_tpts)
        for vr = 1:length(ROIs)
            
            hold on;
            subplot(1,length(ROIs),vr)
            histogram(sort(shuff_Tth(cc,vr,dd,:)))
            line([realT_th(cc,vr,dd)  realT_th(cc,vr,dd)], [0 max(ylim)],'LineWidth',0.75,'color',realT_col(dd,:))
            p_th_twotail(cc,vr,dd) = 2 * min ( mean( shuff_Tth(cc,vr,dd,:) <= realT_th(cc,vr,dd) ), mean ( shuff_Tth(cc,vr,dd,:) >= realT_th(cc,vr,dd)));
            p_th_onetail(cc,vr,dd) = mean( shuff_Tth(cc,vr,dd,:) >= realT_th(cc,vr,dd) );
            title(ROIs{vr})
        end
        
          [pfdr_th_twptail(cc,:,dd), pmask_th_twotail(cc,:,dd)] = fdr(p_th_twotail(cc,:,dd),0.05); % this is two tail
          [pfdr_th_onetail(cc,:,dd), pmask_th_onetail(cc,:,dd)] = fdr(p_th_onetail(cc,:,dd),0.05); %one tail
        
        
    end
    
end


%
% two-way anova perm for each ROI?
 % correct across epochs?

 test =  squeeze(p_th_onetail(2,:,3));
 
 [fdrout, maskout] =  fdr(test,0.05)
 

tmp = [squeeze(p_th_onetail(2,:,1)); squeeze(p_th_onetail(2,:,2)); squeeze(p_th_onetail(2,:,3))]

[fdrout, maskout] =  fdr(tmp,0.05);


%% hold off on other perms for now - 102920 




y_var = fish_zthstore_perm(:,1);
cond_var = fish_zthstore_perm(:,2);
subj_var = fish_zthstore_perm(:,3);
roi_var = fish_zthstore_perm(:,4);
epoch_var = fish_zthstore_perm(:,5);
iter_var = fish_zthstore_perm(:,6);


% 
% for zz =1:iter
%     
%     for vv = 1:length(ROIs)
%         
%    
%         thisroiidx = roi_var ==vv & iter_var == zz;
%         y_shuf = y_var(thisroiidx); %these are already shuffled, see above
%         thisepoch = epoch_var(thisroiidx);
%         thiscond =cond_var(thisroiidx);
%         thissubj=subj_var(thisroiidx);
%         
%         [p_iter{vv,zz},t_iter{vv,zz},~] = anovan(y_shuf,{thissubj,thisepoch,thiscond},'model','full','random',1,'varnames',{'subj','epoch','cond'},'display','off');
%         [f_store_iter_2(vv,zz,:)] = RMAOV2_gh([y_shuf,thisepoch,thiscond,thissubj],0.05);
%         
%         clear thisroiidx thisepoch thiscond thissubj y_shuf
% 
%     end
% end
toc

for zz =1:iter
    
    
    thisidx = iter_var == zz;
    y_shuf = y_var(thisidx); %these are already shuffled, see above
    thisroi = roi_var(thisidx);
    thissubj=subj_var(thisidx);
    
    % [p_iter{vv,zz},t_iter{vv,zz},~] = anovan(y_shuf,{thissubj,thisepoch,thiscond},'model','full','random',1,'varnames',{'subj','epoch','cond'},'display','off');
    % [f_store_iter_2(vv,zz,:)] = RMAOV2_gh([y_shuf,thisepoch,thiscond,thissubj],0.05);
    
    [f_store_iter_1(zz,:)] = RMAOV1_gh([y_shuf,thisroi,thissubj],0.05);
    clear thisroiidx thisepoch thiscond thissubj y_shuf
    
end

% visualize 

% iv_str_rm ={'epoch','cond','epoch*cond'};
% which_effect_rm = [1 2 3];
% exact_store_tmp=[];
% exact_store_rm_2 =nan(length(ROIs),length(which_effect_rm));
% 
% for vv =1:length(ROIs)
%     figure('name','2-way perm; RMAOV2')
%     
%     for ww =1:length(which_effect_rm)
%         extract_vals=[];
%         for ii=1:iter
%             extract_vals = [extract_vals; f_store_iter_2(vv,ii,which_effect_rm(ww))];
%         end
%         hold on;
%         
%         subplot(1,length(which_effect_rm),ww)
%         histogram(extract_vals)
%         line([f_real(vv,which_effect_rm(ww)) f_real(vv,which_effect_rm(ww))], [min(ylim) max(ylim)], 'color','r','linewidth',2,'linestyle','--')
%         title(iv_str_rm{which_effect_rm(ww)},'Interpreter','none');
%        % exact_p = sum(extract_vals >= f_real(vv,which_effect_rm(ww)))/iter;
%         exact_p = 2 * min ( sum(extract_vals >= f_real(vv,which_effect_rm(ww)))/iter, sum(extract_vals <= f_real(vv,which_effect_rm(ww)))/iter );
%         
%         text(max(xlim)-(.5*max(xlim)),max(ylim)-(.25*max(ylim)),sprintf('p=%0.4f',exact_p),'FontSize',9)
%         exact_store_rm_2(vv,ww) = [exact_store_tmp; exact_p];
%         clear exact_p
%         sgtitle(ROIs{vv})
%         title(iv_str_rm{ww})
%         if ww ==1
%             xlabel('T-stat')
%             ylabel('Frequency of T-stat')
%         else
%         end
%         
%     end
% end
% 
% % one-way anova perm
% for xx=1:iter
% for cc =1:length(cond)
% for vv =1:length(ROIs)
% 
%     thisroiidx = roi_var==vv & cond_var==cc & iter_var ==xx;
%     y = y_var(thisroiidx);
%     thissubj = subj_var(thisroiidx);
%     thisepoch = epoch_var(thisroiidx); 
%   
%     [f_real_one_perm(cc,vv,xx)] = RMAOV1_gh([y,thisepoch,thissubj],0.05);
%     clear thisroiidx y thissubj thisepoch 
%   
% end 
% end 
% end 
% 
% % visualize 1--way result 
% 
iv_str_rm ={'epoch'};
which_effect_rm = 1;
exact_store_tmp = [];
exact_store_rm_1 =nan(length(ROIs),length(which_effect_rm));

  
for cc =2 
      figure('name','1-way perm; RMAOV1')
    extract_vals= [];
    for ii= 1:iter
        extract_vals = [extract_vals; f_store_iter_1(ii)];
    end

    
    subplot(1,2,cc)
    hold on;
    histogram(extract_vals)
    line([f_real_one f_real_one], [min(ylim) max(ylim)], 'color','r','linewidth',2,'linestyle','--')
    title(iv_str_rm{which_effect_rm(1)},'Interpreter','none');
    exact_p = sum(extract_vals >= f_real_one)/iter;
    
    text(max(xlim)-(.5*max(xlim)),max(ylim)-(.25*max(ylim)),sprintf('p=%0.4f',exact_p),'FontSize',9)
    exact_store_rm_1 = [exact_store_tmp; exact_p];
    clear exact_p
    title(iv_str_rm{1})
    xlabel('T-stat')
    ylabel('Frequency of T-stat')

    
end


[p_fdr_one_out p_fdrthresh_one] = fdr(exact_store_rm_1,0.05)


% do ttest for each roi,epoch,condition

for xx =1:iter
    for vv =1:length(ROIs)
        for cc =1:length(cond)
            for ee = 1:3
                thisidx = roi_var ==vv & cond_var ==cc & epoch_var ==ee & iter_var == xx;
                y_shuf =y_var(thisidx); %these are already shuffled, see above
                [~,~,~,stats] = ttest(y_shuf);
                shuff_t(vv,cc,ee,xx) = stats.tstat;
                
                clear thisidx y_shuf
            end   
        end
    end
end
fprintf('done')

% pvals

for cc =1:length(cond)
    for ee = 1:3
        for vv =1:length(ROIs)
            
            pval_twotail(cc,vv,ee) = 2 * min(mean(shuff_t(vv,cc,ee,:) <= real_t(vv,cc,ee)), mean(shuff_t(vv,cc,ee,:) >= real_t(vv,cc,ee)));
            pval_onetail(cc,vv,ee) = mean(shuff_t(vv,cc,ee,:) >= real_t(vv,cc,ee));
            
        end
          [pval_fdrthresh_twotail(cc,:,ee), pval_th_twotail_mask(cc,:,ee)] = fdr(pval_th_twotail(cc,:,ee),0.05); % this is two tail
          [pval_fdrthresh_onetail(cc,:,ee), pval_th_onetail_mask(cc,:,ee)] = fdr(pval_th_onetail(cc,:,ee),0.05); %one tail
        
    end
end




fprintf('done')
 
fn2s = sprintf('%s/%s_corr/neurbehcorr_concatROIs_fsacc_%iiterperm_%s.mat',root,task_dir,iter,datestr(now,'mmddyyyyHHMM'));
save(fn2s);
fprintf('saving to %s...\n',fn2s);
end

