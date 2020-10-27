function spDist_neural_behav_corr_sess(subj,sess,ROIs)

root = spDist_loadRoot;


task_dir = 'spDist';

if nargin < 1 || isempty(subj)
    subj = {'CC','KD','AY','MR','XL','SF','EK'};
     

end

if nargin < 2 || isempty(sess)

    sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
  

end

if nargin < 3 || isempty(ROIs)
   % ROIs = {'V1','V2','V3','V3AB','hV4','VO1','VO2','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'};
       ROIs = {'V1','V2','V3','V3AB','hV4','LO1','IPS0','IPS1','IPS2','IPS3','sPCS'}
end


func_suffix = 'surf';

cat_mode = 1; % if 1, look for catSess1Ses...SessN_ files, otherwise, just look for each session in turn

nchan = 8;
which_vox = 0.1; % if > 1 , look for WHICH_VOXvox string; otherwise, look for VE<100*WHICH_VOX>

smooth_by = 1; % if this is 1, use regular files, otherwise, load smooth_by files

myTR = 0.75;

t_range_to_plot = [-inf 12]; % plot b/w these (s)

trn_tpts = 7:15; % if blank, load files w/ no _trn%ito%i, otherwise,

plot_tpts = 7:15; % for a plot where we average reconstructions over a fixed time window

dist_time = 4.5; % onset at 4 s


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
                all_sess_beh(thisbidx,1) = sessidx;
                %all_r_b(startidx:(startidx-1+size(this_scored.ii_sess.trialinfo,1)),1) = (this_data.sess_all*100+this_data.s_all.r_num);
                
                bidx =thisbidx(end)+1;
            end
            %%%%
         
            
            startidx = thisidx(end)+1;
            
            clear data;
 
    end
    
end
all_data_beh.use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, WHICH_EXCL), all_data_beh.s_all.excl_trial, 'UniformOutput',false));
fprintf('pause')
%% behav import
% subj = {'CC','KD','AY','MR','XL','SF','EK'};
% sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}}; %two sessions removed
% 
% 
% WHICH_EXCL = [13 20 21 22]; % see geh spDist_eyeDataNotes.txt on how/why these exclu criteria were chosen. 
% if ismember(WHICH_EXCL,13)
%     which_excl_str ={'broken fix'};
% elseif ismember(WHICH_EXCL,[13 20])
%     which_excl_str ={'broken fix','no sacc'};
% elseif ismember(WHICH_EXCL, [13 20 21])
%     which_excl_str ={'broken fix','no sacc','i sacc too small/short'};
% elseif ismember(WHICH_EXCL, [13 20 22])
%     which_excl_str ={'broken fix','no sacc','i sacc err too lg'};
% elseif ismember(WHICH_EXCL, [13 20 21 22])
%     which_excl_str ={'broken fix','no sacc','i sacc too small/short','i sacc err too lg'};
% else
%     error('which exclusion criteria have you chosen?')
% end
% 
% % first-digit:
% % - 1 - trial-level exclusion (bad drift correction [11], calibration [12], or delay-
% %       fixation break [13]
% % - 2 - primary saccade exclusion (no primary sacc [20]; too small/short [21] bi, large error [22]ei)
% 
% %21 bad initial saccade (duration/amplitude outside range)
% % 22 iniital saccade error
% 
% 
% % concatenate ALL subject data
% all_subj_beh = nan(1000*length(subj),1);
% 
% all_data.beh = [];
% all_data.behn =[];
% startidx = 1;
% bnidx =1;
% for ss = 1:length(subj)
%     for sessidx = 1:length(sess{ss})
%         
%         
%         fn = sprintf('%s/spDist_behav_92220/%s_%s_scored.mat',root,subj{ss},sess{ss}{sessidx});
%         fprintf('Loading scored eye data from %s\n',fn);
%         this_scored = load(fn);
%         
%         this_data.s_all = this_scored.ii_sess;
%         this_data.sess_all = sessidx;
%         
%         this_subj = ss;
%         
%         all_data.beh = cat_struct(all_data.beh,this_data);
%         all_subj_beh(startidx:(startidx-1+size(this_scored.ii_sess.trialinfo,1))) = this_subj;
%         all_r_b(startidx:(startidx-1+size(this_scored.ii_sess.trialinfo,1)),1) = (this_data.sess_all*100+this_data.s_all.r_num);
%         
%         all_data.behn(bnidx:(bnidx-1+(size(this_scored.ii_sess.trialinfo,1)*length(ROIs))),1) = repmat( ~cellfun( @any, cellfun( @(a) ismember(a, WHICH_EXCL), this_data.s_all.excl_trial, 'UniformOutput',false) ), length(ROIs),1);
%         all_sess_b(startidx:(startidx-1+size(this_scored.ii_sess.trialinfo,1)),1) = this_data.sess_all; 
%         
%         bnidx = bnidx + size(this_scored.ii_sess.trialinfo,1)*length(ROIs); 
%         startidx = startidx+size(this_scored.ii_sess.trialinfo,1);
%         clear this_subj this_data;
%     end
% end
% 
% 
% all_subj_beh = all_subj_beh(1:(startidx-1));
% all_data.beh.subj_all = all_subj_beh;
% 
% % determine which trials to include
% % first, narrow based on saccade preprocessing/scoring exclusions
% all_data.beh.use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, WHICH_EXCL), all_data.beh.s_all.excl_trial, 'UniformOutput',false));
% all_data.beh.use_trial(all_data.beh.s_all.f_sacc_err>10) = 0; %exclude trials with errors > 10 deg
% all_data.beh.use_trial(all_data.beh.s_all.i_sacc_err>10) = 0;
% %all_data.beh.use_trial(ismember(all_data.beh.s_all.r_num, [8,12,14]) & all_subj_beh==3) = 0; %3 here specifically refers to subj 3 above (AY), see spDist_eyeDataNotes for 
% % drop trials with very short (< 100 ms) or very long RT (> 1 s)
% %all_data.beh.use_trial(all_data.beh.s_all.i_sacc_rt<0.1 | all_data.beh.s_all.i_sacc_rt>1.0) = 0;
% 


%% align like distractor bins (and flip/average) 1D near only 
% goal here is to align cw/ccw distractor bins and flip one set to match
% - for 0-bin, need to determine which trials are CW/CCW and flip
%   accordingly
% new code 
% flip negative angles
% 

%roi_str= {'V1','V2','V3','V3AB','hV4','LO1','IPS0','IPS1','IPS2','IPS3','sPCS'};

tmprel =  all_angs(:,2) - all_angs(:,1);
this_rel = mod((tmprel+180), 360)-180;
% sign of this matches all_conds(:,6) (relative distractor angle bin)

flipidx = this_rel<0;
all_recons_flipped = all_recons{1}; 
all_recons_flipped(flipidx,:) = fliplr(all_recons_flipped(flipidx,:));

t_range_to_plot_subset = [10.5 12];
tpts_to_plot = (tpts*myTR) >= t_range_to_plot_subset(1) & (tpts*myTR) < t_range_to_plot_subset(2);
sprintf('using delay tpts = %i to %i', t_range_to_plot_subset(1), t_range_to_plot_subset(2))

cond_to_plot =[0];
store_b = nan(length(cond_to_plot),length(ROIs),length(subj),2);
P = nan(length(cond_to_plot),length(ROIs));
H = nan(length(cond_to_plot),length(ROIs));
sessions =[1 2];

figure
test_store =[];
for ff = 1:length(cond_to_plot)
    for vv = 1:length(ROIs)
        subplot(length(cond_to_plot),length(ROIs),(ff-1)*length(ROIs)+vv); hold on;
        thisd = nan(size(all_recons{1},3),size(all_recons_flipped,2),length(subj)); %recon data
        thisb = nan(1,length(subj)); %computed bias
        
        for ss = 1:length(subj)
            
                
               for sh =1:length(sessions)
                   thisidx = all_data_beh.use_trial==1 & all_sess ==sessions(sh) & all_subj==ss & all_ROIs==vv & all_conds(:,1)==2 & all_conds(:,6)==0 ; % %this needs to be further checked..  all_data.behn==1;

                  % thisidx = all_sess ==sessions(sh) & all_subj==ss & all_ROIs==vv & all_conds(:,1)==2 & all_conds(:,6)==0 ; % %this needs to be further checked..  all_data.behn==1;
                   thisd(:,:,ss) = squeeze(mean(all_recons_flipped(thisidx,:,:),1)).';
                   thisb(ss) = atan2d(sum(mean(mean(all_recons_flipped(thisidx,:,tpts_to_plot),1),3).*sind(angs)),...
                       sum(mean(mean(all_recons_flipped(thisidx,:,tpts_to_plot),1),3).*cosd(angs)));
                   store_b(ff,vv,ss,sh) = thisb(ss);
                   test_store = [test_store; thisb(ss) vv ss sh];
               end
         end 

       
    end
    clear thisidx
    clear thisb
    clear thisd
end

% y_tmp = test_store(:,1);
% roi_var_tmp = test_store(:,2);
% subj_var_tmp = test_store(:,3);
% 
% y = y_tmp(~isnan(y_tmp));
% roi_var = roi_var_tmp(~isnan(y_tmp));
% subj_var = subj_var_tmp(~isnan(y_tmp));


[f_out] = RMAOV1_gh([test_store(:,1), test_store(:,2),  test_store(:,3)]);
RMAOV1([test_store(:,1), test_store(:,2),  test_store(:,3)])
%RMAOV1([y, roi_var, subj_var])

for tt =1:size(store_b,1) %condition = 1, near;
for rr =1:length(ROIs)
[H(tt,rr), P(tt,rr), CI, STATS] =  ttest(store_b(tt,rr,:));
t_real =  STATS.tstat;
clear STATS CI
end 
[p_fdr(tt,:) p_masked(tt,:)] = fdr(P(tt,:),0.05);

end
corrected_pval = P <= p_fdr; 
set(gcf,'Position',[15        1052        2542         286])

%% organize beh data
distractor_bins = unique(all_data.beh.s_all.trialinfo(all_data.beh.s_all.trialinfo(:,1)~=1,6));
distractor_spacing = 360/length(distractor_bins);

flip_bins = [1 2 3]; %which distractor bins should we "flip" the Y in order to collapse across CW/CCW distractors?


cond_str = {'No distractor','Distractor'};


params_of_interest = {'f_sacc'};
param_str = {'f sacc'};


subj_col = lines(7); 
mu = nan(2,length(subj)); % no subj performed > 36 runs 
err = nan(2,length(subj));

for ss = 1:length(subj)
    
     
                for sh = 1:length(sessions)                    
                    %%%%%%% NOT EXCLUDING ON BASIS OF USE_TRIAL %%%%%%%%%%
                   % tmpidx =  all_sess_beh==sessions(sh) & all_subj_beh==ss & all_data_beh.s_all.trialinfo(:,1)==2  & all_data_beh.s_all.trialinfo(:,6)==0 & all_data_beh.s_all.trialinfo(:,10) < 0; %negative CW jitter
% 
                    tmpidx = all_data_beh.use_trial==1 &  all_sess_beh==sessions(sh) & all_subj_beh==ss & all_data_beh.s_all.trialinfo(:,1)==2  & all_data_beh.s_all.trialinfo(:,6)==0 & all_data_beh.s_all.trialinfo(:,10) < 0; %negative CW jitter

                    orig_y =  all_data_beh.s_all.(params_of_interest{1})(tmpidx,2);
                    orig_y_flip = orig_y*-1;
                    all_data_beh.s_all.(params_of_interest{1})(tmpidx,2) = orig_y_flip;  % here is where the actual y flip is inserted for CW bins.
                    
                    thisorigidx = all_data_beh.use_trial==1 & all_sess_beh== sessions(sh) & all_subj_beh==ss & all_data_beh.s_all.trialinfo(:,1)==2  & all_data_beh.s_all.trialinfo(:,6)==0; %now, collect all jitters.
                    
                    % distractor bin x param x [radial; tangential] x subj
                    err(:,ss,sh) = nanstd(all_data_beh.s_all.(params_of_interest{1})(thisorigidx,:), [], 1 ); % this wont work here bc we're working with
                    mu(:,ss,sh) = nanmean( all_data_beh.s_all.(params_of_interest{1})(thisorigidx,:),  1 );
                    
                end
          
    
end 

%plot just the behav bias 
figure
subplot(1,length(params_of_interest),1); hold on;
for ss = 1:length(subj)
    thism_near = squeeze(mu(2,ss,:))'; %this is the only condition that is information for this analysis.
    plot(ss,thism_near,'o','MarkerSize',5, 'Color',subj_col(ss,:),'markerfacecolor',subj_col(ss,:));
    hold on;
    plot(ss,nanmean(thism_near),'o','MarkerSize',8, 'Color',subj_col(ss,:),'markerfacecolor',subj_col(ss,:));
    set(gca,'xlim',[0.5 length(subj)+0.5], 'xtick',[1:length(subj)],'xticklabel', [1:length(subj)])
    clear thism_near
end
title('mu')
ylabel('DVA, tangential mean error')
xlabel('subj')
           
figure
for vv=1:length(ROIs)
    subplot(1,length(ROIs),vv); hold on;
    for ss = 1:length(subj)
        thism_near = squeeze(store_b(1,vv,ss,:))'; %this is the only condition that is information for this analysis.
        plot(ss,thism_near,'o','MarkerSize',3, 'Color',subj_col(ss,:),'markerfacecolor',subj_col(ss,:));
        hold on;
        plot(ss,nanmean(thism_near),'o','MarkerSize',8, 'Color',subj_col(ss,:),'markerfacecolor',subj_col(ss,:));
        set(gca,'xlim',[-0.5 length(subj)+0.5], 'xtick',[1:length(subj)],'xticklabel', [1:length(subj)])
        
    end
    ylim([-180 180])
    title(ROIs{vv})
    if vv==1
        ylabel('polar angle')
        xlabel('subj')
    else
    end
    
end
  %group avg, simple
 
  figure
  thism_near_tmp = squeeze(mu(2,:,:)); %this is the only condition that is information for this analysis.
  thism_near = [thism_near_tmp(:,1); thism_near_tmp(:,2)];
  
  plot(1,thism_near,'o','MarkerSize',3, 'Color',subj_col(ss,:),'markerfacecolor',[.5 .5 .5]);
  hold on;
  plot(1,nanmean(thism_near),'o','MarkerSize',8, 'Color',subj_col(ss,:),'markerfacecolor',[.5 .5 .5]);
  tmpe = std(thism_near)/sqrt(length(subj));
  plot(1*[1 1],nanmean(thism_near)+tmpe*[-1 1],'-','LineWidth',1.5,'Color','k');
  
  if vv==1
      ylabel('DVA, mu')
      xlabel('subj')
  else
  end
  
  [h,p]= ttest(thism_near)
  
  title('Bias')
% $ fig 
  figure ;
  for vv=1:length(ROIs)
      
      subplot(1,length(ROIs),vv); hold on;
      thisb_neural_tmp = squeeze(store_b(1,vv,:,:)); %this is the only condition that is information for this analysis.
      thisb_neural = [thisb_neural_tmp(:,1); thisb_neural_tmp(:,2)];
      thisb_behav_tmp = squeeze(mu(2,:,:)); %this is the only condition that is information for this analysis.
      thisb_behav =[thisb_behav_tmp(:,1); thisb_behav_tmp(:,2)];
      
      %thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp));
      %thisb_neural = thisb_neural(~isnan(thisb_behav_tmp)); %need to trim based on behav, nan's are unavoidable
      plot(thisb_neural,thisb_behav,'o','markerface',subj_col(ss,:),'markerfacecolor',subj_col(ss,:))
      lsline
      if vv==1 && ss ==7
          xlabel('Polar angle \circ, Neural bias')
          ylabel('DVA, saccadic tangential mean error (Bias)')
          set(gca,'Xtick',[-180:90:180],'Xticklabel',{'180','90','0','90','180'},'xticklabelrotation',45,'TickDir','out')
      else
          set(gca,'xtick',[-180:90:180],'xticklabel',{'','','','',''},'TickDir','out')
      end
      title(ROIs{vv})
      [rho(vv), pval(vv)] =corr(thisb_neural, thisb_behav)
      clear thisb_neural thisb_behav thisb_behav_tmp thisb_neural_tmp
  end
  
  %do corr using ERR
           figure ;     
           for vv=1:length(ROIs)
        
                   subplot(1,length(ROIs),vv); hold on;
                   thisb_neural_tmp = squeeze(store_b(1,vv,:,:)); %this is the only condition that is information for this analysis.    
                   thisb_behav_tmp = squeeze(err(1,:,:)); %this is the only condition that is information for this analysis.
                   thisb_neural = [thisb_neural_tmp(:,1); thisb_neural_tmp(:,2)];
                   thisb_behav = [thisb_behav_tmp(:,1); thisb_behav_tmp(:,2)];
                   
                   
                   %thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp));
                   %thisb_neural = thisb_neural(~isnan(thisb_behav_tmp)); %need to trim based on behav, nan's are unavoidable

                   plot(thisb_behav,thisb_neural,'o','markerface',subj_col(ss,:),'markerfacecolor',subj_col(ss,:))
                   lsline
                   if vv==1 && ss ==7
                   xlabel('Polar angle \circ, Neural bias')
                   ylabel('DVA, saccadic tangential STD error')
                   set(gca,'ytick',[-180:90:180],'yticklabel',{'180','90','0','90','180'},'yticklabelrotation',45,'TickDir','out')
                   set(gca,'xtick',[-2:0.5:2],'yticklabelrotation',45,'TickDir','out')

                   else
                   set(gca,'ytick',[-180:90:180],'xticklabel',{'','','','',''},'TickDir','out')
                   set(gca,'xtick',[-2:0.5:2],'yticklabelrotation',45,'TickDir','out')

                       
                   end
                   title(ROIs{vv})
                   [rho(vv), pval(vv)] =corr(thisb_neural, thisb_behav)
                   clear thisb_neural thisb_behav
           end
           
%group avg  over runs        
           figure ;     
           for vv=1:length(ROIs)
        
              for ss = 1:length(subj)
                   subplot(1,length(ROIs),vv); hold on;
                   thisb_neural = squeeze(store_b(1,vv,ss,:)); %this is the only condition that is information for this analysis.    
                   thisb_behav= squeeze(mu(2,ss,:)); %this is the only condition that is information for this analysis.
                   thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp));
                   thisb_neural = thisb_neural(~isnan(thisb_behav_tmp)); %need to trim based on behav, nan's are unavoidable

                   plot(thisb_neural,thisb_behav,'o','markerface',subj_col(ss,:),'markerfacecolor',subj_col(ss,:))
                   lsline
                   if vv==1 
                   xlabel('Polar angle \circ, Neural bias')
                   ylabel('DVA, saccadic tangential mean error')
                   set(gca,'Xtick',[-180:90:180],'xticklabel',{'180','90','0','90','180'},'xticklabelrotation',45,'TickDir','out')

                   else
                   set(gca,'Xtick',[-180:90:180],'xticklabel',{'','','','',''},'TickDir','out')
  
                       
                   end
                   title(ROIs{vv})
                   [rho_group(vv), pval_group(vv)] =corr(thisb_neural, thisb_behav)
                   clear thisb_neural thisb_behav
               end
           end
           
            figure ;     
           for vv=1:length(ROIs)
        
              %for ss = 1:length(subj)
                   subplot(1,length(ROIs),vv); hold on;
                   thisb_neural = squeeze(nanmean(store_b(1,vv,:,:),4)); %this is the only condition that is information for this analysis.    
                   thisb_behav= squeeze(nanmean(mu(2,:,:),3))'; %this is the only condition that is information for this analysis.
                  % thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp))';
                   %thisb_neural = thisb_neural(~isnan(thisb_behav_tmp)); %need to trim based on behav, nan's are unavoidable

                   plot(thisb_neural,thisb_behav,'o','markerface',subj_col(ss,:),'markerfacecolor',subj_col(ss,:))
                   lsline
                   if vv==1 
                   xlabel('Polar angle \circ, Neural bias')
                   ylabel('DVA, saccadic tangential mean error')
                   set(gca,'Xtick',[-180:90:180],'xticklabel',{'180','90','0','90','180'},'xticklabelrotation',45,'TickDir','out')

                   else
                   set(gca,'Xtick',[-180:90:180],'xticklabel',{'','','','',''},'TickDir','out')
  
                       
                   end
                   title(ROIs{vv})
                   [rho_group(vv), pval_group(vv)] =corr(thisb_neural, thisb_behav)
                   clear thisb_neural thisb_behav
               %end
           end
           
           figure; imagesc(rho)
           ylabel('ROI')
           set(gca,'Ytick',[1:1:length(ROIs)],'Yticklabel',{'V1','V2','V3','V3AB','hV4','LO1','IPS0','IPS1','IPS2','IPS3','sPCS'})
           title('Trial-by-trial correlations, f sacc bias & neural bias')
           caxis([-1 1])
           colorbar
           fprintf('done')


end


