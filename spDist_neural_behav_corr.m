
function spDist_neural_behav_corr(subj,sess,ROIs)

root = spDist_loadRoot;


task_dir = 'spDist';

if nargin < 1 || isempty(subj)
    subj = {'CC','KD','AY','MR','XL','SF','EK'};
     

end

if nargin < 2 || isempty(sess)

    sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
  

end

if nargin < 3 || isempty(ROIs)
    ROIs = {'V1','V2','V3','V3AB','hV4','VO1','VO2','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'};
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


%% load neural data
startidx = 1;

WHICH_EXCL = [13 20 22]; 
all_data_beh = [];
all_subj_b=[];

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
                this_data.sess_all = sessidx;
                all_data_beh= cat_struct(all_data_beh,this_data);
                %all_subj_b(startidx:(startidx-1+size(this_scored.ii_sess.trialinfo,1))) = ss;
                %all_r_b(startidx:(startidx-1+size(this_scored.ii_sess.trialinfo,1)),1) = (this_data.sess_all*100+this_data.s_all.r_num);
                
            end
            %%%%
         
            
            startidx = thisidx(end)+1;
            
            clear data;
 
    end
    
end
all_data_beh.use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, WHICH_EXCL), all_data_beh.s_all.excl_trial, 'UniformOutput',false));

%% behav import
% subj = {'CC','KD','AY','MR','XL','SF','EK'};
% sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}}; %two sessions removed
% 
% 
% WHICH_EXCL = [13 20 22]; 
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
% all_subj_b = nan(1000*length(subj),1);
% 
% all_data_b = [];
% all_data_bn =[];
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
%         all_data_b = cat_struct(all_data_b,this_data);
%         all_subj_b(startidx:(startidx-1+size(this_scored.ii_sess.trialinfo,1))) = this_subj;
%         all_r_b(startidx:(startidx-1+size(this_scored.ii_sess.trialinfo,1)),1) = (this_data.sess_all*100+this_data.s_all.r_num);
%         
%         %all_data_bn is only used wrt to collecting indices on neural data
%         %keeps track of which trials were thrown out for various reasons,
%         %duplicated per looping subj from ROI
%         all_data_bn(bnidx:(bnidx-1+(size(this_scored.ii_sess.trialinfo,1)*length(ROIs))),1) = repmat( ~cellfun( @any, cellfun( @(a) ismember(a, WHICH_EXCL), this_data.s_all.excl_trial, 'UniformOutput',false) ), length(ROIs),1);
% 
%         bnidx = bnidx + size(this_scored.ii_sess.trialinfo,1)*length(ROIs); 
%         startidx = startidx+size(this_scored.ii_sess.trialinfo,1);
%         clear this_subj this_data;
%     end
% end
% 
% 
% all_subj_b = all_subj_b(1:(startidx-1));
% all_data_b.subj_all = all_subj_b;
% 
% % determine which trials to include
% % first, narrow based on saccade preprocessing/scoring exclusions
% all_data_b.use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, WHICH_EXCL), all_data_b.s_all.excl_trial, 'UniformOutput',false));
% all_data_b.use_trial(all_data_b.s_all.f_sacc_err>8) = 0; %exclude trials with errors > 10 deg
% all_data_b.use_trial(all_data_b.s_all.i_sacc_err>8) = 0;
% %all_data_b.use_trial(ismember(all_data_b.s_all.r_num, [8,12,14]) & all_subj_b==3) = 0; %3 here specifically refers to subj 3 above (AY), see spDist_eyeDataNotes for 
% % drop trials with very short (< 100 ms) or very long RT (> 1 s)
% %all_data_b.use_trial(all_data_b.s_all.i_sacc_rt<0.1 | all_data_b.s_all.i_sacc_rt>1.0) = 0;



%% align like distractor bins (and flip/average) 1D near only 
% goal here is to align cw/ccw distractor bins and flip one set to match
% - for 0-bin, need to determine which trials are CW/CCW and flip
%   accordingly
% new code 
% flip negative angles


tmprel =  all_angs(:,2) - all_angs(:,1);
this_rel = mod((tmprel+180), 360)-180;
% sign of this matches all_conds(:,6) (relative distractor angle bin)

flipidx = this_rel<0;
all_recons_flipped = all_recons{1}; 
all_recons_flipped(flipidx,:) = fliplr(all_recons_flipped(flipidx,:));

t_range_to_plot_subset = [10.5 12];
tpts_to_plot = (tpts*myTR) >= t_range_to_plot_subset(1) & (tpts*myTR) < t_range_to_plot_subset(2);
fprintf('using delay tpts = %i to %i', t_range_to_plot_subset(1), t_range_to_plot_subset(2))

cond_to_plot =[0];
store_b = nan(length(cond_to_plot),length(ROIs),length(subj),36);
P = nan(length(cond_to_plot),length(ROIs));
H = nan(length(cond_to_plot),length(ROIs));


for ff = 1:length(cond_to_plot)
    for vv = 1:length(ROIs)
        thisd = nan(size(all_recons{1},3),size(all_recons_flipped,2),length(subj)); %recon data
        thisb = nan(1,length(subj)); %computed bias
        
        for ss = 1:length(subj)
            
                nruns = 0;
           
                % get one example set of trials
                ru = unique(all_r(all_subj==ss & all_ROIs==vv));
                
                % figure out how many runs we have for that subj, add to total
                nruns = nruns+length(ru);
               for nr =1:nruns

                   thisidx = all_data_beh.use_trial==1 & all_r ==ru(nr) & all_subj==ss & all_ROIs==vv & all_conds(:,1)==2 & all_conds(:,6)==0 ; 
                   thisd(:,:,ss) = squeeze(mean(all_recons_flipped(thisidx,:,:),1)).';
                   thisb(ss) = atan2d(sum(mean(mean(all_recons_flipped(thisidx,:,tpts_to_plot),1),3).*sind(angs)),...
                       sum(mean(mean(all_recons_flipped(thisidx,:,tpts_to_plot),1),3).*cosd(angs)));
                   store_b(ff,vv,ss,nr) = thisb(ss);
               end
         end 

       
    end
    clear thisidx
    clear thisb
    clear thisd
end




%% organize beh data

% % which param should we collect?
% params_of_interest = {'f_sacc'};
% param_str = {'f sacc'};
% 
% mu = nan(2,length(subj),36); % no subj performed > 36 runs 
% err = nan(2,length(subj),36);
% 
% for ss = 1:length(subj)
%     
%                 nruns = 0;
%            
%                 % get one example set of trials
%                 ru = unique(all_r_b(all_subj_b==ss));
%                 
%                 % figure out how many runs we have for that subj, add to total
%                 nruns = nruns+length(ru);
%                 for nr = 1:nruns
%                     %%%%%%% NOT EXCLUDING ON BASIS OF USE_TRIAL %%%%%%%%%%
%                     tmpidx =   all_r_b==ru(nr) & all_subj_b==ss & all_data_b.s_all.trialinfo(:,1)==2  & all_data_b.s_all.trialinfo(:,6)==0 & all_data_b.s_all.trialinfo(:,10) < 0; %negative CW jitter
%                     orig_y =  all_data_b.s_all.(params_of_interest{1})(tmpidx,2);
%                     
%                     if isempty(orig_y)
%                         thisorigidx =  all_r_b==ru(nr) & all_subj_b==ss & all_data_b.s_all.trialinfo(:,1)==2  & all_data_b.s_all.trialinfo(:,6)==0; % skip trying to flip anything if the zero bin trial in this run is not  < 0 jitter 
% %                         if isempty(thisorigidx) %what if the zero bin trial is unusable?
% %                             err(1,ss,nr) = nan;
% %                             mu(:,ss,nr) = nan;
% %                         else
%                             err(1,ss,nr) = nanstd(all_data_b.s_all.(params_of_interest{1})(thisorigidx,2)); 
%                             mu(:,ss,nr) = nanmean( all_data_b.s_all.(params_of_interest{1})(thisorigidx,:),  1 );
%                         %end
%                     else
%                         orig_y_flip = orig_y*-1; % in this run, the zero trial was usuable, and jitter  <0, and we flipped it 
%                         all_data_b.s_all.(params_of_interest{1})(tmpidx,2) = orig_y_flip;  % here is where the actual y flip is inserted for CW bins.
%                         thisorigidx =   all_r_b==ru(nr) & all_subj_b==ss & all_data_b.s_all.trialinfo(:,1)==2  & all_data_b.s_all.trialinfo(:,6)==0; %now, collect all jitters.
%                         err(1,ss,nr) = nanstd(all_data_b.s_all.(params_of_interest{1})(thisorigidx,2)); %lets just take the y-component here. cant STD over 1 number
%                         mu(:,ss,nr) = nanmean( all_data_b.s_all.(params_of_interest{1})(thisorigidx,:),  1 );
%                     end
%       
%                 end
%  
% end 
% 
% 
%     % squeeze over subj, collect each subj, all trials into one container 
%    pval_group=[];
%    rho_group =[];
%     
%        figure ;     
%        for vv=1:length(ROIs)
%            thisb_neural_tmp =[];
%            thisb_behav_tmp =[];
%            for ss = 1:length(subj)
%              
%                thisb_neural_tmp(ss,:) =  squeeze(store_b(1,vv,ss,:));
%                thisb_behav_tmp(ss,:)=  squeeze(mu(2,ss,:))'; 
%            end
%            
%            thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp) & thisb_behav_tmp > -8 &  thisb_behav_tmp < 8); %exclude outliers 
%            thisb_neural = thisb_neural_tmp(~isnan(thisb_behav_tmp) & thisb_behav_tmp > -8 &  thisb_behav_tmp < 8); %use if no exclu 
%     
%            subplot(1,length(ROIs),vv); hold on;
%            scatter(thisb_neural,thisb_behav,30,'filled','MarkerFaceAlpha',.3) 
%            
%            h = lsline;
%            h.LineWidth =2;
%            if vv==1
%                xlabel('Polar angle \circ, Neural bias')
%                ylabel('DVA, saccadic tangential mean error')
%                set(gca,'Xtick',[-180:90:180],'xticklabel',{'180','90','0','90','180'},'xticklabelrotation',45,'TickDir','out')
%                
%            else
%                set(gca,'Xtick',[-180:90:180],'xticklabel',{'','','','',''},'TickDir','out')
%                
%                
%            end
%            title(ROIs{vv})
%            [rho_group(vv), pval_group(vv)] =corr(thisb_neural, thisb_behav);
%            clear thisb_neural thisb_behav h
%            %end
%            
%        end
%            
%             set(gcf,'Position',[15        1052        2542         286])
%             [pfdr(:) pmask(:)]=  fdr(pval_group(:),0.05);
%             
%             all_trials = table(ROIs',rho_group', pval_group',pmask');  
%             all_trials.Properties.VariableNames={'ROIs', 'rho', 'pval','surviveFDR'}       
%            
%      
%            
%            fprintf('done')
%% organize beh data

% which param should we collect?
params_of_interest = {'f_sacc'};
param_str = {'f sacc'};

mu = nan(2,length(subj),36); % no subj performed > 36 runs 
err = nan(2,length(subj),36);

for ss = 1:length(subj)
    
                nruns = 0;
           
                % get one example set of trials
                ru = unique(all_r(all_subj==ss));
                
                % figure out how many runs we have for that subj, add to total
                nruns = nruns+length(ru);
                
                for nr = 1:nruns
                    %%%%%%% NOT EXCLUDING ON BASIS OF USE_TRIAL %%%%%%%%%%
                    tmpidx =  all_data_beh.use_trial==1 & all_ROIs ==1 & all_r==ru(nr) & all_subj==ss & all_data_beh.s_all.trialinfo(:,1)==2  & all_data_beh.s_all.trialinfo(:,6)==0 & all_data_beh.s_all.trialinfo(:,10) < 0; %negative CW jitter
                    %
                    orig_y =  all_data_beh.s_all.(params_of_interest{1})(tmpidx,2);
                    
                    % if isempty(orig_y)
                    %if isempty(thisorigidx) %what if the zero bin trial is unusable?
                    %                             err(1,ss,nr) = nan;
                    %                             mu(:,ss,nr) = nan;
                    %                         else
                    %                             err(1,ss,nr) = nanstd(all_data_beh.s_all.(params_of_interest{1})(thisorigidx,2));
                    %                             mu(:,ss,nr) = nanmean( all_data_beh.s_all.(params_of_interest{1})(thisorigidx,:),  1 );
                    %                         end
                    %                     else
                    orig_y_flip = orig_y*-1; % in this run, the zero trial was usuable, and jitter  <0, and we flipped it
                    all_data_beh.s_all.(params_of_interest{1})(tmpidx,2) = orig_y_flip;  % here is where the actual y flip is inserted for CW bins.
                    thisorigidx = all_data_beh.use_trial==1 &  all_ROIs ==1 & all_r==ru(nr) & all_subj==ss & all_data_beh.s_all.trialinfo(:,1)==2  & all_data_beh.s_all.trialinfo(:,6)==0; %now, collect all jitters.
                    err(1,ss,nr) = nanstd(all_data_beh.s_all.(params_of_interest{1})(thisorigidx,2)); %lets just take the y-component here. cant STD over 1 number
                    mu(:,ss,nr) = nanmean( all_data_beh.s_all.(params_of_interest{1})(thisorigidx,:),  1 );
                end
                
                
 
end 


    % squeeze over subj, collect each subj, all trials into one container 
   pval_group=[];
   rho_group =[];
    
       figure ;     
       for vv=1:length(ROIs)
           thisb_neural_tmp =[];
           thisb_behav_tmp =[];
           for ss = 1:length(subj)
             
               thisb_neural_tmp(ss,:) =  squeeze(store_b(1,vv,ss,:));
               thisb_behav_tmp(ss,:)=  squeeze(mu(2,ss,:))'; 
           end
           
           thisb_behav = thisb_behav_tmp(~isnan(thisb_behav_tmp) & thisb_behav_tmp > -8 &  thisb_behav_tmp < 8); %exclude outliers 
           thisb_neural = thisb_neural_tmp(~isnan(thisb_behav_tmp) & thisb_behav_tmp > -8 &  thisb_behav_tmp < 8); %use if no exclu 
    
           subplot(1,length(ROIs),vv); hold on;
           scatter(thisb_neural,thisb_behav,30,'filled','MarkerFaceAlpha',.3) 
           
           h = lsline;
           h.LineWidth =2;
           if vv==1
               xlabel('Polar angle \circ, Neural bias')
               ylabel('DVA, saccadic tangential mean error')
               set(gca,'Xtick',[-180:90:180],'xticklabel',{'180','90','0','90','180'},'xticklabelrotation',45,'TickDir','out')
               
           else
               set(gca,'Xtick',[-180:90:180],'xticklabel',{'','','','',''},'TickDir','out')
               
               
           end
           title(ROIs{vv})
           [rho_group(vv), pval_group(vv)] =corr(thisb_neural, thisb_behav);
           clear thisb_neural thisb_behav h
           %end
           
       end
           
            set(gcf,'Position',[15        1052        2542         286])
            [pfdr(:) pmask(:)]=  fdr(pval_group(:),0.05);
            
            all_trials = table(ROIs',rho_group', pval_group',pmask');  
            all_trials.Properties.VariableNames={'ROIs', 'rho', 'pval','surviveFDR'}       
           
     
           
           fprintf('done')


end








