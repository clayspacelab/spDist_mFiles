% spDist_neural_bias.m

function spDist_neural_bias_perm(subj,sess,ROIs)

root = spDist_loadRoot;


task_dir = 'spDist';

if nargin < 1 || isempty(subj)
    subj = {'CC','KD','AY','MR','XL','SF','EK'};
     

end

if nargin < 2 || isempty(sess)

    sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
  

end

if nargin < 3 || isempty(ROIs)
    ROIs = {'V1','V2','V3','V3AB','hV4','LO1','IPS0','IPS1','IPS2','IPS3','sPCS'};
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

delay_tpt_range = [3.75 5.25; 7.5 9; 10.5 12];



% for fidelity timecourses
tmpcolors = lines(7);

ROI_colors = [repmat(tmpcolors(1,:),3,1); % V1, V2, V3
    tmpcolors(4,:);             % V3AB
    tmpcolors(1,:);             % hV4
    repmat(tmpcolors(3,:),1,1); % VO1
    repmat(tmpcolors(1,:),2,1); % LO1/2
    
    repmat(tmpcolors(2,:),2,1); % TO1-2
    
    repmat(tmpcolors(5,:),2,1); % IPS0-1
    repmat(tmpcolors(6,:),2,1); % IPS2-3
    tmpcolors(7,:);             % sPCS
    ];             % iPCS (color 1...)


%% load data
startidx = 1;
for ss = 1:length(subj)
    
    for vv = 1:length(ROIs)
        
        if cat_mode == 1
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
            startidx = thisidx(end)+1;
            
            clear data;
            
        else
 
            
            for sess_idx = 1:length(sess{ss})
                % build fn
                fn = sprintf('%swmChoose_reconstructions/%s_%s_%s_%s_%ichan%s%s%s_recon_cv_thruTime1.mat',root,subj{ss},sess{ss}{sess_idx},ROIs{vv},func_suffix,nchan,vox_str,smooth_str,trn_str);
                
                fprintf('loading %s...\n',fn);
                data = load(fn);
                
                
                if vv == 1 && ss == 1
                    % initialize variables...
                    
                    
                    nblankt = length(ROIs)*numel(sess)*size(data.recons,1);
                    
                    all_recons = nan(nblankt,size(data.recons,2),size(data.recons,3));
                    all_conds = nan(nblankt,size(data.c_map,2));
                    
                    all_fidelity = nan(nblankt,size(data.recons,3)); % timecoruse of fidelity
                    
                    
                    all_subj = nan(nblankt,1);
                    all_ROIs = nan(nblankt,1);
                    all_sess = nan(nblankt,1);
                    
                    angs = data.angs;
                    tpts = data.delay_tpts;
                    
                    all_r2 = nan(length(ROIs),length(tpts),length(subj));
                    
                end
                
            
                if sess_idx == 1
                    tmp_r2 = nan(length(tpts),length(sess{ss})); % average acorss sessions...
                end
                
                thisidx = startidx:(startidx+size(data.c_map,1)-1);
                
                
                all_recons(thisidx,:,:) = data.recons;
                all_fidelity(thisidx,:) = squeeze(mean(cosd(angs) .* data.recons,2));
                
                
                all_conds(thisidx,:) = data.c_map;
                
                
                
                all_subj(thisidx) = ss;
                
                
                all_ROIs(thisidx) = vv;
                
                all_sess(thisidx) = sess_idx;
                
                tmp_r2(:,sess_idx) = squeeze(mean(mean(data.r2_all,1),2));
                
                startidx = thisidx(end)+1;
                
                clear data;
                
            end
      
            
        end
    end
    
end


%% align like distractor bins (and flip/average) 1D near only  - TRUTH
% goal here is to align cw/ccw distractor bins and flip one set to match
% - for 0-bin, need to determine which trials are CW/CCW and flip
%   accordingly
% new code 
% flip negative angles
% 

roi_str= {'V1','V2','V3','V3AB','hV4','LO1','IPS0','IPS1','IPS2','IPS3','sPCS'};

cond_colors = [0 0 1;.3 .6 .1 ; .3 .6 .1 ; .3 .6 .1];
%tpts_to_plot = (tpts*myTR) >= t_range_to_plot(1) & (tpts*myTR) <= t_range_to_plot(2);
% look for all trials where <> is < 0, flipLR the reconstruction

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
store_b = nan(length(cond_to_plot),length(ROIs),length(subj));
P = nan(length(cond_to_plot),length(ROIs));
H = nan(length(cond_to_plot),length(ROIs));

figure
test_store =[];
for ff = 1:length(cond_to_plot)
    for vv = 1:length(ROIs)
        subplot(length(cond_to_plot),length(ROIs),(ff-1)*length(ROIs)+vv); hold on;
        thisd = nan(size(all_recons{1},3),size(all_recons_flipped,2),length(subj)); %recon data
        thisb = nan(1,length(subj)); %computed bias
        
        for ss = 1:length(subj)
            thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==2 & all_conds(:,6)==0;
            thisd(:,:,ss) = squeeze(mean(all_recons_flipped(thisidx,:,:),1)).';
            thisb(ss) = atan2d(sum(mean(mean(all_recons_flipped(thisidx,:,tpts_to_plot),1),3).*sind(angs)),...
            sum(mean(mean(all_recons_flipped(thisidx,:,tpts_to_plot),1),3).*cosd(angs)));
            store_b(ff,vv,ss) = thisb(ss);
            test_store = [test_store; thisb(ss) vv ss]; 
        end
    
       hold on
       t = mean(thisd((tpts_to_plot),:,:),1);
       myd_sem = nanstd(t,[],3) /sqrt(length(subj));
       plot(linspace(-180,180,90),mean(mean(thisd(tpts_to_plot,:,:),3),1),'LineWidth',2,'color',cond_colors(ff,:))
       hold on;
       plot(linspace(-180,180,90),mean(mean(thisd(tpts_to_plot,:,:),3))+1.*myd_sem,'-','LineWidth',.2,'color',cond_colors(ff,:))
       hold on;
       plot(linspace(-180,180,90),mean(mean(thisd(tpts_to_plot,:,:),3))-1.*myd_sem,'-','LineWidth',.2,'color',cond_colors(ff,:))       
       btwn_d = [mean(mean(thisd(tpts_to_plot,:,:),3))+1.*myd_sem fliplr(mean(mean(thisd(tpts_to_plot,:,:),3))-1.*myd_sem)];     
       fill([linspace(-180,180,90) fliplr(linspace(-180,180,90))],btwn_d,cond_colors(ff,:),'facealpha',0.4);
       mybias_sem = std(store_b(ff,vv,:)) /sqrt(length(subj));
       line([mean(store_b(ff,vv,:),3) mean(store_b(ff,vv,:),3)], [0 2], 'color','k','linewidth',1,'linestyle','-')
       hold on;
       btwn_b=[mean(store_b(ff,vv,:),3)-1*mybias_sem fliplr(mean(store_b(ff,vv,:),3)+1*mybias_sem)];  
       plot([mean(store_b(ff,vv,:),3)-1*mybias_sem mean(store_b(ff,vv,:),3)-1*mybias_sem], [0 2],'r-.','linewidth',0.5)
       plot([mean(store_b(ff,vv,:),3)+1*mybias_sem mean(store_b(ff,vv,:),3)+1*mybias_sem], [0 2],'r-.','linewidth',0.5)
       fill([btwn_b(1) btwn_b(1) btwn_b(2) btwn_b(2)], [0 2 2 0],'r','facealpha',0.2);
       line([0 0], [-.8 max(ylim)], 'color',[.2 .2 .2],'linewidth',0.5,'linestyle','-')
       plot([min(xlim) max(xlim)], [0 0], '-', 'color', [.2 .2 .2])
       x = linspace(-180,180,90);
       line([0 0], [-.8 max(ylim)], 'color',[.2 .2 .2],'linewidth',0.5,'linestyle','-')
       
       if ff == 1 && vv ==1
           title(ROIs{vv});
           
           set(gca,'XTick',-180:90:180,'XTickLabel',{'-180','90','0','90','180'},'XTickLabelRotation',45,'TickDir','out');
           xlabel('Polar angle (\circ)');
           ylabel('Near Distractor');
       else
            title(ROIs{vv});
            set(gca,'XTick',-180:90:180,'XTickLabel',{'','','','',''},'XTickLabelRotation',45,'TickDir','out')
      
       end 
       ylim([-.8 2])
       xlim([-180 180]);
       
       
    end
    clear thisidx
    clear thisb
    clear thisd
end


%store_b dimens
% 1: near distractor
% 2: ROI
% 3: subj
% do rm anova 


[f_out] = RMAOV1_gh([test_store(:,1), test_store(:,2),  test_store(:,3)])
RMAOV1([test_store(:,1), test_store(:,2),  test_store(:,3)])

for tt =1:size(store_b,1) %condition = 1, near;
for rr =1:length(ROIs)
[H(tt,rr), P(tt,rr), CI, STATS] =  ttest(store_b(tt,rr,:));
t_real(rr) =  STATS.tstat;
%fprintf('%s p=%i, t-score =%i, df=%i,sd=%i,mean =', ROIs{rr}, P(tt,rr),STATS.tstat,STATS.df,STATS.sd,mean(store_b(tt,rr,:),3))
clear STATS CI
end 
[p_fdr(tt,:) p_masked(tt,:)] = fdr(P(tt,:),0.05);

end
corrected_pval = P <= p_fdr; 
set(gcf,'Position',[15        1052        2542         286])
%do RM ANOVA to check for differences among offset conditions
%% shuffle on basis of 0 trial, flip random # (length cw trials), compute bias and recon, store
% align like distractor bins (and flip/average) 1D near only 
% goal here is to align cw/ccw distractor bins and flip one set to match
% - for 0-bin, need to determine which trials are CW/CCW and flip
%   accordingly


cond_colors = [0 0 1;.3 .6 .1 ; .3 .6 .1 ; .3 .6 .1];
%tpts_to_plot = (tpts*myTR) >= t_range_to_plot(1) & (tpts*myTR) <= t_range_to_plot(2);
% look for all trials where <> is < 0, flipLR the reconstruction

tmprel =  all_angs(:,2) - all_angs(:,1);
this_rel = mod((tmprel+180), 360)-180;

t_range_to_plot_subset = [10.5 12];
tpts_to_plot = (tpts*myTR) >= t_range_to_plot_subset(1) & (tpts*myTR) < t_range_to_plot_subset(2);
fprintf('using delay tpts = %i to %i', t_range_to_plot_subset(1), t_range_to_plot_subset(2))

cond_to_plot =[0]; %which trial type?
store_b = nan(length(cond_to_plot),length(ROIs),length(subj));
P = nan(length(cond_to_plot),length(ROIs));
H = nan(length(cond_to_plot),length(ROIs));
bshuff_store=[];
iter = 1; 
data_store =nan(iter,length(angs),length(subj)); 
tic
for xx = 1:iter
 
    
for ff = 1:length(cond_to_plot)
    
    for vv = 1:length(ROIs)
        
        thisdata = nan(size(all_recons{1},3),size(all_recons{1},2),length(subj),iter); %recon data
        thisb = nan(1,length(subj),iter); %computed bias
           %all_recons_tmp = all_recons{1}; 
           %all_recons_flipped = all_recons{1}; 
        for ss = 1:length(subj)
               %thisd(:,:,ss) = squeeze(mean(all_recons_flipped(thisidx,:,:),1)).';

      
            thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==2 & all_conds(:,6)==0; %this is collecting all cw/ccw trials  
            findidx = find(thisidx);
            tmp= randperm(length(findidx))';
            shuff_idx = findidx(tmp);    
            flipidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==2 & all_conds(:,6)==0 & this_rel<0; %collect just cw trials, want to flip this many trials 
            nflip = size(find(flipidx),1); % how many trials are truly cw? a #
            tmp_flip = shuff_idx(1:nflip,1); % flip the correct # of trial labels, but ranomly chosen from the shuffled 0 label idx,  
            all_recons_tmp = all_recons{1}; 
            all_recons_flipped = all_recons{1}; 
            flip_log = zeros(length(flipidx),1); 
            flip_log(tmp_flip,1) = 1; % needs to actually be 'logical' not ones, zeros..
            a =logical(flip_log); 
            all_recons_flipped(a,:,:) = fliplr(all_recons_tmp(a,:,:));
            thisdata = all_recons_flipped(shuff_idx,:,tpts_to_plot);  %now access recon that have been randomly flipped with 0 label shuffled  
            data_store(xx,:,ss) = mean(mean(thisdata,1),3);
            thisb(vv,ss,xx) = atan2d(sum(mean(mean(thisdata,1),3).*sind(angs)), sum(mean(mean(thisdata,1),3).*cosd(angs))); %compute bias using shuffling cw/ccw trials 
         
            bshuff_store = [bshuff_store; xx thisb(ss) vv ss]; 
            clear thisdatat shuff_idx tmpd thisidx all_recons_flipped flip_log tmp_log all_recons_tmp tmp_flip n_flip flip_idx
            
        end
           [~, P(vv,xx), ~, stats] =  ttest(thisb(vv,:,xx));
           t_stat_shuff(vv,xx) = stats.tstat; 
   
    end
    
    clear thisb
 

end
end 

toc
fprintf('shuffling done')


set(gcf,'Position',[15        1052        2542         286])
for vv =1:length(ROIs)
condcompare_pval(vv) = mean(t_stat_shuff(vv,:)>=t_real(1,vv));
end
condcompare_fdrthresh = fdr(condcompare_pval(:),0.05)
end

