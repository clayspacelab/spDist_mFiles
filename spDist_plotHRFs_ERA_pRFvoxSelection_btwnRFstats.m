% spDist_pilot_scanner_plotHRFs_ERA.m
%
% loads trialData files & plots HRFs for distractor/no-distractor trials
% (Note: different # of trials...)
%
% like fidelity, compares matched trial epochs
%
% stats: 3-way shuffled ANOVA, then 2-way shuffled ANOVAs per ROI
%
% TODO: use 'plot' instead of 'text' for significance markers
%

function spDist_plotHRFs_ERA_pRFvoxSelection_btwnRFstats(subj,sess,ROIs,alignment)


task_dir = 'spDist';

root = spDist_loadRoot;

rng(spDist_randSeed);

if nargin < 1 || isempty(subj)
    subj = {'AY','CC','EK','KD','MR','SF','XL'};
end

if nargin < 2 || isempty(sess)
    sess_template = {'spDist1','spDist2'};
    sess = cell(length(subj),1); for ss = 1:length(subj); sess{ss} = sess_template; end 
    clear sess_template
end


if nargin < 3 || isempty(ROIs)
    % all ROIs
    % ROIs = {'V1','V2','V3','V3AB','hV4','LO1','IPS0','IPS1','IPS2','IPS3','sPCS'};
    %  ROIs = {'V1','V2','V3','V3AB','IPS0','IPS1','IPS2','IPS3','sPCS'};
    ROIs = {'V1V2V3','V3AB','hV4','LO1','IPS0IPS1','IPS2IPS3','sPCS'};
    %ROIs = {'V1V2V3','V3AB','IPS0IPS1','IPS2IPS3','sPCS'};
end


if nargin < 4 || isempty(alignment)
   alignment= 'targ_ang_all';
end

if alignment == 'dist_ang_all';
    do_distalign_plot =1;
else
    do_distalign_plot =0;
end

% do stats? (takes ~10 mins)
do_stats = 0;

% inspect # voxels?
make_vox_plot = 0;

% if so, save them?
save_stats = 0; 



func_suffix = 'surf';

ve_thresh = 0.1; 

pol_diff_thresh_in = 15; 
pol_diff_thresh_out = 165; 
min_ecc = 2; 
max_ecc = 15; 


vox_inn =nan(length(subj),2,length(ROIs),360);
vox_outn =nan(length(subj),2,length(ROIs),360);
% 'overall' delay-period (maximizes difference between distractor
% present/absent
% (because we use less than/equal to for tpt selection...)
% (this is based on tr - we want to INCLUDE TRs that span this range)
delay_range = [9 15]; % TR beginning at 7, ending at 15, to match IEM training
%delay_range = [7 15]; 
% these are based on *time*, not TR...
% as in Fig. 1, 4, 5 - 3 epochs for comparing before/during/after distractor
delay_tpt_range = [3.75 5.25; 8.25 9.75; 10.5 12]; %for stats - up for discussion?
% >= n,1 and < n,2

epoch_str = {'PRE','DIST','POST'};

t_markers = [0 4.5 12]; % beginning of delay, beginning of distractor, beginning of response


%% load data
startidx = 1;

for ss = 1:length(subj)
    
    for sess_idx = 1:length(sess{ss})
        
        for vv = 1:length(ROIs)
            
            fn = sprintf('%s/%s_trialData/%s_%s_%s_%s_trialData.mat',root,task_dir,subj{ss},sess{ss}{sess_idx},ROIs{vv},func_suffix);
            fprintf('loading %s...\n',fn);
            data = load(fn);
            
            
            if vv == 1 && ss == 1 && sess_idx == 1
                % initialize variables...
                
                nblank = length(ROIs)*numel(sess)*size(data.dt_all,1);
                
                all_hrfs_in = nan(nblank,size(data.dt_all,3));
                all_hrfs_out = nan(nblank,size(data.dt_all,3));

                all_conds = nan(nblank,size(data.c_all,2));
                
                all_subj = nan(nblank,1);
                all_ROIs = nan(nblank,1);
                
                all_sess = nan(nblank,1);

                
                TR = data.TR;
                which_TRs = data.which_TRs;
                
                
            end
            
            
            thisidx = startidx:(startidx+size(data.dt_all,1)-1);
            
            %%% additions for selecting based on ecc, pol ang diff %%%
            tmp_x = data.rf.x0; %this is in unit?
            tmp_y =data.rf.y0.*-1; % because vista gives -y as upper vf 
            tmp_pol = atan2d(tmp_y,tmp_x); 
            tmp_ecc = data.rf.ecc;
            
            for pp = 1:length(data.targ_ang_all)
                clear align_loc
                align_loc = data.(sprintf('%s',alignment));
                tmprel =  align_loc(pp) - tmp_pol; % align_loc will either be targ or dist 
                pol_diff(1,:)= abs(mod((tmprel+180), 360)-180); %dont care about sign here, just need angular difference magnitude
                
                %vox selection criteria using pol ang sep, max, min ecc 
                which_vox_in = data.rf.ve>=ve_thresh & pol_diff(1,:) <= pol_diff_thresh_in & tmp_ecc > min_ecc & tmp_ecc < max_ecc;
                which_vox_out = data.rf.ve>=ve_thresh & pol_diff(1,:) >= pol_diff_thresh_out & tmp_ecc > min_ecc & tmp_ecc < max_ecc;
                 
                vox_inn(ss,sess_idx,vv,pp) = sum(which_vox_in);
                vox_outn(ss,sess_idx,vv,pp) = sum(which_vox_out);
                
                idx = thisidx(pp);
                all_hrfs_in(idx,:) = squeeze(nanmean(data.dt_allz(pp,which_vox_in,:),2))';
                all_hrfs_out(idx,:) = squeeze(nanmean(data.dt_allz(pp,which_vox_out,:),2))';
                clear pol_diff tmp_rel idx
                
            end
            
            all_conds(thisidx,:) = data.c_all;
            
            all_subj(thisidx) = ss;
            
            all_ROIs(thisidx) = vv;
            
            all_sess(thisidx) = sess_idx;
            
            startidx = thisidx(end)+1;
            
            clear data tmp_y tmp_x tmp_pol tmp_ecc
            
        end
    end
    
end

%% remove baseline
baseline_TRs = which_TRs < 0;
all_hrfs_in = all_hrfs_in - mean(all_hrfs_in(:,baseline_TRs),2);
all_hrfs_out = all_hrfs_out - mean(all_hrfs_out(:,baseline_TRs),2);

%% plot mean HRF activations throughout the delay-period. Sort by RF condition and distractor absent (mag) v distractor present (blue)

%TOP: RFin 
%MIDDLE: RFout 
%BOTTOM: RFin-out

all_HRFs{1} = all_hrfs_in;
all_HRFs{2} = all_hrfs_out;
%all_HRFs{3} = all_hrfs_in - all_hrfs_out;
all_HRFs{3} = ((all_hrfs_in + 100) - (all_hrfs_out +100)) ./ ((all_hrfs_in + 100) + (all_hrfs_out + 100));
cond_colors = spDist_condColors; 

cu = unique(all_conds(:,1));
all_mean_hrf = cell(length(all_HRFs),length(cu),1); 


axhrf = nan(length(all_HRFs),length(ROIs));
mh = nan(length(ROIs),length(t_markers));


% start with just subplots of timeseries (w/ errorbars?)
figure;
for hh = 1:length(all_HRFs)
for vv = 1:length(ROIs)
    
    axhrf(hh,vv) = subplot(length(all_HRFs),length(ROIs),(hh-1)*length(ROIs)+vv); hold on;
    % draw 'baseline'
    plot([which_TRs(1) which_TRs(end)]*TR,[0 0],'k-','LineWidth',0.75);
    
    % draw event markers
    mh(vv,:) = plot(t_markers.*[1;1],[0 .1],'-','Color',[0.7 0.7 0.7],'LineWidth',0.75);
    
    for cc = 1:length(cu)
        
        thisd = nan(length(subj),size(all_HRFs{1},2));
        for ss = 1:length(subj)
            
            thisidx = all_subj == ss & all_ROIs == vv & all_conds(:,1)==cu(cc);
            
            thisd(ss,:) = nanmean(all_HRFs{hh}(thisidx,:),1); %using nanmean here 
            all_mean_hrf{hh,cc}(vv,:,ss) = nanmean(all_HRFs{hh}(thisidx,:),1);
            clear thisidx;
        end
        
   
        % plot, like for reconstructions, such that middle of each
        % datapoint is at middle of TR
        if hh ==1
            plot(which_TRs*TR + TR/2,nanmean(thisd,1),'-','LineWidth',1.0,'Color',cond_colors(cc,:));

        elseif hh ==2
            plot(which_TRs*TR + TR/2,nanmean(thisd,1),'--','LineWidth',1.0,'Color',cond_colors(cc,:)); % dashed for RFout
        elseif hh ==3
            plot(which_TRs*TR + TR/2,nanmean(thisd,1),':','LineWidth',1.0,'Color',cond_colors(cc,:));

        end

        btwn_fill = [nanmean(thisd,1)+1.*nanstd(thisd,[],1)/sqrt(length(subj)) fliplr( nanmean(thisd,1)-1.*nanstd(thisd,[],1)/sqrt(length(subj)) )]; % geh test
        
        fill([which_TRs*TR + TR/2 fliplr(which_TRs*TR + TR/2)],btwn_fill,cond_colors(cc,:),'linestyle','none','facealpha',0.3);
 
        clear thisd;
    end
    
    title(ROIs{vv});
    if vv==1 && hh ==1
        ylabel({'Voxels in RF'; 'BOLD Z-score'});
        xticks([0 4.5 12])
        set(gca,'XTickLabel',[0 4.5 12])
        xlabel('Time (s)');
        yticks([-.2:.2:1.4])
    elseif vv==1 && hh ==2
        ylabel({'Voxels out of RF ';'BOLD Z-score'});
        set(gca,'xTickLabel',[]);
    elseif vv==1 && hh ==3
        ylabel({'Voxels in RF -out of RF';'BOLD Z-score'});
        xticks([0 4.5 12])
        set(gca,'XTickLabel',[0 4.5 12])
        xlabel('Time (s)');
    else
        yticks([-.2:.2:1.4])
        xticks([0 4.5 12])
        set(gca,'YTickLabel',[]);
        set(gca,'xTickLabel',[]);
    end
    
    
    
end
myy = cell2mat(get(axhrf(hh,:),'YLim'));
set(axhrf(hh,:),'YLim',[min(myy(:,1)) max(myy(:,2))],'XLim',[which_TRs(1) which_TRs(end)]*TR,'TickDir','out');
set(mh,'YData',[min(myy(:,1)) max(myy(:,2))]);
if length(ROIs)==7
set(gcf,'Position',[ 157         697        2041         594])
else
set(gcf,'Position',[ 157         697        2541         594])
end
%match_ylim(get(gcf,'Children'));
end


%% plot # voxels per trial, subject, ROI, RF-in & RF-out
if make_vox_plot ==1
    
subj_col =lines(7);
myd = {};
myd{1} = vox_inn;
myd{2} = vox_outn;

dplot = nan(360,1);
d_store=[];

figure(4)
for ss=1:length(subj)
    
    for vv =1:length(ROIs)
        
        
        dtmpi = squeeze(myd{1}(ss,:,vv,:))';
        di = dtmpi(~isnan(dtmpi));
        dtmpo = squeeze(myd{2}(ss,:,vv,:))';
        do = dtmpo(~isnan(dtmpo));
        
        figure(4)
        subplot(length(subj),length(ROIs),(ss-1)*length(ROIs)+vv);
        hold on;
        imagesc([di(:) do(:)])
        xlim([0 3])
        clear dtmpi dtmpo
        
        
        
        if ss==1
            
            title(ROIs{vv})
        elseif ss==1 && vv==1
            ylabel('# voxels per trial')
            xlabel(sprintf('%s',subj{ss}))
        else
        end
        
    end
end


end
%% plot within dist condition & between RF condition
HRF_str ={'RFin','RFout','RFin-out'};
%HRF_str ={'RFin','RFout','RFin-out','RF mod'};

% store something that's ROI x time x subj for each condition
cu = unique(all_conds(:,1));

axhrf = nan(length(all_HRFs),length(ROIs));
mh = nan(length(ROIs),length(t_markers));

condstr = {'Dist Absent','Dist Present'}; 

 figure;
for cc = 1:length(cu)
   
for vv = 1:length(ROIs)
    axhrf(cc,vv) = subplot(length(cu),length(ROIs),(cc-1)*length(ROIs)+vv); hold on;
    % draw 'baseline'
    plot([which_TRs(1) which_TRs(end)]*TR,[0 0],'k-','LineWidth',0.75);
    
    % draw event markers
    mh(vv,:) = plot(t_markers.*[1;1],[0 .1],'-','Color',[0.7 0.7 0.7],'LineWidth',0.75);
    
    for hh = 1:2 
        
        thisd = nan(length(subj),size(all_HRFs{1},2));
        for ss = 1:length(subj)
            
            thisidx = all_subj == ss & all_ROIs == vv & all_conds(:,1)==cu(cc);% & floor(all_conds_task(:,1)/10)==which_conds(cc);
            
            thisd(ss,:) = nanmean(all_HRFs{hh}(thisidx,:),1); %using nanmean here 
            clear thisidx;
        end
        
   
        % plot, like for reconstructions, such that middle of each
        % datapoint is at middle of TR
        if hh ==1
          lh(hh,cc) = plot(which_TRs*TR + TR/2,nanmean(thisd,1),'-','LineWidth',1.0,'Color',cond_colors(cc,:));
        else
          lh(hh,cc) = plot(which_TRs*TR + TR/2,nanmean(thisd,1),'--','LineWidth',1.0,'Color',cond_colors(cc,:));

        end

        btwn_fill = [nanmean(thisd,1)+1.*nanstd(thisd,[],1)/sqrt(length(subj)) fliplr( nanmean(thisd,1)-1.*nanstd(thisd,[],1)/sqrt(length(subj)) )]; % geh test
        
        fill([which_TRs*TR + TR/2 fliplr(which_TRs*TR + TR/2)],btwn_fill,cond_colors(cc,:),'linestyle','none','facealpha',0.3);
 
        
        
        clear thisd;
    end
    
    title(ROIs{vv});
    
    if vv==1 && hh ==2
        ylabel({'Voxels in vs out RF'; 'BOLD Z-score'});
        xlabel('Time (s)');
        xticks([0 4.5 12])
        set(gca,'XTickLabel',[0 4.5 12])
    else
        xticks([0 4.5 12])
        set(gca,'XTickLabel',[0 4.5 12])

    end
%    
%     if vv==length(ROIs) && hh==2
%         legend(lh,sprintf('%s RFin',condstr{cc}),sprintf('%s RFout',condstr{cc}),'location','NorthEast');
%     else
%     end
        
    
end
myy = cell2mat(get(axhrf(cc,:),'YLim'));
set(axhrf(cc,:),'YLim',[min(myy(:,1)) max(myy(:,2))],'XLim',[which_TRs(1) which_TRs(end)]*TR,'TickDir','out');
set(mh,'YData',[min(myy(:,1)) max(myy(:,2))]);
if length(ROIs)==7
set(gcf,'Position',[ 429         451        2041         471])
else
set(gcf,'Position',[ 429         451        2541         471])  
end
 %legend(lh,sprintf('%s RFin',condstr{cc}),sprintf('%s RFout',condstr{cc}),'location','NorthEast');
%match_ylim(get(gcf,'Children'));
end
%% in 'dist_ang_all' alignment case, do the following plot 
if do_distalign_plot==1
    
    all_HRFs{1} = all_hrfs_in;
    all_HRFs{2} = all_hrfs_out;
    all_HRFs{3} = all_hrfs_in - all_hrfs_out;
    cond_colors = spDist_condColors;
    HRF_str ={'RFin','RFout','RFin-out'};
    % store something that's ROI x time x subj for each condition
    cu = unique(all_conds(:,1));
    
    axhrf = nan(length(all_HRFs),length(ROIs));
    mh = nan(length(ROIs),length(t_markers));
    
    condstr = {'Dist Absent','Dist Present'};
    
    figure;
    for cc = 2 %need only one condition here
        
        for vv = 1:length(ROIs)
            axhrf(cc,vv) = subplot(1,length(ROIs),vv); hold on;
            % draw 'baseline'
            plot([which_TRs(1) which_TRs(end)]*TR,[0 0],'k-','LineWidth',0.75);
            
            % draw event markers
            mh(vv,:) = plot(t_markers.*[1;1],[0 .1],'-','Color',[0.7 0.7 0.7],'LineWidth',0.75);
            
            for hh = 1:2 % only care about RFin, RFout here, therefore dont loop over length(all_mean_HRF)
                
                thisd = nan(length(subj),size(all_HRFs{1},2));
                for ss = 1:length(subj)
                    
                    thisidx = all_subj == ss & all_ROIs == vv & all_conds(:,1)==cu(cc);% & floor(all_conds_task(:,1)/10)==which_conds(cc);
                    
                    thisd(ss,:) = nanmean(all_HRFs{hh}(thisidx,:),1); %using nanmean here
                    clear thisidx;
                end
                
                
                % plot, like for reconstructions, such that middle of each
                % datapoint is at middle of TR
                if hh ==1
                    lh(hh,cc) = plot(which_TRs*TR + TR/2,nanmean(thisd,1),'-','LineWidth',1.0,'Color',cond_colors(cc,:));

                else
                     lh(hh,cc) = plot(which_TRs*TR + TR/2,nanmean(thisd,1),'--','LineWidth',1.0,'Color',cond_colors(cc,:));
                end
                
                btwn_fill = [nanmean(thisd,1)+1.*nanstd(thisd,[],1)/sqrt(length(subj)) fliplr( nanmean(thisd,1)-1.*nanstd(thisd,[],1)/sqrt(length(subj)) )];
                
                fill([which_TRs*TR + TR/2 fliplr(which_TRs*TR + TR/2)],btwn_fill,cond_colors(cc,:),'linestyle','none','facealpha',0.3);
                
                
                
                clear thisd btwn_fill fill
            end
            
            title(ROIs{vv});
            if vv==1 && hh ==2
                ylabel({'Voxels in vs out RF'; 'BOLD Z-score'});
                xlabel('Time (s)');
                xticks([0 4.5 12])
                set(gca,'XTickLabel',[0 4.5 12])
            else
                xticks([0 4.5 12])
                set(gca,'XTickLabel',[0 4.5 12])
                
            end
            
            
            
        end
        
        myy = cell2mat(get(axhrf(cc,:),'YLim'));
        set(axhrf(cc,:),'YLim',[min(myy(:,1)) max(myy(:,2))],'XLim',[which_TRs(1) which_TRs(end)]*TR,'TickDir','out');
        set(mh,'YData',[min(myy(:,1)) max(myy(:,2))]);
        if length(ROIs)==7
        set(gcf,'Position',[ 429         451        2041         223])
        else
         set(gcf,'Position',[  -71         407        2541         223])   
        end
        
    end
else
end
%% bar graph of mean delay period activity
% for hh =1:length(all_HRFs)
%     figure;
%     sgtitle(HRF_str{hh})
%     plot([0 length(ROIs)+1],[0 0],'k--');
%     
%     offsets = linspace(-0.15,0.15,length(all_mean_hrf));
%     
%     for cc = 1:length(cu)
%         % NOTE: here, we use lte rather than lt...
%         this_TR_range = which_TRs>=delay_tpt_range(1,1)&which_TRs<=delay_tpt_range(2,2);
%         
%         all_mean_delay = squeeze(mean(all_mean_hrf{hh,cc}(:,this_TR_range,:),2)); % ROI x subj
%         hold on;
%         for vv = 1:length(ROIs)
%             thise = std(all_mean_delay(vv,:),[],2)/sqrt(length(subj));
%             thism = mean(all_mean_delay(vv,:),2);
%             
%             plot(vv*[1 1]+offsets(cc),thism+thise*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(cc,:));
%             plot(vv+offsets(cc),thism,'o','Color',cond_colors(cc,:),'MarkerFaceColor','w','MarkerSize',8,'LineWidth',1.5);
%         end
%         
%         
%     end
%     
%     
%     set(gca,'XTick',1:length(ROIs),'XTickLabel',ROIs,'XTickLabelRotation',-45,'FontSize',14,'TickDir','out','Box','off');
%     ylabel({'BOLD Z-score';'Mean delay period activation'});
% end
%% plot, mean delay period activity RF in vs RF out, within distractor condition
RFmodfig = figure('Name','RFmod');

HRF_str ={'RFin','RFout','RfModIdx'};
store_mean =nan(2,2,length(ROIs));

for cc =1:length(cu)
    
    figure;
    sgtitle(condstr{cc})
    plot([0 length(ROIs)+1],[0 0],'k--');
    
    for hh = 1:length(all_mean_hrf)

        offsets = linspace(-0.15,0.15,length(all_mean_hrf));
        
        % NOTE: here, we use lte rather than lt...
        %   this_TR_range = which_TRs >=delay_tpt_range(1,1) & which_TRs<delay_tpt_range(2,2);
        
        this_TR_range = (which_TRs*TR)>=delay_tpt_range(1,1) & (which_TRs*TR)<delay_tpt_range(3,2);
        
        all_mean_delay = squeeze(mean(all_mean_hrf{hh,cc}(:,this_TR_range,:),2)); % ROI x subj
      
        for vv = 1:length(ROIs)
            thise = std(all_mean_delay(vv,:),[],2)/sqrt(length(subj));
            thism = mean(all_mean_delay(vv,:),2);
            if hh==1
                plot(vv*[1 1]+offsets(cc),thism+thise*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(cc,:));
                hold on;
                plot(vv+offsets(cc),thism,'o','Color',cond_colors(cc,:),'MarkerFaceColor',cond_colors(cc,:),'MarkerSize',8,'LineWidth',1.5);
                
            elseif hh==2
                plot(vv*[1 1]+offsets(cc),thism+thise*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(cc,:));
                plot(vv+offsets(cc),thism,'o','Color',cond_colors(cc,:),'MarkerFaceColor','w','MarkerSize',8,'LineWidth',1.5);
           
            elseif hh==3
                set(0, 'CurrentFigure', RFmodfig)
                hold on; 
                plot(vv*[1 1]+offsets(cc),thism+thise*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(cc,:));
                plot(vv+offsets(cc),thism,'o','Color',cond_colors(cc,:),'MarkerFaceColor','k','MarkerSize',8,'LineWidth',1.5);
                
            end
            
            store_mean(cc,hh,vv) = thism;

            
        end

        
        
    end
    
    
    set(gca,'XTick',1:length(ROIs),'XTickLabel',ROIs,'XTickLabelRotation',-45,'FontSize',14,'TickDir','out','Box','off');
    ylabel({'BOLD Z-score';'Mean delay period activation'});
end


% create a container for all labels ( ROI x DIST COND x RF COND x SUBJ) 
% corresponding data using all TRs in our window - from 3.25 to 12

data_allt   = nan(length(subj)*length(ROIs)*size(delay_tpt_range,1)*length(all_mean_hrf),1);
labels_allt = nan(length(subj)*length(ROIs)*size(delay_tpt_range,1)*length(all_mean_hrf),4); % subj, ROI, epochs, RFs

idx = 1;
for hh = 1:length(all_mean_hrf) % RF-in X RF-out X RF-in minus RF-out
    for ss = 1:length(subj)
        for vv = 1:length(ROIs)
            for cc = 1:length(cu)
     
                this_TR_range = (which_TRs*TR) >= delay_tpt_range(1,1) & (which_TRs*TR) < delay_tpt_range(3,2);
                
                
                data_allt(idx)     = mean(all_mean_hrf{hh,cc}(vv,this_TR_range,ss));
                
                
                labels_allt(idx,:) = [vv cc hh ss]; % ROI COND RF EPOCH SUBJ
                idx = idx+1;
                
            end
        end
    end
end


% REGRESSION ON RFIN, DISTRACTOR ABSENT (OR PRESENT) %SUBJECTS!!!
oneway_delayROI_realF = nan(2,1);
pval_roi= nan(2,1);
fitresult = cell( 3, 1 );
which_RF = 1;
for cc =1 %:length(cu)
    
    clear thisidx
    thisidx =  labels_allt(:,2)==cc & labels_allt(:,3)==which_RF; % change this for 1 or 3 depending - are we using mod? 
    myd = data_allt(thisidx);
    myrois = labels_allt(thisidx,1);
    [oneway_delayROI_realF(cc),pval_roi(cc)] = RMAOV1_gh([data_allt(thisidx) labels_allt(thisidx,[1 4])]);
    RMAOV1([data_allt(thisidx) labels_allt(thisidx,[1 4])]) % use OG version to print results to command line
    
    % Initialize arrays to store fits and goodness-of-fit.
 
    
    % Fit: 'distabsent_RFin'.
    [xData, yData] = prepareCurveData(  myrois, myd  );
    
    % Set up fittype and options. %USE THE 
   % ft = fittype( {'(sin(x-pi))', '((x-10)^2)', '1'}, 'independent', 'x', 'dependent', 'y', 'coefficients', {'a', 'b', 'c'} );
    ft = fittype( 'poly1' );
    % Fit model to data.
    gof = struct( 'sse', cell( 3, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
    [fitresult{cc}, gof(1)] = fit( xData, yData, ft );
    
    gof_store(cc) = gof(1).rsquare;
    % Plot fit with data.
    figure( 'Name', condstr{cc} );
    h = plot( fitresult{cc}, xData, yData );
    legend( h, 'amplitude ', condstr{cc}, 'Location', 'NorthEast' );
    % Label axes
    xlabel ROIs
    ylabel Amplitude
    grid on
    title(sprintf('R-square = %0.2f \n ',gof_store(cc)))
  
    clear myd myrois gof
    
end

% SAMPLE MEANS
% REGRESSION ON RFIN, DISTRACTOR ABSENT (OR PRESENT) ON INDIVIDUAL SUBJECTS
oneway_delayROIsubj_realF = nan(length(cu),length(subj));
pval_roisubj = nan(length(cu),length(subj));
betas = nan(length(cu),length(subj));
rsq_store = nan(length(cu),length(subj)); 
fitresult = cell(2,7);
which_RF = 1;

for cc =1 %:length(cu)
   for ss =1:length(subj)
       
    clear thisidx myd myrois xData yData
    thisidx =  labels_allt(:,2)==cc & labels_allt(:,3)==which_RF & labels_allt(:,4)==ss;  
    myd = data_allt(thisidx);
    myrois = labels_allt(thisidx,1);
    
    % Initialize arrays to store fits and goodness-of-fit.
   
    gof = struct( 'sse', cell( 3, 1 ), ...
        'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
    
    % Fit: 'distabsent_RFin'.
    [xData, yData] = prepareCurveData(  myrois, myd  );
    
    % Set up fittype and options. %USE THE 
    ft = fittype( 'poly1' );
    % Fit model to data.
    [fitresult{cc,ss}, gof(1)] = fit( xData, yData, ft);
    betas(cc,ss) = fitresult{cc,ss}.p1;
    rsq_store(cc,ss) = gof(1).rsquare;
    % Plot fit with data.
    figure( 'Name',sprintf('subj %s %s',subj{ss},condstr{cc}) );
    h = plot( fitresult{cc,ss}, xData, yData);
    legend( h, 'amplitude ', condstr{cc}, 'Location', 'NorthEast' );
    % Label axes
    xlabel ROIs
    ylabel Amplitude
    grid on
    title(sprintf('R-square = %0.2f \n ',rsq_store(cc,ss)))
    clear thisiidx xData yData
   end 
    [~,p_store(cc),~,stats] = ttest(betas(cc,:));
    realT(cc) = stats.tstat;
end

% shuffle roi labels for regression, within subj


iter=1000;
fitresult_shuf = cell( 2,length(subj),iter);
betas_shuf = nan(length(subj),iter);
labels_shuf = labels_allt;
which_RF = 1; % on which RF condition do we want to perform regression? 

for cc =1 % only performing this for distractor absent, built to add distractor present if needs change
        tic;
        
        for ii = 1:iter
            
            for ss = 1:length(subj)
                % cond, epoch, subj % ROI RF COND EPOCH SUBJ
                subjidx = labels_allt(:,2)==cc & labels_allt(:,3)==which_RF & labels_allt(:,4)==ss;
                thisidx = find(subjidx);
                shufidx = randperm(length(thisidx))';
                labels_shuf(subjidx,1) = labels_allt(thisidx(shufidx),1); %shuffle ROI labels for this subj
   
                myd = data_allt(subjidx); % maintain original data indices
                myrois = labels_shuf(subjidx,1); % shuffle only this ss labels 
                
                gof = struct( 'sse', cell( 3, 1 ), ...
                    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
                
                % Fit: 'distabsent_RFin'.
                [xData, yData] = prepareCurveData(  myrois, myd  );
                
                % Set up fittype and options. %USE THE
                ft = fittype( 'poly1' );
                % Fit model to data.
                [fitresult_shuf{cc,ss,ii}, gof(1)] = fit( xData, yData, ft);
                betas_shuf(ss,ii) = fitresult_shuf{cc,ss,ii}.p1;
                rsq_store_shuf(ss,ii) = gof(1).rsquare;
                
                clear xData yData gof my myrois subjidx thisidx shufidx ft
            end
 
            [~,p_store,~,stats] =ttest(betas_shuf(:,ii));
            shufT(ii) = stats.tstat;

        end
         toc
end
    

figure; hist(shufT,50)   
hold on; plot([realT realT], [ylim],'r','linewidth',2)

exact_P = 2 * min(mean(shufT(:) > realT), mean(shufT(:) < realT))

 


  delay_2way_realF = nan(2,3);
  
  % IV1 = ROI, IV2 = RF 
  for cc =1 %:length(cu)
      
      clear thisidx
      thisidx =  labels_allt(:,2)==cc & labels_allt(:,3)~=3;
      [delay_2way_realF(cc,:)] = RMAOV2_gh([data_allt(thisidx) labels_allt(thisidx,[1 3 4])]);
      RMAOV2([data_allt(thisidx) labels_allt(thisidx,[1 3 4])])
      
  end
  
  %%%
   
    % do shuffled 2-way ANOVA for each ROI
    iter=1000;
    delay_2way_shufF = nan(iter,3);
    fprintf('Shuffling - 2 way ANOVAs\n');
    
        % get data for this ROI
        thisidx =  labels_allt(:,2)==1 & labels_allt(:,3)~=3;
        thisdata   = data_allt(thisidx);
        thislabels = labels_allt(thisidx,[1 3 4]); % ROI RF SUBJ
        
        tic;
        for ii = 1:iter
            
            data_shuf = nan(size(thisdata));
            
            for ss = 1:length(subj)
                
                subjidx = thislabels(:,3)==ss;
                thisidx = find(subjidx);
                shufidx = randperm(length(thisidx)); 
                data_shuf(subjidx) = thisdata(thisidx(shufidx)); 

            end
             delay_2way_shufF(ii,:) = RMAOV2_gh([data_shuf thislabels]); %shuffling data, not labels ... 
        end   
        
        toc;
   
  
    

     % calculate p-values for 2-way ANOVA
     delay_2way_labels = {'ROI','RF','ROI x RF'};
     delay_2way_pval = nan(1,3);
for ee =1:3
     delay_2way_pval(ee) = mean( delay_2way_realF(1,ee) <=  delay_2way_shufF(:,ee)'); % 3 pvalues: COND, EPOCH, COND x EPOCH
end 
    
    % FDR thresholds
%     fdr_thresh = nan(3,size(epoch_2way_pval,3)); % # of RF conds x # main effects x # ROIs
%     fdr_mask = nan(3,size(epoch_2way_pval,3));
% 
%     for ii = 1:size(epoch_2way_pval,2) % this should be the # of effects, ie 3, correct across ROI
%         [fdr_thresh(ii,:), fdr_mask(ii,:)] = fdr(squeeze(epoch_2way_pval(ii,:)),0.05);
%     end
%  
%   
  %%%

  % ONE-WAY ANOVA PER ROI, WITH RF AS FACTOR - DO THIS 
  
  delay_1way_realF = nan(2,length(ROIs));
  delay_1way_realP = nan(2,length(ROIs));
    for cc =1:length(cu)
        for vv = 1:length(ROIs)
            clear thisidx
            thisidx = labels_allt(:,1)==vv & labels_allt(:,2)==cc & labels_allt(:,3)~=3; 
            [ delay_1way_realF(cc,vv),  delay_1way_realP(cc,vv)] = RMAOV1_gh([data_allt(thisidx) labels_allt(thisidx,[3 4])]);
        end
        figure;
        plot( delay_1way_realF(cc,:),'ko-') % take a look at F-stats
        title(condstr{cc})
        set(gca,'XTick',1:length(ROIs),'XTickLabel',ROIs,'XTickLabelRotation',-45,'FontSize',14,'TickDir','out','Box','off');
        ylabel({'F-values, one-way ANOVA, '})
    end
  
 % do shuffled 1-way ANOVA for each ROI
    
    delay_1way_shufF = nan(length(ROIs),iter);
    fprintf('Shuffling subj-wise RF labels- 1 way ANOVAs\n');
     
    for vv = 1:length(ROIs)
        
        fprintf('ROI %s\n',ROIs{vv});
        clear thisidx thisdata thislabels
        % get data for this ROI
        thisidx = labels_allt(:,1)==vv & labels_allt(:,2)==1 & labels_allt(:,3)~=3; % for this ROI, no distractor, not RFmod
        thisdata   = data_allt(thisidx);
        thislabels = labels_allt(thisidx,[3 4]); % RF SUBJ
        
        tic;
        for ii = 1:iter
            
            data_shuf = nan(size(thisdata));
            
            for ss = 1:length(subj)
                
                subjidx = thislabels(:,2)==ss;
                thisidx = find(subjidx);
                shufidx = randperm(length(thisidx)); 
                data_shuf(subjidx) = thisdata(thisidx(shufidx)); 
                
                clear subjidx thisidx shufidx 
            end
            
             delay_1way_shufF(vv,ii) = RMAOV1_gh([data_shuf thislabels]);
        end
        toc;
    end

    

     % calculate p-values for 2-way ANOVA
%     epoch_2way_labels = {'condition','epoch','condition x epoch'};
    delay_1way_pval = nan(1,length(ROIs));
% only care about distractor absent 
    for vv = 1:length(ROIs)
        delay_1way_pval(1,vv) = mean(delay_1way_realF(1,vv) <= delay_1way_shufF(vv,:)) % 
    end

    
    % FDR thresholds
    fdr_thresh = nan(1,size(delay_1way_pval,2)); % # of RF conds x # main effects x # ROIs
    fdr_mask = nan(1,size(delay_1way_pval,2));


    [fdr_thresh(:), fdr_mask(:)] = fdr(delay_1way_pval(:),0.05);

 %%%
%% mean epoch activation 

hoffset=0.15;
RF_str = {'RFin','RFout','RFinminusRFout'};
for hh =1:length(all_HRFs)
    clear thisd thisd_store
fh.(sprintf('%s',RF_str{hh})) = figure;
 sgtitle(HRF_str{hh})
for vv = 1:length(ROIs)
    
    subplot(1,length(ROIs),vv); hold on;
    
    % plot a zero line
    plot([0 size(delay_tpt_range,1)+1],[0 0],'k--');

        % epoch x subj
        thisd = nan(size(delay_tpt_range,1),length(subj));
        thisd_store = nan(2,size(delay_tpt_range,1),length(subj));

       for cc = 1:length(cu)
           
           for tt = 1:size(delay_tpt_range,1)
          
               
               this_TR_range = (which_TRs*TR)>=delay_tpt_range(tt,1) & (which_TRs*TR)<delay_tpt_range(tt,2);
               
               thisd_store(cc,tt,:) = squeeze(mean(all_mean_hrf{hh,cc}(vv,this_TR_range,:),2))';
               thisd(tt,:) = mean(all_mean_hrf{hh,cc}(vv,this_TR_range,:),2);
               
               % all_mean_hrf: ROIs x tpts x subj
               %all_mean_delay = squeeze(mean(all_mean_hrf{cc}(vv,this_TR_range,:),2)); % ROI x subj
           end
             
               thise = std(thisd,[],2)/sqrt(length(subj));
               thism = mean(thisd,2);
               
               
               plot((1:size(delay_tpt_range,1)) .* [1; 1],thism.'+thise.'.*[-1; 1],'-','LineWidth',1.0,'Color',cond_colors(cc,:));
               hold on;
               lh(cc,:) = plot(1:size(delay_tpt_range,1),thism,'o-','Color',cond_colors(cc,:),'MarkerFaceColor',cond_colors(cc,:),'MarkerSize',7,'LineWidth',1.0);
         
       end
           
    
       for tt = 1:size(delay_tpt_range,1)

         hold on;
         
         plot(tt+hoffset*[-1;1],squeeze(thisd_store(:,tt,:)),'-','LineWidth',0.25,'Color',[0.4 0.4 0.4]);

         plot(tt+hoffset*-1,squeeze(thisd_store(1,tt,:)),'.','markersize',10,'Color',cond_colors(1,:));
         
         plot(tt+hoffset*1,squeeze(thisd_store(2,tt,:)),'.','markersize',10,'Color',cond_colors(2,:));
           
  
       end
    
      
    hold off;
    if length(ROIs) ==7
        set(gcf,'Position',[ 429         451        2041         305])
    else
        set(gcf,'Position',[ 429         451        2541         305])
    end
    set(gca,'XTick',1:size(delay_tpt_range,1),'XTickLabel',[],'XTickLabelRotation',-45,'FontSize',12,'TickDir','out','Box','off','YTick',-0.5:.5:2.25,'XLim',[0.5 size(delay_tpt_range,1)+0.5]);
    if vv == 1
        set(gca,'XTickLabel',epoch_str,'XTickLabelRotation',-45,'YTickLabel',-.5:.5:1.5);
        ylabel('BOLD Z-score');
    else
    end
    title(ROIs{vv});
    ylim([-0.5 1.5])
%     if vv==length(ROIs) && cc ==2
%     legend(lh,sprintf('%s RFin',condstr{1}),sprintf('%s RFin',condstr{2}),'location','NorthEast');
%     else
%     end
    
end
%match_ylim(get(gcf,'Children'))
end

fprintf('pause')


%% new

hoffset=0.15;
RF_str = {'RFin','RFout','RFinminusRFout'};
RF_str = {'RFin','RFout'};
for cc =1:2 %length(all_HRFs)
 figure   
    for vv = 1:length(ROIs)
        
        subplot(1,length(ROIs),vv); hold on;
        
        % plot a zero line
        plot([0 size(delay_tpt_range,1)+1],[0 0],'k--');
        
        % epoch x subj
        thisd = nan(size(delay_tpt_range,1),length(subj));
        thisd_store = nan(2,size(delay_tpt_range,1),length(subj));
        
        for hh = 1:2
            
            for tt = 1:size(delay_tpt_range,1)
                
                
                this_TR_range = (which_TRs*TR)>=delay_tpt_range(tt,1) & (which_TRs*TR)<delay_tpt_range(tt,2);
                
                thisd_store(hh,tt,:) = squeeze(mean(all_mean_hrf{hh,cc}(vv,this_TR_range,:),2))';
                thisd(tt,:) = mean(all_mean_hrf{hh,cc}(vv,this_TR_range,:),2);
                
                % all_mean_hrf: ROIs x tpts x subj
                %all_mean_delay = squeeze(mean(all_mean_hrf{cc}(vv,this_TR_range,:),2)); % ROI x subj
            end
            
            thise = std(thisd,[],2)/sqrt(length(subj));
            thism = mean(thisd,2);
            
            if hh==1
                
                plot((1:size(delay_tpt_range,1)) .* [1; 1],thism.'+thise.'.*[-1; 1],'-','LineWidth',1.0,'Color',cond_colors(cc,:));
                hold on;
                lh(cc,:) = plot(1:size(delay_tpt_range,1),thism,'o-','Color',cond_colors(cc,:),'MarkerFaceColor',cond_colors(cc,:),'MarkerSize',7,'LineWidth',1.0);
            elseif hh==2
                plot((1:size(delay_tpt_range,1)) .* [1; 1],thism.'+thise.'.*[-1; 1],'-','LineWidth',1.0,'Color',cond_colors(cc,:));
               
                lh(cc,:) = plot(1:size(delay_tpt_range,1),thism,'o-','Color',cond_colors(cc,:),'MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.0);
            end
            
            %
            %        for tt = 1:size(delay_tpt_range,1)
            %
            %          hold on;
            %
            %          plot(tt+hoffset*[-1;1],squeeze(thisd_store(:,tt,:)),'-','LineWidth',0.25,'Color',[0.4 0.4 0.4]);
            %
            %          plot(tt+hoffset*-1,squeeze(thisd_store(1,tt,:)),'.','markersize',10,'Color',cond_colors(1,:));
            %
            %          plot(tt+hoffset*1,squeeze(thisd_store(2,tt,:)),'.','markersize',10,'Color',cond_colors(2,:));
            %
            %
            %        end
            
            
         
            if length(ROIs) ==7
                set(gcf,'Position',[ 429         451        2041         305])
            else
                set(gcf,'Position',[ 429         451        2541         305])
            end
            set(gca,'XTick',1:size(delay_tpt_range,1),'XTickLabel',[],'XTickLabelRotation',-45,'FontSize',12,'TickDir','out','Box','off','YTick',-0.5:.5:2.25,'XLim',[0.5 size(delay_tpt_range,1)+0.5]);
            if vv == 1
                set(gca,'XTickLabel',epoch_str,'XTickLabelRotation',-45,'YTickLabel',-.5:.5:1.5);
                ylabel('BOLD Z-score');
            else
            end
            title(ROIs{vv});
            ylim([-0.5 1.5])
            %     if vv==length(ROIs) && cc ==2
            %     legend(lh,sprintf('%s RFin',condstr{1}),sprintf('%s RFin',condstr{2}),'location','NorthEast');
            %     else
            %     end
            
        end
        %match_ylim(get(gcf,'Children'))
    end
end
    fprintf('pause')
    
    



%% let's do some stats!
%%% GEH modifications to TCS code, see spDist_plotHRFs_ERA_pRFvoxSelection.m for original construction %%%

if do_stats == 1
% GOAL pr conversation with CC on 01/21/21

% therefore, hold hh (RF case) constant by idx'ing based on hh ==1, e.g.
% all_HRFs is organized by RFin, RFout, RFinout, loop over these cases for
% stats
 
    data_all   = nan(length(subj)*length(ROIs)*size(delay_tpt_range,1)*length(all_mean_hrf),1);
    labels_all = nan(length(subj)*length(ROIs)*size(delay_tpt_range,1)*length(all_mean_hrf),5); % subj, ROI, epochs, RFs
    idx = 1;
  for hh = 1:2 %:length(all_mean_hrf) % RF-in X RF-out X RF-in minus RF-out 
    for ss = 1:length(subj)
        for vv = 1:length(ROIs)
            for cc = 1:length(cu)
                for tt = 1:size(delay_tpt_range,1)
                    
                    this_TR_range = (which_TRs*TR)>=delay_tpt_range(tt,1) & (which_TRs*TR)<delay_tpt_range(tt,2);
                    
                    data_all(idx)     = mean(all_mean_hrf{hh,cc}(vv,this_TR_range,ss));
        
                    
                    labels_all(idx,:) = [vv cc hh tt ss]; % ROI COND RF EPOCH SUBJ
                    idx = idx+1;
                end
            end
        end
    end
  end 
    % real ANOVAs
    
    % 3-way with interaction
    % returns: IV1, IV2, IV3, IV1xIV2, IV1xIV3, IV2xIV3, IV1xIV2xIV3
    % subj must be last column of X
    
    %what are our main effects? we want the same RF case (e.g. RFin, hh ==1), but want cc,tt,vv to vary 
    %IV output is the following: ROI, RF,EPOCH
    
    %here, hold hh (which RF condition, constant. so ANOVA will be
    %performed WITHIN RF cond%%
    
    epoch_3way_realF= nan(length(all_mean_hrf),7); % HRF conditions by # main effects and interactions 
    
    for cc =1:length(cu)
        clear idx;
        idx = labels_all(:,2)== cc; 
        epoch_3way_realF(cc,:) = RMAOV33_gh([data_all(idx) labels_all(idx,1) labels_all(idx,3) labels_all(idx,4) labels_all(idx,5)]); % re: labels_all; 1 =ROI; 2 = COND; 4 =EPOCH; 5 =SUBJ
    end 
    
    
    % 2-way for each ROI
    % factors IV1=RF, IV2 =epoch
    % output: Fvals: IV1, IV2, IV1xIV2
    epoch_2way_realF = nan(2,3,length(ROIs));
    for cc =1:length(cu)
        for vv = 1:length(ROIs)
            clear thisidx
            thisidx = labels_all(:,1)==vv & labels_all(:,2)==cc;
            epoch_2way_realF(cc,:,vv) = RMAOV2_gh([data_all(thisidx) labels_all(thisidx,[3 4 5])]);
        end
                figure;
        plot(epoch_2way_realF(cc,:),'ko-')
        title(condstr{cc})
        set(gca,'XTick',1:length(ROIs),'XTickLabel',ROIs,'XTickLabelRotation',-45,'FontSize',14,'TickDir','out','Box','off');
        ylabel({'F-values, one-way ANOVA, '})
    end


    
    % want to see the F value for RF as the effect over all ROIs
    epoch_1way_realF = nan(2,length(ROIs));
    for cc =1:length(cu)
        for vv = 1:length(ROIs)
            clear thisidx
            thisidx = labels_all(:,1)==vv & labels_all(:,2)==cc;
            [epoch_1way_realF(cc,vv),~] = RMAOV1_gh([data_all(thisidx) labels_all(thisidx,[3 5])]);
        end
        figure;
        plot(epoch_1way_realF(cc,:),'ko-')
        title(condstr{cc})
        set(gca,'XTick',1:length(ROIs),'XTickLabel',ROIs,'XTickLabelRotation',-45,'FontSize',14,'TickDir','out','Box','off');
        ylabel({'F-values, one-way ANOVA, '})
    end
    
          oneway3_realF = nan(2,length(ROIs));
     %ROI is the factor
    for cc =1:length(cu)
       % for vv = 1:length(ROIs)
       clear thisidx
       thisidx =  labels_all(:,2)==cc & labels_all(:,4)==3 & labels_all(:,3)==1;
       [oneway3_realF(cc),~] = RMAOV1_gh([data_all(thisidx) labels_all(thisidx,[1 5])]);

    end 
    
   figure;
    for cc =1 %:length(cu)
        for vv = 1:length(ROIs)
            clear thisidx
            thisidx = labels_all(:,1)==vv & labels_all(:,2)==cc & labels_all(:,3)==1 & labels_all(:,4)==3;
            myd(cc,vv) = mean(data_all(thisidx));
        end
        hold on;
         mysem = std(myd(cc,:),[],2)/sqrt(length(subj))
         plot(myd(cc,:),'ko-');
         hold on;
         plot([1:length(ROIs)],myd(cc,:)+1.*mysem,'k')
         plot([1:length(ROIs)],myd(cc,:)-1.*mysem,'k')
    end
       
    set(gca,'XTick',1:length(ROIs),'XTickLabel',ROIs,'XTickLabelRotation',-45,'FontSize',14,'TickDir','out','Box','off');
    ylabel({'BOLD Z-score';'Mean epoch 3 period activation'})
    
    
    
       figure;
    for cc =1 %:length(cu)
        for vv = 1:length(ROIs)
            clear thisidx
            thisidx = labels_all(:,1)==vv & labels_all(:,2)==cc & labels_all(:,3)==3 & labels_all(:,4)==3;
            myd(cc,vv) = mean(data_all(thisidx));
        end
        hold on;
         mysem = std(myd(cc,:),[],2)/sqrt(length(subj))
         plot(myd(cc,:),'ko-');
         hold on;
         plot([1:length(ROIs)],myd(cc,:)+1.*mysem,'k')
         plot([1:length(ROIs)],myd(cc,:)-1.*mysem,'k')
    end
       
    set(gca,'XTick',1:length(ROIs),'XTickLabel',ROIs,'XTickLabelRotation',-45,'FontSize',14,'TickDir','out','Box','off');
    ylabel({'BOLD Z-score, Rf in minus Rf out';'Mean epoch 3 period activation'})
    
    

    
    niter = 1000;
    
    % do shuffled 3-way ANOVA (shuffle labels w/in subj 1000x)
    
    epoch_3way_shufF = nan(2,niter,7);
    fprintf('Shuffling (3-way ANOVA...)\n');
    
    tic;
    for ii = 1:niter
        
        data_shuf = nan(size(data_all));
        
        for ss = 1:length(subj)
            subjidx = labels_all(:,5)==ss;
            thisidx = find(subjidx);
            shufidx = randperm(length(thisidx));
            data_shuf(subjidx) = data_all(thisidx(shufidx));
            
            clear thisidx subjidx shufidx;
        end
        
        
        for hh =1:length(all_mean_hrf)
        clear idx;
        idx = labels_all(:,3)== hh; 
        epoch_3way_shufF(hh,ii,:) = RMAOV33_gh([data_shuf(idx) labels_all(idx,1) labels_all(idx,2) labels_all(idx,4) labels_all(idx,5)]);
        end 
    end
    toc;
    
    
    
    % do shuffled 2-way ANOVA for each ROI
    
    epoch_2way_shufF = cell(2, length(ROIs));
    fprintf('Shuffling - 2 way ANOVAs\n');
    
   for hh=1:length(all_mean_hrf)
       
    for vv = 1:length(ROIs)
        
        fprintf('ROI %s\n',ROIs{vv});
        
        % get data for this ROI
        thisidx = labels_all(:,1)==vv & labels_all(:,3)==hh;
        thisdata   = data_all(thisidx);
        thislabels = labels_all(thisidx,[2 4 5]); % cond, epoch, subj % ROI RF COND EPOCH SUBJ
        
        tic;
        for ii = 1:niter
            
            data_shuf = nan(size(thisdata));
            
            for ss = 1:length(subj)
                
                subjidx = thislabels(:,3)==ss;
                thisidx = find(subjidx);
                shufidx = randperm(length(thisidx)); 
                data_shuf(subjidx) = thisdata(thisidx(shufidx)); 
               
                
            end
             epoch_2way_shufF{hh,vv}(ii,:) = RMAOV2_gh([data_shuf thislabels]);
        end
        toc;
    end
   end
    
    % calculate p-values for 3-way ANOVA
    
    for hh=1:length(all_mean_hrf)
        epoch_3way_pval(hh,:) = mean(epoch_3way_realF(hh,:) <= squeeze(epoch_3way_shufF(hh,:,:)),1);
    end
    
     epoch_3way_labels = {'ROI','RF','epoch','ROI x RF','ROI x epoch','RF x epoch','ROI x RF x epoch'};
    
     % calculate p-values for 2-way ANOVA
    epoch_2way_labels = {'condition','epoch','condition x epoch'};
    epoch_2way_pval = nan(2,3,length(ROIs));
    
    for hh =1:length(all_mean_hrf)
    for vv = 1:length(ROIs)
        epoch_2way_pval(hh,:,vv) = mean(epoch_2way_realF(hh,:,vv) <= epoch_2way_shufF{hh,vv}); % 3 pvalues: COND, EPOCH, COND x EPOCH
    end
    end
    
    % FDR thresholds
    fdr_thresh = nan(3,3,size(epoch_2way_pval,3)); % # of RF conds x # main effects x # ROIs
    fdr_mask = nan(3,3,size(epoch_2way_pval,3));
    for hh=1:length(all_mean_hrf)
    for ii = 1:size(epoch_2way_pval,2) % this should be the # of effects, ie 3, correct across ROI
        [fdr_thresh(hh,ii,:), fdr_mask(hh,ii,:)] = fdr(squeeze(epoch_2way_pval(hh,ii,:)),0.05);
    end
    end 
    % save stats, if necessary
    
    if save_stats == 1
        fn2s = sprintf('%s/spDist_stats/n%i_HRFstats_shuf_%iIter_%s.mat',root,length(subj),niter,datestr(now,30));
        fprintf('Saving to %s\n',fn2s);
        save(fn2s,'subj','ROIs','epoch_3way_realF','epoch_3way_shufF','epoch_2way_realF','epoch_2way_shufF','epoch_3way_pval','epoch_2way_pval','epoch_3way_labels','epoch_2way_labels','fdr_thresh');
    end
    
    
    % add stats markers to figures
    % cond: +, epoch: o, interaction: x
    
    %stats_markers = {'+','\circ','\times'}; %RF EPOCH RFxEPOCH 
     stats_markers = {'C','E','X'}
    stats_pos = [-1.30 -.03; -1.0 -.04; -.75 -.03]; % in plotted units, x,y, relative to top right
    
    signif_color = [0 0 0];
    trend_color  = [0.5 0.5 0.5];
for hh =1:length(all_mean_hrf) 
    figure(fh.(sprintf('%s',RF_str{hh})));
    for vv = 1:length(ROIs)
        
        subplot(1,length(ROIs),vv); hold on;
        
        tmpxlim = get(gca,'XLim');
        tmpylim = get(gca,'YLim');
        
        topright = [tmpxlim(2) tmpylim(2)];
        
        % main effect of cond; epoch, interaction
        for ii = 1:3
            if epoch_2way_pval(hh,ii,vv) <= fdr_thresh(hh,ii,vv)
                text(topright(1)+stats_pos(ii,1),topright(2)+stats_pos(ii,2),stats_markers{ii},'FontSize',18,'Color',signif_color);
            elseif epoch_2way_pval(hh,ii,vv) <= 0.05
                text(topright(1)+stats_pos(ii,1),topright(2)+stats_pos(ii,2),stats_markers{ii},'FontSize',18,'Color',trend_color);
            end
        end
        
    end
end 
end 

fprintf('\n stats complete\n')

return