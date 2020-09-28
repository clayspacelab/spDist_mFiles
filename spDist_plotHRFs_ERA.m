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

function spDist_plotHRFs_ERA(subj,sess,ROIs)


task_dir = 'spDist';

root = spDist_loadRoot;

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
    % ROIs = {'V1','V2','V3','V3AB','VO1','VO2','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'};
    
    % for manuscript
    ROIs = {'V1','V2','V3','V3AB','hV4','LO1','IPS0','IPS1','IPS2','IPS3','sPCS'};
end

% do stats? (takes ~10 mins)
do_stats = 1;

% if so, save them?
save_stats = 1;



func_suffix = 'surf';

ve_thresh = 0.1; 

% 'overall' delay-period (maximizes difference between distractor
% present/absent
% (because we use less than/equal to for tpt selection...)
% (this is based on tr - we want to INCLUDE TRs that span this range)
delay_range = [9 15]; % TR beginning at 7, ending at 15, to match IEM training

% these are based on *time*, not TR...
% as in Fig. 1, 4, 5 - 3 epochs for comparing before/during/after distractor
delay_tpt_range = [3.75 5.25; 7.5 9; 10.5 12]; %for stats - up for discussion?
% >= n,1 and < n,2

epoch_str = {'Pre-distractor','Distractor','Post-distractor'};

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
                
                all_hrfs = nan(nblank,size(data.dt_all,3));
                all_conds = nan(nblank,size(data.c_all,2));
                
                all_subj = nan(nblank,1);
                all_ROIs = nan(nblank,1);
                
                all_sess = nan(nblank,1);
                
                
                TR = data.TR;
                which_TRs = data.which_TRs;
                
            end
            
            
            thisidx = startidx:(startidx+size(data.dt_all,1)-1);
            
            which_vox = data.rf.ve>=ve_thresh;
            
            all_hrfs(thisidx,:) = squeeze(mean(data.dt_allz(:,which_vox,:),2));
            
            all_conds(thisidx,:) = data.c_all;
            
            all_subj(thisidx) = ss;
            
            all_ROIs(thisidx) = vv;
            
            all_sess(thisidx) = sess_idx;
            
            startidx = thisidx(end)+1;
            
            clear data;
            
        end
    end
    
end

%% remove baseline
baseline_TRs = which_TRs < 0;
all_hrfs = all_hrfs - mean(all_hrfs(:,baseline_TRs),2);


%% plot data



cond_colors = spDist_condColors; %[0.7100 0.2128 0.4772; 0 0 1;];

% store something that's ROI x time x subj for each condition
cu = unique(all_conds(:,1));
all_mean_hrf = cell(length(cu),1); % mapping, cued, chooose 

axhrf = nan(1,length(ROIs));
mh = nan(length(ROIs),length(t_markers));



% start with just subplots of timeseries (w/ errorbars?)
figure;
for vv = 1:length(ROIs)
    
    axhrf(1,vv) = subplot(1,length(ROIs),vv); hold on;
    
    % draw 'baseline'
    plot([which_TRs(1) which_TRs(end)]*TR,[0 0],'k-','LineWidth',0.75);
    
    % draw event markers
    mh(vv,:) = plot(t_markers.*[1;1],[0 .1],'-','Color',[0.7 0.7 0.7],'LineWidth',0.75);
    
    for cc = 1:length(cu)
        
        thisd = nan(length(subj),size(all_hrfs,2));
        for ss = 1:length(subj)
            
            thisidx = all_subj == ss & all_ROIs == vv & all_conds(:,1)==cu(cc);% & floor(all_conds_task(:,1)/10)==which_conds(cc);
            
            thisd(ss,:) = mean(all_hrfs(thisidx,:),1);
            all_mean_hrf{cc}(vv,:,ss) = mean(all_hrfs(thisidx,:),1);
            clear thisidx;
        end
        
        
        
        plot(which_TRs*TR,mean(thisd,1),'-','LineWidth',1.5,'Color',cond_colors(cc,:));
        plot(which_TRs*TR,mean(thisd,1)+std(thisd,[],1)/sqrt(length(subj)),':','LineWidth',0.75,'Color',cond_colors(cc,:));
        plot(which_TRs*TR,mean(thisd,1)-std(thisd,[],1)/sqrt(length(subj)),':','LineWidth',0.75,'Color',cond_colors(cc,:));
        clear thisd;
    end
    
    title(ROIs{vv});
    if vv == 1
        ylabel('BOLD Z-score');
        xlabel('Time (s)');
    else
        set(gca,'YTickLabel',[]);
    end
    
end

myy = cell2mat(get(axhrf,'YLim'));
set(axhrf,'YLim',[min(myy(:,1)) max(myy(:,2))],'XLim',[which_TRs(1) which_TRs(end)]*TR,'TickDir','out');
set(mh,'YData',[min(myy(:,1)) max(myy(:,2))]);

%set(gcf,'Position',[ 66        1180        2279         158])
%legend(condstr,'location','best');


%% also individual subj...
figure;
ax_subjhrf = nan(length(subj),length(ROIs)); % axis handles
for ss = 1:length(subj)
    for vv = 1:length(ROIs)
        
        ax_subjhrf(ss,vv) = subplot(length(subj),length(ROIs),(ss-1)*length(ROIs)+vv);
        hold on;
        
        for cc = 1:length(all_mean_hrf)
            plot(which_TRs*TR,all_mean_hrf{cc}(vv,:,ss),'-','LineWidth',1.25,'Color',cond_colors(cc,:));
            
        end
        
        if ss == 1
            title(ROIs{vv});
        end
        
        if vv == 1
            ylabel(sprintf('%s (BOLD Z-score)',subj{ss}));
        else 
            set(gca,'YTicKLabel',[]);
        end
        
        if ss == length(subj) && vv == 1
            xlabel('Time (s)');
            
        end
        
        
    end
end

match_ylim(ax_subjhrf(:));
set(ax_subjhrf(:),'XLim',[which_TRs(1) which_TRs(end)]*TR,'TickDir','out');

%% a top-down plot like fidelity
%cond_str = {'Mapping','R2-cued','R2-choose'};
cond_str = {'No distractor','Distractor'};
figure; 
for cc = 1:length(all_mean_hrf)
    subplot(1,length(all_mean_hrf),cc);
    
    this_TRs = which_TRs;
    imagesc(this_TRs*TR, 1:length(ROIs),mean(all_mean_hrf{cc},3));
    colormap viridis;
    set(gca,'YTick',1:length(ROIs),'YTickLabel',ROIs,'TickDir','out','Box','off','FontSize',14);
    xlabel('Time (s)');
    title(cond_str{cc});
    %tmp_clim = get(gca,'CLim');
    %set(gca,'CLim',[-1 1] * max(abs(tmp_clim)),'TickDir','out','box','off','FontSize',14);
end
match_clim(get(gcf,'Children'));

%% bar graph of mean delay period activity
figure;

plot([0 length(ROIs)+1],[0 0],'k--');

offsets = linspace(-0.15,0.15,length(all_mean_hrf));

for cc = 1:length(all_mean_hrf)
    % NOTE: here, we use lte rather than lt...
    this_TR_range = which_TRs>=delay_range(1)&which_TRs<=delay_range(2);
    
    all_mean_delay = squeeze(mean(all_mean_hrf{cc}(:,this_TR_range,:),2)); % ROI x subj
    hold on;
    for vv = 1:length(ROIs)
        thise = std(all_mean_delay(vv,:),[],2)/sqrt(length(subj));
        thism = mean(all_mean_delay(vv,:),2);
        
        plot(vv*[1 1]+offsets(cc),thism+thise*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(cc,:));
        plot(vv+offsets(cc),thism,'o','Color',cond_colors(cc,:),'MarkerFaceColor','w','MarkerSize',8,'LineWidth',1.5);
    end
    

end


set(gca,'XTick',1:length(ROIs),'XTickLabel',ROIs,'XTickLabelRotation',-45,'FontSize',14,'TickDir','out','Box','off');
title('Mean delay period activation');
ylabel('BOLD Z-score');


%% bar graph of mean activation during each epoch
fh_epoch = figure;


for vv = 1:length(ROIs)
    
    subplot(1,length(ROIs),vv); hold on;
    
    % plot a zero line
    plot([0 size(delay_tpt_range,1)+1],[0 0],'k--');
    
    
    
    
    for cc = 1:length(all_mean_hrf)
        
        % epoch x subj
        thisd = nan(size(delay_tpt_range,1),length(subj));
        
        for tt = 1:size(delay_tpt_range,1)
        
            this_TR_range = (which_TRs*TR)>=delay_tpt_range(tt,1) & (which_TRs*TR)<delay_tpt_range(tt,2);
            
            thisd(tt,:) = mean(all_mean_hrf{cc}(vv,this_TR_range,:),2);
            
            % all_mean_hrf: ROIs x tpts x subj
            %all_mean_delay = squeeze(mean(all_mean_hrf{cc}(vv,this_TR_range,:),2)); % ROI x subj
            
            
        end
        
        thise = std(thisd,[],2)/sqrt(length(subj));
        thism = mean(thisd,2);

        
        plot((1:size(delay_tpt_range,1)) .* [1; 1],thism.'+thise.'.*[-1; 1],'-','LineWidth',1.0,'Color',cond_colors(cc,:));
        plot(1:size(delay_tpt_range,1),thism,'o-','Color',cond_colors(cc,:),'MarkerFaceColor',cond_colors(cc,:),'MarkerSize',3,'LineWidth',1.0);
    end
    
    hold off;
    
    
    set(gca,'XTick',1:size(delay_tpt_range,1),'XTickLabel',[],'XTickLabelRotation',-45,'FontSize',12,'TickDir','out','Box','off','YTick',-0.5:.25:1,'XLim',[0.5 size(delay_tpt_range,1)+0.5]);
    if vv == 1
        %set(gca,'XTickLabel',epoch_str,'XTickLabelRotation',-45,'YTickLabel',-.5:.25:.1);
        %set(gca,'XTickLabel',[],'XTickLabelRotation',-45,'YTickLabel',-.5:.25:.1);
        ylabel('BOLD Z-score');
    else
        set(gca,'YTickLabel',[]);
    end
    title(ROIs{vv});
    
end


match_ylim(get(gcf,'Children'));


%% let's do some stats!
%
% start with: 3-way ANOVA, incl 3-way interaction, against shuffled null
%
% then, 2-way ANOVAs per ROI, against shuffled null

% organize data into tall format expected by ANOVA functions


if do_stats == 1
    
    
    data_all   = nan(length(subj)*length(ROIs)*size(delay_tpt_range,1)*length(all_mean_hrf),1);
    labels_all = nan(length(subj)*length(ROIs)*size(delay_tpt_range,1)*length(all_mean_hrf),4); % subj, ROI, conds, epochs
    idx = 1;
    
    for ss = 1:length(subj)
        for vv = 1:length(ROIs)
            for cc = 1:length(all_mean_hrf)
                for tt = 1:size(delay_tpt_range,1)
                    
                    this_TR_range = (which_TRs*TR)>=delay_tpt_range(tt,1) & (which_TRs*TR)<delay_tpt_range(tt,2);
                    
                    data_all(idx)     = mean(all_mean_hrf{cc}(vv,this_TR_range,ss));
                    
                    
                    labels_all(idx,:) = [vv cc tt ss];
                    idx = idx+1;
                end
            end
        end
    end
    
    % real ANOVAs
    
    % 3-way with interaction
    % returns: IV1, IV2, IV3, IV1xIV2, IV1xIV3, IV2xIV3, IV1xIV2xIV3
    % subj must be last column of X
    epoch_3way_realF = RMAOV33_gh([data_all labels_all]);
    
    % 2-way for each ROI
    epoch_2way_realF = nan(length(ROIs),3);
    for vv = 1:length(ROIs)
        thisidx = labels_all(:,1)==vv;
        epoch_2way_realF(vv,:) = RMAOV2_gh([data_all(thisidx) labels_all(thisidx,2:end)]);
    end
    
    
    % seed RNG
    rng(spDist_randSeed);
    
    niter = 1000;
    
    % do shuffled 3-way ANOVA (shuffle labels w/in subj 1000x)
    
    epoch_3way_shufF = nan(7,niter);
    fprintf('Shuffling (3-way ANOVA...)\n');
    tic;
    for ii = 1:niter
        
        data_shuf = nan(size(data_all));
        
        for ss = 1:length(subj)
            subjidx = labels_all(:,4)==ss;
            thisidx = find(subjidx);
            shufidx = randperm(length(thisidx));
            
            
            data_shuf(subjidx) = data_all(thisidx(shufidx));
            
            clear thisidx subjidx shufidx;
        end
        
        epoch_3way_shufF(:,ii) = RMAOV33_gh([data_shuf labels_all]);
        
    end
    toc;
    
    
    
    % do shuffled 2-way ANOVA for each ROI
    
    epoch_2way_shufF = cell(length(ROIs),1);
    fprintf('Shuffling - 2 way ANOVAs\n');
    for vv = 1:length(ROIs)
        
        fprintf('ROI %s\n',ROIs{vv});
        
        % get data for this ROI
        thisdata   = data_all(labels_all(:,1)==vv);
        thislabels = labels_all(labels_all(:,1)==vv,[2 3 4]); % cond, epoch, subj
        
        epoch_2way_shufF{vv} = nan(3,niter);
        tic;
        for ii = 1:niter
            
            data_shuf = nan(size(thisdata));
            
            for ss = 1:length(subj)
                
                subjidx = thislabels(:,3)==ss;
                thisidx = find(subjidx);
                shufidx = randperm(length(thisidx));
                
                data_shuf(subjidx) = thisdata(thisidx(shufidx));
                
                epoch_2way_shufF{vv}(:,ii) = RMAOV2_gh([data_shuf thislabels]);
                
            end
        end
        toc;
    end
    
    
    % calculate p-values for 3-way ANOVA
    epoch_3way_pval = mean(epoch_3way_realF.' <= epoch_3way_shufF,2);
    epoch_3way_labels = {'ROI','condition','epoch','ROI x condition','ROI x epoch','condition x epoch','ROI x condition x epoch'};
    
    
    % calculate p-values for 2-way ANOVA
    epoch_2way_labels = {'condition','epoch','condition x epoch'};
    
    epoch_2way_pval = nan(length(ROIs),3);
    for vv = 1:length(ROIs)
        epoch_2way_pval(vv,:) = mean(epoch_2way_realF(vv,:) <= epoch_2way_shufF{vv}.',1);
    end
    
    
    % FDR thresholds
    fdr_thresh = nan(1,size(epoch_2way_pval,2));
    for ii = 1:size(epoch_2way_pval,2)
        fdr_thresh(ii) = fdr(epoch_2way_pval(:,ii),0.05);
    end
    
    % save stats, if necessary
    
    if save_stats == 1
        fn2s = sprintf('%s/spDist_stats/n%i_HRFstats_shuf_%iIter_%s.mat',root,length(subj),niter,datestr(now,30));
        fprintf('Saving to %s\n',fn2s);
        save(fn2s,'epoch_3way_realF','epoch_3way_shufF','epoch_2way_realF','epoch_2way_shufF','epoch_3way_pval','epoch_2way_pval','epoch_3way_labels','epoch_2way_labels','fdr_thresh');
    end
    
    
    % add stats markers to figures
    % cond: +, epoch: o, interaction: x
    
    stats_markers = {'+','\circ','\times'};
    stats_pos = [-1.30 -.03; -1.0 -.04; -.75 -.03]; % in plotted units, x,y, relative to top right
    
    signif_color = [0 0 0];
    trend_color  = [0.5 0.5 0.5];
    
    figure(fh_epoch);
    for vv = 1:length(ROIs)
        
        subplot(1,length(ROIs),vv); hold on;
        
        tmpxlim = get(gca,'XLim');
        tmpylim = get(gca,'YLim');
        
        topright = [tmpxlim(2) tmpylim(2)];
        
        % main effect of cond; epoch, interaction
        for ii = 1:3
            if epoch_2way_pval(vv,ii) <= fdr_thresh(ii)
                text(topright(1)+stats_pos(ii,1),topright(2)+stats_pos(ii,2),stats_markers{ii},'FontSize',14,'Color',signif_color);
            elseif epoch_2way_pval(vv,1) <= 0.05
                text(topright(1)+stats_pos(ii,1),topright(2)+stats_pos(ii,2),stats_markers{ii},'FontSize',14,'Color',trend_color);
            end
        end
        
    end
end

return