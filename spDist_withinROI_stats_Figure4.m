

function spDist_withinROI_stats_Figure4(subj,sess,ROIs)

root = spDist_loadRoot;%'/share/data/spDist/'; %updated root

task_dir = 'spDist';

if nargin < 1 || isempty(subj)
    subj = {'AY','CC','EK','KD','MR','SF','XL'};
end

if nargin < 2 || isempty(sess)
    sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
end

if nargin < 3 || isempty(ROIs)
    ROIs = {'V1','V2','V3','V3AB','hV4','LO1','IPS0','IPS1','IPS2','IPS3','sPCS'};
end

%roi_str = {'V1','V2','V3','V3AB','hV4','LO1','IPS0','IPS1','IPS2','IPS3','sPCS'};
func_suffix = 'surf';

plot_hists = 0; % whether to plot histograms of shuffled stat values or not...

nchan = 8;
which_vox = 0.1; % if > 1 , look for WHICH_VOXvox string; otherwise, look for VE<100*WHICH_VOX>

smooth_by = 1; % if this is 1, use regular files, otherwise, load smooth_by files

myTR = 0.75;

%t_range_to_plot = [-inf 12]; % plot b/w these (s)

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


delay_tpt_range = [3.75 5.25; 8 9.5; 10.5 12]; %for stats - up for discussion?

% number of iterations for permutation...
iter = 1000;


% seed random number generator:     
rng(spDist_randSeed);



%% load data
startidx = 1;
for ss = 1:length(subj)
    
    for vv = 1:length(ROIs)
        
        
        % just one file to load
        fn = sprintf('%s/%s_reconstructions/%s_%s_%s_%s_%ichan%s%s%s_recon_thruTime1.mat',root,task_dir,subj{ss},horzcat(sess{ss}{:}),ROIs{vv},func_suffix,nchan,vox_str,smooth_str,trn_str);
        
        fprintf('loading %s...\n',fn);
        data_ru = load(fn);
        
        
        if vv == 1 && ss == 1
            % initialize variables...
            
            
            nblankt = length(ROIs)*size(data_ru.recons{1},1);
            all_recons = cell(size(data_ru.recons));
            for aa = 1:length(data_ru.recons)
                all_recons{aa} = nan(nblankt,size(data_ru.recons{aa},2),size(data_ru.recons{aa},3));
            end
            
            
            all_conds = nan(nblankt,size(data_ru.c_all,2));
            all_angs = nan(nblankt,size(data_ru.a_all,2));
            
            all_fidelity = nan(nblankt,size(data_ru.recons{1},3),length(data_ru.recons)); % timecoruse of fidelity for each alignment condition
            
            all_subj = nan(nblankt,1);
            all_ROIs = nan(nblankt,1);
            all_sess = nan(nblankt,1);
            
            all_r = nan(nblankt,1);
            
            
            angs = data_ru.angs;
            tpts = data_ru.delay_tpts;
            
            
        end
        
        
        
        thisidx = startidx:(startidx+size(data_ru.c_all,1)-1);
        
        for aa = 1:length(all_recons)
            all_recons{aa}(thisidx,:,:) = data_ru.recons{aa};
            all_fidelity(thisidx,:,aa) = squeeze(mean(cosd(angs) .* data_ru.recons{aa},2));
        end
        
        
        all_conds(thisidx,:) = data_ru.c_all;
        all_angs(thisidx,:) = data_ru.a_all;
        
        
        all_subj(thisidx) = ss;
        
        
        all_ROIs(thisidx) = vv;
        
        all_sess(thisidx) = data_ru.sess_all;
        
        % run and session combined for shuffled t-test within epoch
        all_r(thisidx) = data_ru.sess_all*100+data_ru.r_all;
        
        
        startidx = thisidx(end)+1;
        
        clear data;
        
        
    end
    
end


%% which tpts are we plotting throughout?

delay_tpts = cell(size(delay_tpt_range,1),1);
for dd = 1:size(delay_tpt_range,1)
    delay_tpts{dd} = (tpts*myTR) >= delay_tpt_range(dd,1) & (tpts*myTR) < delay_tpt_range(dd,2);
end

%% do this, over average delay epochs, for each condition
cu = [1 2];
dist_colors = spDist_condColors;%     [0.7100 0.2128 0.4772; 0 0 1;]; %first is RED, second is blue, beware its opposite matlab (lines(2)). this color scheme is in alignment w precue colors in actual exp (mag - mgs, cyan -dist)
cond_colors = dist_colors;

reconcomp = figure('name','reconcomp');

for dd = 1:length(delay_tpts)
    for vv = 1:length(ROIs)
        subplot(length(delay_tpts),length(ROIs),(dd-1)*length(ROIs)+vv); hold on;
        
        for cc = 1:length(cu)
            thisdata = nan(length(subj),length(angs));
            for ss = 1:length(subj)
                thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==cu(cc);
                thisdata(ss,:) = mean(mean(all_recons{1}(thisidx,:,delay_tpts{dd}),3),1);
                
            end
            
            % added sqrt here (TCS)
            my_sem = std(thisdata,[],1)/sqrt(length(subj));
            % TCS: plot actual values, not min-removed
            h(cc) = plot(linspace(-180,180,90),mean(thisdata,1),'-','LineWidth',0.5,'color',cond_colors(cc,:));
            % orig:
            %h(cc) = plot(linspace(-180,180,90),mean(thisdata,1) - min(mean(thisdata,1)),'-','LineWidth',0.5,'color',cond_colors(cc,:));
            hold on;
            %plot(linspace(-180,180,90),(mean(thisdata,1) - min(mean(thisdata,1)))+1.*my_sem,'-','LineWidth',.1,'color',cond_colors(cc,:),'HandleVisibility','off')
            
            %plot(linspace(-180,180,90),(mean(thisdata,1) - min(mean(thisdata,1)))-1.*my_sem,'-','LineWidth',.1,'color',cond_colors(cc,:),'HandleVisibility','off')
            
            %color in btwn +/- SEM
            % TCS: plot actual values, not min removed
            btwn_fill = [(mean(thisdata,1))+1.*my_sem fliplr((mean(thisdata,1))-1.*fliplr(my_sem))]; % TCS: fliplr also for my_sem?
            
            % orig:
            %btwn_fill = [(mean(thisdata,1) - min(mean(thisdata,1)))+1.*my_sem fliplr((mean(thisdata,1) - min(mean(thisdata,1)))-1.*my_sem)];
            fill([linspace(-180,180,90) fliplr(linspace(-180,180,90))],btwn_fill,cond_colors(cc,:),'linestyle','none','facealpha',0.3);
            
            
            %line([0 0], [-.1 2.25], 'color',[.2 .2 .2],'linewidth',0.5,'linestyle','-')
            line([min(xlim) max(xlim)], [0 0], 'color',[0 0 0],'linewidth',0.5,'linestyle','-')
            
            clear my_sem;
        end
        
        % things to do for every axis:
        set(gca,'TickDir','out','XTick',-180:90:180,'XTickLabel',[],'TickLength',1.5*get(gca,'TickLength'));
        
        
        % TCS: here, I added specific XTicks (and ensured XLims are
        % -180/180)
        if vv == 1 && dd== 1
            ylabel(sprintf('Pre-distractor \n %0.01f to %0.01f s',delay_tpt_range(dd,1),delay_tpt_range(dd,2)));
            set(gca,'XTick',-180:90:180,'XTickLabel',{'-180','','0','','180'});
            
        elseif vv == 1 && dd== 2
            ylabel(sprintf('During distractor \n %0.01f to %0.01f s',delay_tpt_range(dd,1),delay_tpt_range(dd,2)));
            set(gca,'XTick',-180:90:180,'XTickLabel',{'-180','','0','','180'});
            
        elseif vv == 1 && dd== 3
            ylabel(sprintf('Post-distractor \n %0.01f to %0.01f s',delay_tpt_range(dd,1),delay_tpt_range(dd,2)));
            if dd == length(delay_tpts)
                xlabel('Position (\circ)');
                set(gca,'XTick',-180:90:180,'XTickLabel',{'-180','','0','','180'});
                
            end
        end
        if dd == 1
            title(ROIs{vv});
        end
        
        if vv ~= 1
            set(gca,'YTickLabel',[]);
        end
        
        xlim([-180 180]);
        hold off;
    end
end
match_ylim(get(gcf,'Children'));
set(gcf,'position', [ 23         245        2386         453])



%% t-test between condition during each epoch against shuffled null
%
% to shuffle:
% - compute n_runs fidelity for each condition for each subj
% - shuffle labels of *these data* 1000 times, then average across
% condition/subj, compute null T distribution; compute p against this
% (two-tailed) (TODO)


% figure out how many total runs... (hack...)
nruns = 0;
for ss = 1:length(subj)
    % get one example set of trials
    ru = unique(all_r(all_subj==ss & all_ROIs==vv));
    
    % figure out how many runs we have for that subj, add to total
    nruns = nruns+length(ru);
    clear ru;
end


% store 'real' T (ROI x epoch)
fidelity_condcompare_realT = nan(length(ROIs),size(delay_tpt_range,1));

% store iter shuffled T's (ROI x epoch x iter)
fidelity_condcompare_shufT = nan(length(ROIs),size(delay_tpt_range,1),iter);

fprintf('Computing null T-dist for condition comparison w/ in epochs\n');

for vv = 1:length(ROIs)
    
    for tt = 1:size(delay_tpt_range,1)
    
        fprintf('%s, Epoch %i\n',ROIs{vv},tt);
        tic;
        thisd = nan(length(subj),2);
        for ss = 1:length(subj)
            for cc = 1:2
                thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==cc;
                thisd(ss,cc) = mean( mean( all_fidelity(thisidx,delay_tpts{tt},1) , 2), 1 ); % mean across tpts (inner), trials (outer)
                clear thisidx;
            end
        end
        
        % compute real T
        [~,~,~,thisstats] = ttest(thisd(:,1),thisd(:,2));
        fidelity_condcompare_realT(vv,tt) = thisstats.tstat; clear thisstats thisd;
        
        % create a data structure we can shuffle: for each subj/ROI/cond,
        % n_runs datapoints
        data_ru   = nan(nruns,1);
        labels_ru = nan(nruns,2); % condition, subj
        
        tmpidx = 1; % where to put this run....
        
        for ss = 1:length(subj)
            ru = unique(all_r(all_subj==ss & all_ROIs==vv));
            for rr = 1:length(ru)
                for cc = 1:2
                    thisidx = all_r==ru(rr) & all_subj==ss & all_ROIs==vv & all_conds(:,1)==cc;
                    data_ru(  tmpidx) = mean( mean( all_fidelity(thisidx,delay_tpts{tt},1) , 2), 1 );
                    labels_ru(tmpidx,:) = [cc ss];
                    tmpidx = tmpidx+1;
                end
            end
        end
        
        
        % iter times..
        for ii = 1:iter
            
            tmpdata = nan(size(data_ru));
            
            % within each subj, shuffle condition labels
            for ss = 1:length(subj)
                % logical indices for this subject, across runs
                thisidx = labels_ru(:,2)==ss;
                
                % list of numerical indices - plug data in here
                subjidx = find(thisidx);
                
                % used to shuffle subjidx
                shufidx = randperm(length(subjidx));
                
                % tmpdata will have data, shuffled wrt condition label, for
                % each subj
                tmpdata(subjidx) = data_ru(subjidx(shufidx)); 
                
                clear thisidx subjidx shufidx;
            end
            
            
            % generate data structure for T-test (average across
            % runs within condition label within each subj)
            thisd = nan(length(subj),2);
            for ss = 1:length(subj)
                for cc = 1:2
                    thisidx = labels_ru(:,1)==cc & labels_ru(:,2)==ss;
                    thisd(ss,cc) = mean( tmpdata(thisidx), 1 ); % mean across *shuffled* data
                    clear thisidx;
                end
            end
            
            % compute shuffled T
            [~,~,~,thisstats] = ttest(thisd(:,1),thisd(:,2));
            fidelity_condcompare_shufT(vv,ii) = thisstats.tstat; clear thisstats
            
            clear thisd tmpdata;
            
        end
        toc;
    end
    
end


% compute p-values (ROIs x epochs)
fidelity_condcompare_pval = 2 * min(mean(fidelity_condcompare_shufT<=fidelity_condcompare_realT,3), mean(fidelity_condcompare_shufT>=fidelity_condcompare_realT,3));
fidelity_condcompare_fdrthresh = fdr(fidelity_condcompare_pval(:),0.05);

%% fidelity line plot per ROI & true 2-way anova stat

% TCS: replace w/ loading colors from file
dist_colors = spDist_condColors;% [0.7100 0.2128 0.4772; 0 0 1;];
cond_colors = dist_colors;

thisd_store =[];

% TCS: should add option for single-subj lines here too...

plot_indiv = 1; % if 1, plot individual subj lines in line plot, otherwise, don't.
hoffset = 0.15; % +/- this much 

fidlines =  figure('Name','Fidelity lines');

h=[];

for vv = 1:length(ROIs)
    thisd = []; %recon data
    
    % TCS: moved from an inner loop to here
    subplot(1,length(ROIs),vv);hold on;

    for dd =1:length(delay_tpts)
        
        
        for cc = 1:length(cu)
            for ss = 1:length(subj)
                
                
                thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==cu(cc);
                thisd(dd,cc,ss) = mean(mean(all_fidelity(thisidx,delay_tpts{dd},1)));
                thisd_store = [thisd_store; mean(mean(all_fidelity(thisidx,delay_tpts{dd},1))) dd vv cc ss];
                
            end
            % TCS: added sqrt here
            my_sem(dd,cc) = std(thisd(dd,cc,:))/sqrt(length(subj));
        end
        
        % if we want, this is a good place to plot lines connecting indiv
        % subj here
        if plot_indiv == 1
            plot(hoffset*[-1;1],squeeze(thisd(1,:,:)),'-','LineWidth',0.25,'Color',[0.4 0.4 0.4]);
        end
        
        
    end
    
    for cc =1:length(cu)
        % TCS: matched these to HRFs (dot makes it slightly easier to see
        % value when overlapping...)
        h(cc) = plot([1 2 3],mean(thisd(:,cc,:),3)','o-','color',cond_colors(cc,:),'linewidth',1,'MarkerSize',3,'MarkerFaceColor',cond_colors(cc,:));
 
        for dd =1:length(delay_tpts)
            plot([dd dd],[mean(thisd(dd,cc,:),3)+1.*my_sem(dd,cc) mean(thisd(dd,cc,:),3)-1.*my_sem(dd,cc)], '-','LineWidth',1,'color',cond_colors(cc,:))
        end
    end
    clear thisidx
    clear thisd
    
    xlim([0.5 3.5])
    if vv ==1
        title(ROIs{vv})
        
        ylabel('WM target Fidelity')
        set(gca,'Xtick',[0 1 2 3 4],'Xticklabel',{'','Epoch 1','Epoch 2','Epoch 3',''},'XTickLabelRotation',45,'TickDir','out');
        
    else
        set(gca,'Xtick',[0 1 2 3 4],'Xticklabel',{'','','','',''},'XTickLabelRotation',45,'TickDir','out');
        title(ROIs{vv})
    end
    
    ylim([0 0.65])
    
    
    
end

y = thisd_store(:,1);
epoch_var = thisd_store(:,2);
roi_var =thisd_store(:,3);
cond_var =thisd_store(:,4);
subj_var = thisd_store(:,5);

h_roi = nan(length(ROIs),3);
p_roi = nan(length(ROIs),3);
t_stats = cell(length(ROIs),3);
p_fdr =nan(length(ROIs),3);
p_masked =nan(length(ROIs),3);


% t-test across epoch per condition
for vv = 1:length(ROIs)
    
    for ee =1:3 %epoch var
        thisidx_nod = roi_var ==vv & epoch_var == ee & cond_var ==1;
        thisidx_dist = roi_var ==vv & epoch_var == ee & cond_var ==2;
        [h_roi(vv,ee), p_roi(vv,ee),~,t_stats{vv,ee}] = ttest(y(thisidx_nod),y(thisidx_dist));
        clear thisidx_nod thisidx_dist
    end
    if vv == length(ROIs) && ee ==3
        [p_fdr, p_masked] = fdr(p_roi(:),0.05); % TCS: added (:) to p_roi
        corr_pval_ind = p_roi(p_masked);
    else
    end
    
end
clear vv ee
% perm ttest
for zz = 1:iter
    
    for vv = 1:length(ROIs)
        
        for ee =1:3 %epoch var
            y_shuf= y;
            
            for sbj = 1:length(subj)
                thisidx = roi_var ==vv & subj_var==sbj & epoch_var ==ee; % dont select by condition
                tmp_y = y(thisidx);
                shuff_idx = randperm(size(tmp_y,1))';
                thisidx_val = find(thisidx);
                y_shuf(thisidx_val,:) = tmp_y(shuff_idx,:);
            end
            
            thisidx_nod = roi_var ==vv & epoch_var == ee & cond_var ==1;
            thisidx_dist = roi_var ==vv & epoch_var == ee & cond_var ==2;
            [h_roi_iter(vv,ee,zz), p_roi_iter(vv,ee,zz),~,t_stats_iter(vv,ee,zz)] = ttest(y_shuf(thisidx_nod),y(thisidx_dist));
            clear thisidx_nod thisidx_dist y_shuf
        end
        
        if vv == length(ROIs) && ee ==3
            [p_fdr_iter(:,:,zz), p_masked_iter(:,:,zz)] = fdr(p_roi_iter(:,:,zz),0.05);
        else
        end
        
    end
    clear vv ee
end

% TCS: why are some struct arrays and others cell arrays?

%% plot results of shuffling (hist)
if plot_hists == 1
    figure('name','ttest perm')
    for vv =1:length(ROIs)
        for ee =1:3
            subplot(1,length(ROIs),vv); hold on;
            histogram([t_stats_iter(vv,ee,:).tstat])
            line([t_stats{vv,ee}.tstat t_stats{vv,ee}.tstat], [min(ylim) max(ylim)], 'color','r','linewidth',2,'linestyle','--')
            exact_p = sum([t_stats_iter(vv,ee,:).tstat] >= t_stats{vv,ee}.tstat)/iter;
            text(max(xlim)-(.5*max(xlim)),max(ylim)-(.25*max(ylim)),sprintf('%0.3f',exact_p),'FontSize',9)
        end
    end
end

%% add significance markers to figure
% this is for parametric stats:
% for vv= 1:length(ROIs)
%     for ee =1:3
%         figure(reconcomp); hold on;
%         subplot(length(delay_tpts),length(ROIs),(ee-1)*length(ROIs)+vv); hold on;
%         
%         if p_roi(vv,ee) > p_fdr && p_roi(vv,ee) <= 0.05
%             text(max(xlim)-(.3*max(xlim)),max(ylim)-(.3*max(ylim)),'*','color',[.5 .5 .5],'fontsize',30)
%         elseif  p_roi(vv,ee) <= p_fdr
%             text(max(xlim)-(.3*max(xlim)),max(ylim)-(.3*max(ylim)),'*','color','k','fontsize',30)
%             
%         end
%     end
% end

% for shuffled stats (TCS)
for vv= 1:length(ROIs)
    for ee =1:size(delay_tpt_range,1)
        
        figure(reconcomp); hold on;
        subplot(length(delay_tpts),length(ROIs),(ee-1)*length(ROIs)+vv); hold on;
        
        if fidelity_condcompare_pval(vv,ee) > fidelity_condcompare_fdrthresh && fidelity_condcompare_pval(vv,ee) <= 0.05
            text(max(xlim)-(.3*max(xlim)),max(ylim)-(.3*max(ylim)),'*','color',[.5 .5 .5],'fontsize',20)
        elseif  fidelity_condcompare_pval(vv,ee) <= fidelity_condcompare_fdrthresh
            text(max(xlim)-(.3*max(xlim)),max(ylim)-(.3*max(ylim)),'*','color','k','fontsize',20)
        end
    end
end

%% shuffled ANOVAs (geh)
P_truth ={};
T_truth = {};
F_store_truth_2=[];

for vv = 1:length(ROIs)
    
    thisroiidx = roi_var ==vv;
    thisy = y(thisroiidx);
    thisepoch = epoch_var(thisroiidx);
    thiscond =cond_var(thisroiidx);
    thissubj=subj_var(thisroiidx);
    [P_truth{vv},T_truth{vv},~] = anovan(thisy,{thissubj,thisepoch,thiscond},'model','full','random',1,'varnames',{'subj','epoch','cond'},'display','off');
    [F_store_truth_2(vv,:)] = RMAOV2_gh([thisy,thisepoch,thiscond,thissubj],0.05);
    
    clear thisroiidx thisepoch thiscond thissubj thisy
    
    
    %[p_fdr(tt,:) p_masked(tt,:)] = fdr(P_roi(vv,:),0.05); dont need to
    %correct? confirm w/ TCS
    
end



match_xlim(get(gcf,'Children'));
legend(h, {'No distractor', 'Distractor'})
set(gcf,'position', [ 23         245        2386         453])

% %% 3-way ANOVA permutation
% % get fidelity , over average delay epochs, for each condition, permutation test
% %%%% this takes ~2.7 hrs minutes to run!!!!! (anovan)
%
% cu =[1 2];
% thisfide_store = [];
% the_y_store= [];
%
% for dd = 1:length(delay_tpts)
%     for vv = 1:length(ROIs)
%         for cc = 1:length(cu)
%             for ss = 1:length(subj)
%                 thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==cu(cc);
%                 thisfide_store(dd,vv,cc,ss) = mean(mean(all_fidelity(thisidx,delay_tpts{dd},1),2));  %targ aligned
%                 the_y_store = [the_y_store;   thisfide_store(dd,vv,cc,ss)  dd vv cc ss];
%             end
%
%         end
%
%
%     end
% end
%
% % sort y_store into individual columns for feasibility/ transparency
% y_3 = the_y_store(:,1);
% epoch_var_3 = the_y_store(:,2);
% roi_var_3 =the_y_store(:,3);
% cond_var_3 =the_y_store(:,4);
% subj_var_3 = the_y_store(:,5);
%
% P_truth_3= cell(1,1);
% T_truth_3= cell(1,1);
% F_store_truth_3 =nan(1,7);
%
% %perform the true ANOVA once
% thisy = y_3;
% thisepoch = epoch_var_3;
% thiscond =cond_var_3;
% thissubj=subj_var_3;
% thisroi =roi_var_3;
% [P_truth_3,T_truth_3,~] = anovan(thisy,{thissubj,thisepoch,thiscond,thisroi},'model','full','random',1,'varnames',{'subj','epoch','cond','roi'},'display','off');
% [F_store_truth_3] = RMAOV33_gh([thisy,thisepoch,thiscond,thisroi,thissubj],0.05);
%
% clear thisroiidx thisepoch thiscond thissubj thisy
%
% %[p_fdr(tt,:) p_masked(tt,:)] = fdr(P_roi(vv,:),0.05); dont need to
% %correct? confirm w/ TCS
%
% % now compute n_iter and store output, compare to true F-values for each
% % factor and each interaction
%
% iter = 1000;
% fprintf(sprintf('computing %i 3-way ANOVAs',iter))
%
% T_iter_3 =cell(iter,1);
% P_iter_3 =cell(iter,1);
% F_store_iter_3 =nan(iter,7);
%
%
% % 3-way anova for each ROI
% tic
%
% for xx=1:iter
%
%       y_shuf=y_3;  %start with empty each iter
%
%         for sbj = 1:length(subj)
%             thisidx =  subj_var==sbj;
%             tmp_y = y(thisidx);
%             shuff_idx = randperm(size(tmp_y,1))';
%             thisidx_val = find(thisidx);
%             y_shuf(thisidx_val,:) = tmp_y(shuff_idx,:);
%             clear thisidx shuff_idx thisidx_val tmp_y
%         end
%
%       thisepoch = epoch_var_3; %maintain og idx for all IVs
%       thiscond =cond_var_3;
%       thissubj=subj_var_3;
%       thisroi =roi_var_3;
%
%    [P_iter_3{xx},T_iter_3{xx},~] = anovan(y_shuf,{thissubj,thisepoch,thiscond,thisroi},'model','full','random',1,'varnames',{'subj','epoch','cond','roi'},'display','off');
%    [F_store_iter_3(xx,:)] = RMAOV33_gh([y_shuf,thisepoch,thiscond,thisroi,thissubj],0.05);
%
%    clear thisidx thisepoch thiscond thissubj thisy tmpy shuffidx thisroi
%
%
% %[p_fdr(tt,:) p_masked(tt,:)] = fdr(P_roi(vv,:),0.05); dont need to
% %correct?
%
% end
%
% toc
%
% %plot F distributions
% %which F vals do we want? idx is as follows {iter,1}{iv_str_3,6} %1 is
% %fixed because we store iters vertically(1 column). 6 is fixed bc this is the Fval col
% iv_str_3 ={'Source';'subj';'epoch';'cond';'roi';'subj*epoch';...
%    'subj*cond';'subj*roi';'epoch*cond'; 'epoch*roi';'cond_roi';...
%    'subj*epoch*cond'; 'subj*epoch*roi'; 'subj*cond*roi';'epoch*cond*roi';'subj*epoch*cond*roi';...
%     'Error';'Total'};
%
%
% col_idx = 6; % F val col
%
% %what factor do we care about? concurs w iv_str_3
% which_effect = [3 4 5 9 10 11 15]; %3 epoch 4 cond 5 roi 9 epoch*cond 10 epoch*cond 11 cond*roi 15 epohc*cond*roi
% exact_store_tmp=[];
%
% exact_store_3=[];
% figure('name','3way perm; anovan')
% for ww =1:length(which_effect)
%    extract_vals=[];
%
% for ii=1:iter
%    extract_vals = [extract_vals; T_iter_3{ii}{which_effect(ww),col_idx}];
% end
% hold on;
%
% subplot(1,length(which_effect),ww)
% histogram(extract_vals)
% line([T_truth_3{which_effect(ww),col_idx} T_truth_3{which_effect(ww),col_idx}], [min(ylim) max(ylim)], 'color','r','linewidth',2,'linestyle','--')
% title(iv_str_3{which_effect(ww)},'Interpreter','none');
% exact_p = sum(extract_vals >= T_truth_3{which_effect(ww),col_idx})/iter;
% text(max(xlim)-(.5*max(xlim)),max(ylim)-(.25*max(ylim)),sprintf('%i',exact_p),'FontSize',9)
% exact_store_3(ww) = [exact_store_tmp; exact_p];
% clear exact_p
%
% if ww ==1
%     xlabel('T-stat')
%     ylabel('Frequency of T-stat')
% else
% end
%
% end
%
% ta_3 = table(exact_store_3(:,1),exact_store_3(:,2),exact_store_3(:,3),exact_store_3(:,4),exact_store_3(:,5),exact_store_3(:,6),exact_store_3(:,7));
% ta_3.Properties.VariableNames={'Epoch','Cond','ROI','EpochCond','EpochROI','CondROI','EpochCondROI'};
%
% % make analagous plot for RMAOV33 to compare perms
%
% %what factor do we care about? concurs w iv_str_3
% which_effect_rm_3 = [1 2 3 4 5 6 7]; %3 epoch 4 cond 5 roi 9 epoch*cond 10 epoch*cond 11 cond*roi 15 epohc*cond*roi
% rm_str ={'epoch','cond','roi','epochcond','epochroi','condroi','epochcondroi'};
% exact_store_tmp_rm=[];
% exact_store_rm3=[];
%
% figure('name','3way perm; RMOV33')
% for ww =1:length(which_effect_rm_3)
%    extract_vals=[];
%
% for ii=1:iter
%    extract_vals = [extract_vals; F_store_iter_3(ii,which_effect_rm_3(ww))];
% end
% hold on;
%
% subplot(1,length(which_effect_rm_3),ww)
% histogram(extract_vals)
% %
% line([F_store_truth_3(which_effect_rm_3(ww)) F_store_truth_3(which_effect_rm_3(ww))], [min(ylim) max(ylim)], 'color','r','linewidth',2,'linestyle','--')
% title(rm_str{which_effect_rm_3(ww)},'Interpreter','none');
% exact_p = sum(extract_vals >= F_store_truth_3(which_effect_rm_3(ww)))/iter;
% text(max(xlim)-(.5*max(xlim)),max(ylim)-(.25*max(ylim)),sprintf('p=%0.4f',exact_p),'FontSize',9)
% exact_store_rm3(ww) = [exact_store_tmp_rm; exact_p];
% clear exact_p
%
% if ww ==1
%     xlabel('T-stat')
%     ylabel('Frequency of T-stat')
% else
% end
%
% end
%
% %%%% TMP COM
%% 2-way ANOVA permutation
% get fidelity , over average delay epochs, for each condition, permutation test
%%%% this takes ~24 minutes to run!!!!!

cu =[1 2];
thisfide_store = [];
the_y_store= [];

for dd = 1:length(delay_tpts)
    for vv = 1:length(ROIs)
        for cc = 1:length(cu)
            for ss = 1:length(subj)
                thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==cu(cc);
                thisfide_store(dd,vv,cc,ss) = mean(mean(all_fidelity(thisidx,delay_tpts{dd},1),2));  %targ aligned
                the_y_store = [the_y_store;   thisfide_store(dd,vv,cc,ss)  dd vv cc ss];
            end
            
        end
        
        
    end
end


% sort y_store into individual columns for transparency
y = the_y_store(:,1);
epoch_var = the_y_store(:,2);
roi_var =the_y_store(:,3);
cond_var =the_y_store(:,4);
subj_var = the_y_store(:,5);


fprintf(sprintf('computing %i 2-way ANOVAs per ROI',iter))

T_iter =cell(length(ROIs),iter);
P_iter =cell(length(ROIs),iter);
F_store_iter_2 = nan(length(ROIs),iter,3);
% two-way anova for each ROI
tic

for xx=1:iter
    for vv = 1:length(ROIs)
        y_shuf =y;
        thisroiidx = roi_var ==vv;
        for sbj = 1:length(subj)
            thisidx = roi_var ==vv & subj_var==sbj;
            tmp_y = y(thisidx);
            shuff_idx = randperm(size(tmp_y,1))';
            thisidx_val = find(thisidx);
            y_shuf(thisidx_val,:) = tmp_y(shuff_idx,:);
        end
        thisepoch = epoch_var(thisroiidx);
        thiscond =cond_var(thisroiidx);
        thissubj=subj_var(thisroiidx);
        
        [P_iter{vv,xx},T_iter{vv,xx},~] = anovan(y_shuf(thisroiidx),{thissubj,thisepoch,thiscond},'model','full','random',1,'varnames',{'subj','epoch','cond'},'display','off');
        [F_store_iter_2(vv,xx,:)] = RMAOV2_gh([y_shuf(thisroiidx),thisepoch,thiscond,thissubj],0.05);
        
        clear thisroiidx thisepoch thiscond thissubj thisy tmpy shuffidx
        
        
        %[p_fdr(tt,:) p_masked(tt,:)] = fdr(P_roi(vv,:),0.05); dont need to
        %correct?
        
    end
end
toc

%plot F distributions; first, for anovan outputs
%which F vals do we want? idx is as follows {iter,1}{iv_str,6} %1 is
%fixed because we store iters vertically(1 column). 6 is fixed bc this is
%the Fval col. iv_str is defined from the vertical organization of the
%output
iv_str ={'Source';'subj';'epoch';'cond';'subj*epoch';...
    'subj*cond';'epoch*cond';...
    'subj*epoch*cond';...
    'Error';'Total'};


col_idx = 6; % F-val col

%what factor do we care about? concurs w iv_str(vs)
which_effect = [3 4 7]; %3 epoch 4 cond 7epoch*cond
exact_store_tmp=[];

for vv =1:length(ROIs)
    figure('name','2-way perm;anovan')
    for ww =1:length(which_effect)
        extract_vals=[];
        for ii=1:iter
            extract_vals = [extract_vals; T_iter{vv,ii}{which_effect(ww),col_idx}];
        end
        hold on;
        
        subplot(1,length(which_effect),ww)
        histogram(extract_vals)
        line([T_truth{vv}{which_effect(ww),col_idx} T_truth{vv}{which_effect(ww),col_idx}], [min(ylim) max(ylim)], 'color','r','linewidth',2,'linestyle','--')
        title(iv_str{which_effect(ww)},'Interpreter','none');
        exact_p = sum(extract_vals >= T_truth{vv}{which_effect(ww),col_idx})/iter;
        text(max(xlim)-(.5*max(xlim)),max(ylim)-(.25*max(ylim)),sprintf('%i',exact_p),'FontSize',9)
        exact_store(vv,ww) = [exact_store_tmp; exact_p];
        clear exact_p
        
        if ww ==1
            xlabel('T-stat')
            ylabel('Frequency of T-stat')
        else
        end
        
    end
end
ta = table(ROIs',exact_store(:,1),exact_store(:,2),exact_store(:,3));
ta.Properties.VariableNames={'ROIs','Epoch','Cond','EpochCond'};

% plot for RMAOV2 output
iv_str_rm ={'epoch','cond','epoch*cond'};
which_effect_rm = [1 2 3];
exact_store_tmp=[];
exact_store_rm_2 =nan(length(ROIs),length(which_effect_rm));

for vv =1:length(ROIs)
    figure('name','2-way perm; RMAOV2')
    
    for ww =1:length(which_effect_rm)
        extract_vals=[];
        for ii=1:iter
            extract_vals = [extract_vals; F_store_iter_2(vv,ii,which_effect_rm(ww))];
        end
        hold on;
        
        subplot(1,length(which_effect_rm),ww)
        histogram(extract_vals)
        line([F_store_truth_2(vv,which_effect_rm(ww)) F_store_truth_2(vv,which_effect_rm(ww))], [min(ylim) max(ylim)], 'color','r','linewidth',2,'linestyle','--')
        title(iv_str_rm{which_effect_rm(ww)},'Interpreter','none');
        exact_p = sum(extract_vals >= F_store_truth_2(vv,which_effect_rm(ww)))/iter;
        text(max(xlim)-(.5*max(xlim)),max(ylim)-(.25*max(ylim)),sprintf('p=%0.4f',exact_p),'FontSize',9)
        exact_store_rm_2(vv,ww) = [exact_store_tmp; exact_p];
        clear exact_p
        
        if ww ==1
            xlabel('T-stat')
            ylabel('Frequency of T-stat')
        else
        end
        
    end
end
%% plot sigs on the subplots -- anovan
sig_colors = lines(3);

%correct them pvals
p_fdr_perm= nan(3);
p_masked_perm= nan(length(ROIs),3);

for ee =1:3 %epoch var
    [p_fdr_perm(ee), p_masked_perm(:,ee)] = fdr(exact_store(:,ee),0.05);
    
    
end

sig_mrkr ={'o','+','x'};
y_mod = [.1 .15 .2];
for vv = 1:length(ROIs)
    for ee =1:3
        figure(fidlines); hold on;
        subplot(1,length(ROIs),vv)
        
        
        if exact_store(vv,ee) > p_fdr_perm(ee) && exact_store(vv,ee) <= 0.05
            text(max(xlim)-(.2*max(xlim)),max(ylim)-(y_mod(ee)*max(ylim)),sprintf('%s',sig_mrkr{ee}),'color',[.5 .5 .5],'fontsize',15)
        elseif  exact_store(vv,ee) <= p_fdr_perm(ee)
            text(max(xlim)-(.2*max(xlim)),max(ylim)-(y_mod(ee)*max(ylim)),sprintf('%s',sig_mrkr{ee}),'color','k','fontsize',15)
        else
        end
        
        
    end
    
end

% %perform again for rmaov2
% for vv = 1:length(ROIs)
%     for ee =1:3
%         figure(fidlines); hold on;
%         subplot(1,length(ROIs),vv)
%         if exact_store_rm_2(vv,ee) >0.01 && exact_store_rm_2(vv,ee) <= 0.05
%             text(max(xlim)-(.7*max(xlim)),max(ylim)-(y_mod(ee)*max(ylim)),'*','color',sig_colors(ee,:))
%         elseif exact_store_rm_2(vv,ee) >0.001 && exact_store_rm_2(vv,ee) <= 0.01
%             text(max(xlim)-(.7*max(xlim)),max(ylim)-(y_mod(ee)*max(ylim)),'**','color',sig_colors(ee,:))
%         elseif exact_store_rm_2(vv,ee) <= 0.001
%             text(max(xlim)-(.7*max(xlim)),max(ylim)-(y_mod(ee)*max(ylim)),'***','color',sig_colors(ee,:))
%
%
%         end
%     end
%
% end

fprintf('2-way computations are done. ... onto 3-way')
% tmp

% %% plot sigs on the subplots - 3-way  ... TODO: i dont think this is necessary / would not make sense on the figure 4 ROI plot
% % sig_colors = lines(3);
%  y_mod =[.1 .15 .2];
%
% for vv = 1:length(ROIs)
%     for ee =1:3
%         figure(linep); hold on;
%         subplot(1,length(ROIs),vv)
%         if exact_store(vv,ee) >0.01 && exact_store(vv,ee) <= 0.05
%             text(max(xlim)-(.7*max(xlim)),max(ylim)-(y_mod(ee)*max(ylim)),'*','color',sig_colors(ee,:))
%         elseif exact_store(vv,ee) >0.001 && exact_store(vv,ee) <= 0.01
%             text(max(xlim)-(.7*max(xlim)),max(ylim)-(y_mod(ee)*max(ylim)),'**','color',sig_colors(ee,:))
%         elseif exact_store(vv,ee) <= 0.001
%             text(max(xlim)-(.7*max(xlim)),max(ylim)-(y_mod(ee)*max(ylim)),'***','color',sig_colors(ee,:))
%
%
%         end
%     end
%
% end

return
