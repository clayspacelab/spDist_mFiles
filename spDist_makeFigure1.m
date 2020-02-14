
%spDist_plotEyeData.m
% in the actual task, cyan = distractor, magenta = no distractor 
% dependencies: misc_util, RMAOV1 

root = spDist_loadRoot;

%load raw subject data 
subj = {'KD','CC','AY','MR','XL','EK','SF'};
sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist2'},{'spDist1','spDist2'},{'spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}}; %two sessions removed


WHICH_EXCL = [13 20 22]; % see geh spDist_eyeDataNotes.txt on how/why these exclu criteria were chosen. 
if ismember(WHICH_EXCL,13)
    which_excl_str ={'broken fix'};
elseif ismember(WHICH_EXCL,[13 20])
    which_excl_str ={'broken fix','no sacc'};
elseif ismember(WHICH_EXCL, [13 20 21])
    which_excl_str ={'broken fix','no sacc','i sacc too small/short'};
elseif ismember(WHICH_EXCL, [13 20 22])
    which_excl_str ={'broken fix','no sacc','i sacc err too lg'};
elseif ismember(WHICH_EXCL, [13 20 21 22])
    which_excl_str ={'broken fix','no sacc','i sacc too small/short','i sacc err too lg'};
else
    error('which exclusion criteria have you chosen?')
end

% first-digit:
% - 1 - trial-level exclusion (bad drift correction [11], calibration [12], or delay-
%       fixation break [13]
% - 2 - primary saccade exclusion (no primary sacc [20]; too small/short [21] bi, large error [22]ei)

%21 bad initial saccade (duration/amplitude outside range)
% 22 iniital saccade error


% concatenate ALL subject data
all_subj = nan(1000*length(subj),1);

all_data = [];

startidx = 1;

for ss = 1:length(subj)
    for sessidx = 1:length(sess{ss})
        
        
        fn = sprintf('%s/spDist_behav_82719/%s_%s_scored.mat',root,subj{ss},sess{ss}{sessidx});
        fprintf('Loading scored eye data from %s\n',fn);
        this_scored = load(fn);
        
        this_data.s_all = this_scored.ii_sess;
        this_data.sess_all = sessidx;
        
        this_subj = ss;
        
        all_data = cat_struct(all_data,this_data);
        all_subj(startidx:(startidx-1+size(this_scored.ii_sess.trialinfo,1))) = this_subj;
        
        startidx = startidx+size(this_scored.ii_sess.trialinfo,1);
        
        clear this_subj this_data;
    end
end


all_subj = all_subj(1:(startidx-1));
all_data.subj_all = all_subj;

% determine which trials to include
% first, narrow based on saccade preprocessing/scoring exclusions
all_data.use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, WHICH_EXCL), all_data.s_all.excl_trial, 'UniformOutput',false));
all_data.use_trial(all_data.s_all.f_sacc_err>10) = 0; %exclude trials with errors > 10 deg
all_data.use_trial(all_data.s_all.i_sacc_err>10) = 0;
all_data.use_trial(ismember(all_data.s_all.r_num, [8,12,14]) & all_subj==3) = 0; %3 here specifically refers to subj 3 above (AY), see spDist_eyeDataNotes for 
% drop trials with very short (< 100 ms) or very long RT (> 1 s)
all_data.use_trial(all_data.s_all.i_sacc_rt<0.1 | all_data.s_all.i_sacc_rt>1.0) = 0;

%% organize data
distractor_bins = unique(all_data.s_all.trialinfo(all_data.s_all.trialinfo(:,1)~=1,6));
distractor_spacing = 360/length(distractor_bins);

flip_bins = [1 2 3]; %which distractor bins should we "flip" the Y in order to collapse across CW/CCW distractors?

%initialize storage cells
all_err = cell(3,1);
all_mu = cell(3,1);
all_rt = cell(3,1);

cond_str = {'No distractor','Distractor'};
cond_colors = lines(3);

params_of_interest = {'f_sacc'};
param_str = {};


% first, error for no-distractor trials
tmp_err = nan(1,length(params_of_interest),2,length(subj)); % initial, final
tmp_mu  = nan(1,length(params_of_interest),2,length(subj)); % initial, final

for ss = 1:length(subj)
    thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==1 & all_data.use_trial==1;
    for pp = 1:length(params_of_interest)
        % distractor bin x param x [radial; tangential] x subj
        tmp_err(1,pp,:,ss) = nanstd( all_data.s_all.(params_of_interest{pp})(thisidx,:), [], 1 );
        tmp_mu(1,pp,:,ss)  = nanmean( all_data.s_all.(params_of_interest{pp})(thisidx,:), 1 );
    end
end

all_err{1} = tmp_err; % 1st cell of all err and all mu = no distractor
all_mu{1} = tmp_mu; %note: the first index of all_ is NO DIST

clear tmp_err tmp_mu;

tmp_err = nan(1,length(params_of_interest),2,length(subj)); 
tmp_mu = nan(1,length(params_of_interest),2,length(subj));
tmp_cw = nan(1,length(params_of_interest),2,length(subj));
tmp_ccw = nan(1,length(params_of_interest),2,length(subj));

for ss = 1:length(subj)
    
    for bb =1 %placeholder, not relevant here bc we're collecting bin zero, but want to keep dimens same size as when we collect others 
       
        for pp = 1:length(params_of_interest)
            tmpidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & all_data.s_all.trialinfo(:,6)==0 & all_data.s_all.trialinfo(:,10) < 0; %negative CW jitter
            orig_y =  all_data.s_all.(params_of_interest{pp})(tmpidx,2);
            orig_y_flip = orig_y*-1;
            all_data.s_all.(params_of_interest{pp})(tmpidx,2) = orig_y_flip;  % here is where the actual y flip is inserted for CW bins.
            
            thisorigidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & all_data.s_all.trialinfo(:,6)==0; %now, collect all jitters.
 
            % distractor bin x param x [radial; tangential] x subj
            tmp_err(bb,pp,:,ss) = nanstd( all_data.s_all.(params_of_interest{pp})(thisorigidx,:), [], 1 );
            tmp_mu(bb,pp,:,ss) = nanmean( all_data.s_all.(params_of_interest{pp})(thisorigidx,:),  1 );
           
        end
    end
end

all_err{2} = tmp_err;
all_mu{2} = tmp_mu;
clear tmp_err tmp_mu;

%put other bins in
flip_bins =[1 2 3];
tmp_err = nan(length(flip_bins),length(params_of_interest),2,length(subj));
tmp_mu = nan(length(flip_bins),length(params_of_interest),2,length(subj));

for ss = 1:length(subj)
    for bb = 1:length(flip_bins)        
        for pp = 1:length(params_of_interest)
            tmpidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & all_data.s_all.trialinfo(:,6)==(flip_bins(bb)*-1);
            orig_y =  all_data.s_all.(params_of_interest{pp})(tmpidx,2);
            orig_y_flip = orig_y*-1;
            all_data.s_all.(params_of_interest{pp})(tmpidx,2) = orig_y_flip;  % here is where the actual y flip is inserted for CW bins.
            
            thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & (all_data.s_all.trialinfo(:,6)==flip_bins(bb) | all_data.s_all.trialinfo(:,6)==flip_bins(bb)*-1);
            % distractor bin x param x [radial; tangential] x subj
            tmp_err(bb,pp,:,ss) = nanstd( all_data.s_all.(params_of_interest{pp})(thisidx,:), [], 1 );
            tmp_mu(bb,pp,:,ss) = nanmean( all_data.s_all.(params_of_interest{pp})(thisidx,:),  1 );
            
        end
    end
end

all_err{3} = tmp_err;
all_mu{3} = tmp_mu;
clear tmp_err tmp_mu;

% first, error for no-distractor trials
tmp_rt = nan(1,length(params_of_interest),1,length(subj)); % initial, final
params_of_interest = {'i_sacc_rt'};
param_str = {'RT'};

for ss = 1:length(subj)
    
    thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==1 & all_data.use_trial==1;

    for pp = 1:length(params_of_interest)
        % distractor bin x param x rt x subj
        tmp_rt(1,pp,1,ss)  = nanmean( all_data.s_all.(params_of_interest{pp})(thisidx),1);

    end
    
end

all_rt{1} = tmp_rt;
clear tmp_rt
tmp_rt = nan(1,length(params_of_interest),1,length(subj)); % initial, final

for ss = 1:length(subj)
    
    thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & all_data.s_all.trialinfo(:,6)==0; % cw

    for pp = 1:length(params_of_interest)
        % distractor bin x param x rt x subj
        tmp_rt(1,pp,1,ss)  = nanmean( all_data.s_all.(params_of_interest{pp})(thisidx),1); %cw
    
    end
    
end

all_rt{2} = tmp_rt;
flip_bins =[1 2 3];
tmp_rt = nan(length(flip_bins),length(params_of_interest),1,length(subj)); % initial, final
for ss = 1:length(subj)
    for bb = 1:length(flip_bins)
        thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & (all_data.s_all.trialinfo(:,6)==flip_bins(bb) | all_data.s_all.trialinfo(:,6)==flip_bins(bb)*-1);
        
        for pp = 1:length(params_of_interest)
            % distractor bin x param x [radial; tangential] x subj
            tmp_rt(bb,pp,:,ss) = nanmean( all_data.s_all.(params_of_interest{pp})(thisidx,:),  1 );
            
        end
    end
end

all_rt{3} = tmp_rt;
%% Figure 1B : Example participant eye-trace 

dist_colors = [0.7100 0.2128 0.4772; 0 0 1;]; %1 is red, no distractor, 2 is blue (near distractor), 3 is green (far distractor)
figure;
hold on;
xdat_to_plot = [1 2 3 4 5 6 7 8 9 10]; % these are the annotations from our eye-date marking the distinct epochs of the task
dist_cond =[1 2];
%use subj 6, trial # 23 no d and #12 d.
for cc = 1:length(dist_cond)
    thisidx = find(all_data.use_trial==1 & all_data.subj_all ==6 & all_data.s_all.trialinfo(:,1) == dist_cond(cc));
    if cc ==1
        which_tri = 23;
    else
        which_tri = 12;
    end
    for tt = which_tri
        figure(1)
        %transform the cartesian X,Y coords we have for the trace and put
        %it in polar
        [tmpth,tmpr] = cart2pol(all_data.s_all.X{thisidx(tt)},all_data.s_all.Y{thisidx(tt)});
        
        % change th, keeping r the same, based on th of all_data.s_all.targ
        [adjth,~] = cart2pol(all_data.s_all.targ(thisidx(tt),1),all_data.s_all.targ(thisidx(tt),2));
        
        [aligned_x,~] = pol2cart(tmpth-adjth,tmpr); % we're not using y here, just x
        
        aligned_x = aligned_x(ismember(all_data.s_all.XDAT{thisidx(tt)},xdat_to_plot));
        if cc ==2
            aligned_x = aligned_x.*-1; % for visualization purposes, flip the trace of the distractor condition to be mirrored over the x-axis
        else
        end
        
        this_t = (1:length(aligned_x))/500; %Time in seconds (how many samples were taken / frequency of recording device)
        subplot(2,1,1)
        plot(this_t,aligned_x,'-','LineWidth',2,'Color',dist_colors(cc,:));
        ylim([-18 18]);
        xlim([0 15])
        xticks([0 1.5 6 13.5])
        xticklabels({'-1.5','0','4.5','12'})
        set(gca,'YTick',[-18:3:18],'TickDir','out');
        hold on;
        
        
        if cc ==1
            %for this subject, this trial, this epoch, turn XDAT marker on
            epoch_1_condc = all_data.s_all.XDAT{thisidx(tt)} ==1 ; % store the entire epoch in a variable 
            epoch_2_targ = all_data.s_all.XDAT{thisidx(tt)} ==2;
            epoch_3_d1 = all_data.s_all.XDAT{thisidx(tt)}==3;
            epoch_4_dist = all_data.s_all.XDAT{thisidx(tt)} ==4;
            epoch_5_d2 = all_data.s_all.XDAT{thisidx(tt)} ==5;
            epoch_6_go = all_data.s_all.XDAT{thisidx(tt)}==6;
            epoch_7_fb = all_data.s_all.XDAT{thisidx(tt)}==7;
            epoch_8_iti = all_data.s_all.XDAT{thisidx(tt)}==8;
            subplot(2,1,2)
            hold on;
            plot(this_t,epoch_1_condc+6,'k-','LineWidth',1);
            plot(this_t,epoch_2_targ+4,'k-','LineWidth',1);
            plot(this_t,epoch_4_dist+2,'k-','LineWidth',1);
            plot(this_t,epoch_6_go,'k-','LineWidth',1);
            xlim([0 15])
            xticks([0 1.5 6 13.5])
            xticklabels({'-1.5','0','4.5','12'})
            yticks([1 3 5 7])
            yticklabels({'Response','Distractor', 'Target','Condition Cue'})
            xlabel('Time relative to delay onset (s)');
            ylabel('Trial Epoch')
        else
        end     
    end

end

xlabel('Time relative to delay onset (s)');
ylabel('Eye position (DVA)');

%% Figure 1C : Aligned Saccadic Endpoints 

dist_colors = [0 0 1; 0.7100 0.2128 0.4772; .3 .6 .1 ; .3 .6 .1 ; .3 .6 .1 ];

params_of_interest = {'f_sacc'};
param_str = {'f saccade'};

figure;
to_plot = {'f_sacc'}; % what fields do we want to plot?
dist_bins = [-1 0 1 2 3]; %why is this like this? -1 means NO DISTRACTOR

alph = [.4 .4 .4 .4 .4]; %translucence 

for pp = 1:length(to_plot)
    
    for dd =1:length(dist_bins)
        
        if dist_bins(dd) == -1
            thisidx = all_data.s_all.trialinfo(:,1)==1 & all_data.use_trial==1;
            x = all_data.s_all.(to_plot{pp})(thisidx,1);
            y = all_data.s_all.(to_plot{pp})(thisidx,2);
            y_save = nanmean(y);
            
        elseif dist_bins(dd) == 0
            
            thisidx = all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & all_data.s_all.trialinfo(:,6)==0;
            x = all_data.s_all.(to_plot{pp})(thisidx,1);
            y = all_data.s_all.(to_plot{pp})(thisidx,2);
            y_save = nanmean(y);
        else
            thisidx =  all_data.use_trial==1 & (all_data.s_all.trialinfo(:,6) ==dist_bins(dd) | all_data.s_all.trialinfo(:,6) ==dist_bins(dd)*-1); %collect the numbered distractor bins, [1 2 3] = CW, [-1 -2 -3] = CCW. we have alreayd flipped the y-val, we're just collect the bins here
            
            x = all_data.s_all.(to_plot{pp})(thisidx,1);
            y = all_data.s_all.(to_plot{pp})(thisidx,2);
            y_save = nanmean(y);
                 
        end
        
        if dd ==1
            subplot(1,3,1)
            title('No Distractor')
        elseif dd==2
            subplot(1,3,2)
             title('Near Distractor')
        else
            subplot(1,3,3)
             title('Far Distractor')
        end
        hold on;
        scatter(x,y,30,dist_colors(dd,:),'filled','MarkerFaceAlpha',.5)
        scatter(mean(x),mean(y),30,'k','filled')
        hold on;
        plot([3 15], [0 0],'--','linewidth',1,'color',[.1 .1 .1])
        clear tmpy
        clear thisidx
        ylim([-6 6])
        xlim([3 15])
        xticks([5 10 15])
        yticks([-5 0 5])
        xticklabels([5 10 15])
        axis equal
    end
end

set(gcf,'Renderer','painters')
match_ylim(get(gcf,'Children'))

%% Figure 1D : Precision
% directly compare error for no distractor, near distractor, far distractor
% trials (subplot for each param); averaged over radial/tang...

params_of_interest = {'f_sacc'}; %we're using the final saccade
param_str = {'final sacc'};
cond_colors = [0.7100    0.2128    0.4772;0 0 1;0 0 1;];
cond_str = {'No distractor','Near Distractor','Far Distractor'};

figure;

for pp = 1:length(params_of_interest)
    
    subplot(1,length(params_of_interest),pp); hold on;
    
    % distractor cond x subj
    thise = nan(length(all_err),length(subj)); %collect this error
    for ii = 1:length(all_err)
        thise(ii,:) = mean(mean(all_err{ii}(:,pp,:,:),3),1); % mean over radial/tangential; distractor bin; for n_bins, then mean over these bins
    end
    
    plot(1:size(thise,1),thise,'-','Color',[0.3 0.3 0.3],'LineWidth',.5);
    plot(1:size(thise,1),thise,'-','Color',[0.5 0.5 0.5],'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1);
    for ii = 1:length(all_err)
        hold on;
        tmpe = std(thise(ii,:))/sqrt(length(subj));
        plot(ii*[1 1],mean(thise(ii,:))+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(ii,:));
        plot(ii,mean(thise(ii,:)),'o','LineWidth',1.5,'Color',cond_colors(ii,:),'MarkerSize',10,'MarkerFaceColor',cond_colors(ii,:));
    end
    
    xlim([0 length(all_err)+1]);
    
    set(gca,'XTick',1:length(all_err),'XTickLabel',cond_str','XTickLabelRotation',45,'TickDir','out');
    xlabel('Condition');
    if pp == 1
        ylabel('Precision (avg sd, \circ)');
    end
    title(param_str{pp});
    
    the_y = [thise(1,:)';thise(2,:)'; thise(3,:)';];
    the_iv =[ones(length(thise(1,:)'),1); 2*ones(length(thise(1,:)'),1); 3*ones(length(thise(1,:)'),1)];
    subj = [1 2 3 4 5 6 7]';
    the_subj =[subj;subj;subj];
    x = [the_y the_iv the_subj];
    RMAOV1(x,0.05) %one-way RM ANOVA
    

end

match_ylim(get(gcf,'Children'));
%% Figure 1E : RT

cond_colors = [0.7100    0.2128    0.4772;0 0 1;0 0 1;0 0 1;0 0 1;]; %colors for illustrator 
params_of_interest = {'i_sacc_rt'};
param_str = {'RT'};

figure;

for pp = 1:length(params_of_interest)
    
    subplot(1,length(params_of_interest),pp); hold on;                                                                                    
    
    % nbins x nsubj
    thisnod_rt = squeeze(all_rt{1}(1,pp,1,:))';
    thisnear_dist = squeeze(all_rt{2}(:,pp,1,:))';
    thisfar_dist = mean(squeeze(all_rt{3}(:,pp,1,:)));
    thisrt= [thisnod_rt; thisnear_dist; thisfar_dist];
    plot(1:length(flip_bins),thisrt,'-','color',[.2 .2 .2]);
    
    for ii = 1:size(thisrt,1)
        tmpe = std(thisrt(ii,:))/sqrt(length(subj));
        plot(ii*[1 1],mean(thisrt(ii,:))+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(ii,:));
        plot(ii,mean(thisrt(ii,:)),'o','LineWidth',1.5,'color',cond_colors(ii,:),'markerfacecolor',cond_colors(ii,:),'MarkerSize',10);
    end
    title(param_str{pp});
    xlabel('Distractor Condition');
   
    if pp == 1
        ylabel('RT');
        
    end
    set(gca,'XTick',[1 2 3],'Xticklabels',({'No Distractor', 'Near Distractor','Far Distractor'}), 'TickDir','out','XTickLabelRotation',45)
    the_RT = [thisrt(1,:)'; thisrt(2,:)'; thisrt(3,:)']; %use for just dist
    the_ivRT =[ones(length(thisrt(1,:)'),1); 2*ones(length(thisrt(2,:)'),1); 3*ones(length(thisrt(3,:)'),1);]; %use for just dist
    the_subj =[subj;subj;subj;];
    x = [the_RT the_ivRT the_subj];
    RMAOV1(x,0.05)
end
match_ylim(get(gcf,'Children'));
xlim([0 size(thisrt,1)+1]);


%% Figure 1F : Bias
params_of_interest = {'f_sacc'};
param_str = {'final sacc'};

figure;

for pp = 1:length(params_of_interest)
    
    subplot(1,length(params_of_interest),pp); hold on;
    
    % nbins x nsubj
    thism_near = squeeze(all_mu{2}(:,pp,2,:))'; %this is the only condition that is information for this analysis.
    plot(1,thism_near,'o','MarkerSize',5,'Markerfacecolor',[0.3 0.3 0.3], 'Color',[0.3 0.3 0.3]);
    plot(1,mean(thism_near,2),'o','MarkerSize',20,'Color',cond_colors(2,:),'MarkerFaceColor',cond_colors(2,:));
    hold on;
    errorbar(1, mean(thism_near), [std(thism_near)/sqrt(length(subj))],'linewidth',2,'Color',cond_colors(2,:))
    title(param_str{pp});
    
   
    ylabel('Memory Bias, toward distractor (\circ)');
end
xlabel('Near Distractor');

[h_bias_near p_bias_near,CI,STATS] = ttest(thism_near(1,:)') %t-test to determine if bias is different than zero

match_ylim(get(gcf,'Children'));