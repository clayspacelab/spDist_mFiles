
%spDist_plotEyeData.m
% in the actual task, cyan = distractor, magenta = no distractor 
% dependencies: misc_util, RMAOV1 

%root = spDist_loadRoot;
root = '/share/data/spDist/';

%load raw subject data 
subj = {'KD','CC','AY','MR','XL','EK','SF'};
%sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist2'},{'spDist1','spDist2'},{'spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}}; %two sessions removed

sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}}; %two sessions removed

WHICH_EXCL = [13 20 21 22]; % see geh spDist_eyeDataNotes.txt on how/why these exclu criteria were chosen. 
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
        

       % fn = sprintf('%s/spDist_behav_82719/%s_%s_scored.mat',root,subj{ss},sess{ss}{sessidx});
       fn = sprintf('%s/spDist_behav_92220/%s_%s_scored.mat',root,subj{ss},sess{ss}{sessidx});
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
%drop trials with very short (< 100 ms) or very long RT (> 1 s)
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
tmp_errpol = nan(1,length(params_of_interest),2,length(subj));
tmp_mupol = nan(1,length(params_of_interest),2,length(subj));

for ss = 1:length(subj)
    thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==1 & all_data.use_trial==1;
    for pp = 1:length(params_of_interest)
        % distractor bin x param x [radial; tangential] x subj
        tmp_err(1,pp,:,ss) = nanstd( all_data.s_all.(params_of_interest{pp})(thisidx,:), [], 1 );
        tmp_mu(1,pp,:,ss)  = nanmean( all_data.s_all.(params_of_interest{pp})(thisidx,:), 1 );
        % new
        [dtheta_err, drad_err]= cart2pol(all_data.s_all.(params_of_interest{pp})(thisidx,1), all_data.s_all.(params_of_interest{pp})(thisidx,2));
        tmp_errpol(1,pp,:,ss) = nanstd([drad_err rad2deg(dtheta_err)]);
        [dtheta_mu, drad_mu] = cart2pol(all_data.s_all.(params_of_interest{pp})(thisidx,1), all_data.s_all.(params_of_interest{pp})(thisidx,2));
        tmp_mupol(1,pp,:,ss) = nanmean([drad_mu rad2deg(dtheta_mu)]); %put rad as x FIRST, theta as Y second, in keeping with later indexing
    end
end

all_err{1} = tmp_err; % 1st cell of all err and all mu = no distractor
all_mu{1} = tmp_mu; %note: the first index of all_ is NO DIST
all_mupol{1} = tmp_mupol;
all_errpol{1} = tmp_errpol;

clear tmp_err tmp_mu tmp_errpol tmp_mupol;

tmp_err = nan(1,length(params_of_interest),2,length(subj)); 
tmp_mu = nan(1,length(params_of_interest),2,length(subj));
tmp_cw = nan(1,length(params_of_interest),2,length(subj));
tmp_ccw = nan(1,length(params_of_interest),2,length(subj));
tmp_errpol = nan(1,length(params_of_interest),2,length(subj));
tmp_mupol = nan(1,length(params_of_interest),2,length(subj));

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

            [dtheta_err, drad_err]= cart2pol(all_data.s_all.(params_of_interest{pp})(thisorigidx,1), all_data.s_all.(params_of_interest{pp})(thisorigidx,2));
            tmp_errpol(bb,pp,:,ss) = nanstd([drad_err rad2deg(dtheta_err)]);
            [dtheta_mu, drad_mu] = cart2pol(all_data.s_all.(params_of_interest{pp})(thisorigidx,1), all_data.s_all.(params_of_interest{pp})(thisorigidx,2));
            tmp_mupol(bb,pp,:,ss) = nanmean([drad_mu rad2deg(dtheta_mu)]); %put rad as x FIRST, theta as Y second, in keeping with later indexing 

            clear dtheta_err drad_err drad_err dtheta_err
        end
    end
end

all_err{2} = tmp_err;
all_mu{2} = tmp_mu;
all_mupol{2} = tmp_mupol;
all_errpol{2} = tmp_errpol;

clear tmp_err tmp_mu tmp_mupol tmp_errpol;

%put other bins in
flip_bins =[1 2 3];
tmp_err = nan(length(flip_bins),length(params_of_interest),2,length(subj));
tmp_mu = nan(length(flip_bins),length(params_of_interest),2,length(subj));
tmp_errpol = nan(1,length(params_of_interest),2,length(subj));
tmp_mupol = nan(1,length(params_of_interest),2,length(subj));

for ss = 1:length(subj)
    for bb = 1:length(flip_bins)        
        for pp = 1:length(params_of_interest)
            tmpidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & all_data.s_all.trialinfo(:,6)==(flip_bins(bb)*-1);
            orig_y =  all_data.s_all.(params_of_interest{pp})(tmpidx,2);
            orig_y_flip = orig_y*-1;
            all_data.s_all.(params_of_interest{pp})(tmpidx,2) = orig_y_flip;  % the y flip is inserted for CW bins.
            
            thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & (all_data.s_all.trialinfo(:,6)==flip_bins(bb) | all_data.s_all.trialinfo(:,6)==flip_bins(bb)*-1);
            % distractor bin x param x [radial; tangential] x subj
            tmp_err(bb,pp,:,ss) = nanstd( all_data.s_all.(params_of_interest{pp})(thisidx,:), [], 1 );
            tmp_mu(bb,pp,:,ss) = nanmean( all_data.s_all.(params_of_interest{pp})(thisidx,:),  1 );
            
            % new
            [dtheta_err, drad_err]= cart2pol(all_data.s_all.(params_of_interest{pp})(thisidx,1), all_data.s_all.(params_of_interest{pp})(thisidx,2));
            tmp_errpol(bb,pp,:,ss) = nanstd([drad_err rad2deg(dtheta_err)]);
            [dtheta_mu, drad_mu] = cart2pol(all_data.s_all.(params_of_interest{pp})(thisidx,1), all_data.s_all.(params_of_interest{pp})(thisidx,2));
            tmp_mupol(bb,pp,:,ss) = nanmean([drad_mu rad2deg(dtheta_mu)]); %put rad as x FIRST, theta as Y second, in keeping with later indexing 
            
            clear dtheta_err drad_err drad_err dtheta_err
            
        end
    end
end

all_err{3} = tmp_err;
all_mu{3} = tmp_mu;
all_mupol{3} = tmp_mupol;
all_errpol{3} = tmp_errpol;
clear tmp_err tmp_mu tmp_errpol tmpmupol;

% create a new container for _alldist 
tmp_err = nan(1,length(params_of_interest),2,length(subj)); % initial, final
tmp_mu  = nan(1,length(params_of_interest),2,length(subj)); % initial, final
tmp_errpol = nan(1,length(params_of_interest),2,length(subj));
tmp_mupol = nan(1,length(params_of_interest),2,length(subj));

for ss = 1:length(subj)
    thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1;
    for pp = 1:length(params_of_interest)
        % distractor bin x param x [radial; tangential] x subj
        tmp_err(1,pp,:,ss) = nanstd( all_data.s_all.(params_of_interest{pp})(thisidx,:), [], 1 );
        tmp_mu(1,pp,:,ss)  = nanmean( all_data.s_all.(params_of_interest{pp})(thisidx,:), 1 );
        [dtheta_err, drad_err]= cart2pol(all_data.s_all.(params_of_interest{pp})(thisidx,1), all_data.s_all.(params_of_interest{pp})(thisidx,2));
        tmp_errpol(1,pp,:,ss) = nanstd([drad_err rad2deg(dtheta_err)]);
        [dtheta_mu, drad_mu] = cart2pol(all_data.s_all.(params_of_interest{pp})(thisidx,1), all_data.s_all.(params_of_interest{pp})(thisidx,2));
        tmp_mupol(1,pp,:,ss) = nanmean([drad_mu rad2deg(dtheta_mu)]); %put rad as x FIRST, theta as Y second, in keeping with later indexing 
        
        clear dtheta_err drad_err drad_err dtheta_err
    end
end

all_err_alldist{1} = tmp_err;
all_mu_alldist{1} = tmp_mu;
all_errpol_alldist{1} = tmp_errpol;
all_mupol_alldist{1} = tmp_mupol;
clear tmp_err tmp_mu tmp_errpol tmpmupol;

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
    
    thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & all_data.s_all.trialinfo(:,6)==0; 

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

% new container for _alldist

tmp_rt = nan(1,length(params_of_interest),1,length(subj)); % initial, final
params_of_interest = {'i_sacc_rt'};
param_str = {'RT'};

for ss = 1:length(subj)
    
    thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1;

    for pp = 1:length(params_of_interest)
        % distractor bin x param x rt x subj
        tmp_rt(1,pp,1,ss)  = nanmean( all_data.s_all.(params_of_interest{pp})(thisidx),1);

    end
    
end

all_rt_alldist{1} =tmp_rt; 



%% no dist, near, far 
% directly compare error for no distractor, near distractor
% trials (subplot for each param); averaged over radial/tang...

params_of_interest = {'f_sacc'}; %we're using the final saccade
param_str = {'final sacc'};
cond_colors = [0.7100    0.2128    0.4772;0 0 1;0 0 1;]; % 1 = red, no dist; 2 = blue, dist 
cond_str = {'Distractor Absent','Near Dist','Far Dist'};

figure;

for pp = 1:length(params_of_interest)
    
    subplot(1,length(params_of_interest),pp); hold on;
    
    % distractor cond x subj
    thise = nan(length(all_err),length(subj)); %collect this error
    for ii = 1:length(all_err)
     
        thise(ii,:) = mean(mean(all_err{ii}(:,pp,:,:),3),1); % mean over radial/tangential; distractor bin; for n_bins, then mean over these bin
    
    end
 
    plot(1:size(thise,1),thise,'-','Color',[0.3 0.3 0.3],'LineWidth',.5);
    plot(1:size(thise,1),thise,'-','Color',[0.5 0.5 0.5],'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1);
    for ii = 1:length(all_err)
        hold on;
        tmpe = std(thise(ii,:))/sqrt(length(subj));
        plot(ii*[1 1],mean(thise(ii,:))+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(ii,:));
        plot(ii,mean(thise(ii,:)),'o','LineWidth',1.5,'Color',cond_colors(ii,:),'MarkerSize',10,'MarkerFaceColor',cond_colors(ii,:));
    clear tmpe
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
    [fval,pval] = RMAOV1_gh(x,0.05) %one-way RM ANOVA
    RMAOV1(x,0.05)
    text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',pval),'color','k','fontsize',15) %filled


end
clear tmpe
match_ylim(get(gcf,'Children'));
set(gcf,'position',[ 549   724   499   571])
%% Figure 1D : Precision, Distractor Absent vs. Distractor Present (all distractor bins collapsed) 

% directly compare avg sd for no distractor, near distractor
% trials, final sacc only ; averaged over radial/tang...

params_of_interest = {'f_sacc'}; 
param_str = {'final sacc'};
cond_colors = [0.7100    0.2128    0.4772;0 0 1;]; %red, blue
cond_str = {'Distractor Absent','Distractor Present'};

figure;

for pp = 1:length(params_of_interest) 
    
    subplot(1,length(params_of_interest),pp); hold on;
    
    % distractor cond x subj
    thise = nan(2,length(subj)); %collect this error
    for ii = 1:2
     if ii ==1
        thise(ii,:) = squeeze(mean(all_err{1}(:,pp,:,:),3)); % mean over radial/tangential; distractor bin; for n_bins, then mean over these bin
     elseif ii ==2
        thise(ii,:) = squeeze(mean(all_err_alldist{1}(:,pp,:,:),3)); % this "all_err_alldist" container is used bc all_err separates near vs. far dists. Simplified w this {} where all dist trials are within a single cell 
     end 
    end
 
    plot(1:size(thise,1),thise,'-','Color',[0.3 0.3 0.3],'LineWidth',.5);
    plot(1:size(thise,1),thise,'-','Color',[0.5 0.5 0.5],'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1);
    for ii = 1:2
        hold on;
        tmpe = std(thise(ii,:))/sqrt(length(subj));
        plot(ii*[1 1],mean(thise(ii,:))+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(ii,:));
        plot(ii,mean(thise(ii,:)),'o','LineWidth',1.5,'Color',cond_colors(ii,:),'MarkerSize',10,'MarkerFaceColor',cond_colors(ii,:));
    end
    
    xlim([0 3]);
    
    set(gca,'XTick',2,'XTickLabel',cond_str','XTickLabelRotation',45,'TickDir','out');
    xlabel('Condition');
    if pp == 1
        ylabel('Precision (avg sd, \circ)');
    end
    title(param_str{pp});
    set(gca,'XTick',1:2,'XTickLabel',cond_str','XTickLabelRotation',45,'TickDir','out');

    the_y = [thise(1,:)';thise(2,:)';];
    the_iv =[ones(length(thise(1,:)'),1); 2*ones(length(thise(1,:)'),1);];
    subj = [1 2 3 4 5 6 7]';
    the_subj =[subj;subj;];
    x = [the_y the_iv the_subj];
    [fval,pval] = RMAOV1_gh(x,0.05); %one-way RM ANOVA
    RMAOV1(x,0.05) 
    
   [p, h]= ttest(thise(1,:)',thise(2,:)')
    text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',pval),'color','k','fontsize',15) %filled

%ylim([0 1.5])
end
set(gcf,'position',[ 549   724   499   571])
match_ylim(get(gcf,'Children'));
%%ll bins 
% 
% %directly compare avg sd for no distractor, near distractor
% %trials, final sacc only ; averaged over radial/tang...
% 
params_of_interest = {'f_sacc'}; 
param_str = {'final sacc'};
cond_colors = [0.7100    0.2128    0.4772;0 0 1;0 0 1;0 0 1; 0 0 1;]; %red, blue
cond_str = {'No distractor','Bin 0 Dist','Bin 1 Dist','Bin 2 Dist','Bin 3 Dist'};

figure;

for pp = 1:length(params_of_interest) 
    
    subplot(1,length(params_of_interest),pp); hold on;
    
    % distractor cond x subj
    thise = nan(5,length(subj)); %collect this error
    
    for ii = 1:length(all_errpol)
       if ii == 1 || ii ==2
           thise(ii,:) = squeeze(mean(all_err{ii}(:,pp,:,:),3))'; % mean over radial/tangential; distractor bin; for n_bins, then mean over these bins
       else
           thise([ii:5],:)= squeeze(mean(all_err{ii}(:,pp,:,:),3)); % mean over radial/tangential; distractor bin; for n_bins, then mean over these bins
           
       end
    end
    
    plot(1:size(thise,1),thise,'-','Color',[0.3 0.3 0.3],'LineWidth',.5);
   % plot(1:size(thise,1),thise,'-','Color',[0.5 0.5 0.5],'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1);
    
    for ii = 1:5
        hold on;
        tmpe = std(thise(ii,:))/sqrt(length(subj));
        plot(ii*[1 1],mean(thise(ii,:))+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(ii,:));
        plot(ii,mean(thise(ii,:)),'o','LineWidth',1.5,'Color',cond_colors(ii,:),'MarkerSize',10,'MarkerFaceColor',cond_colors(ii,:));
    end
    
    xlim([0 6]);
    end 
    set(gca,'XTick',2,'XTickLabel',cond_str','XTickLabelRotation',45,'TickDir','out');
    xlabel('Condition');
    if pp == 1
        ylabel('Precision (avg sd,\circ)');
    end
    title(param_str{pp});
    set(gca,'XTick',1:5,'XTickLabel',cond_str','XTickLabelRotation',45,'TickDir','out');

    the_y = [thise(1,:)';thise(2,:)'; thise(3,:)'; thise(4,:)'; thise(5,:)';];
    the_iv =[ones(length(thise(1,:)'),1); 2*ones(length(thise(1,:)'),1); 3*ones(length(thise(1,:)'),1); 4*ones(length(thise(1,:)'),1); 5*ones(length(thise(1,:)'),1)];
    subj = [1 2 3 4 5 6 7]';
    the_subj =[subj;subj;subj;subj;subj];
    x = [the_y the_iv the_subj];
    subj = [1 2 3 4 5 6 7]';
    the_subj =[subj;subj;subj;subj;subj];
    x = [the_y the_iv the_subj];
    [fval,pval] = RMAOV1_gh(x,0.05); %one-way RM ANOVA
    RMAOV1(x,0.05) 
    
   [p, h]= ttest(thise(1,:)',thise(2,:)')
    text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',pval),'color','k','fontsize',15) %filled

%ylim([0 1.5])

set(gcf,'position',[ 549   724   499   571])
match_ylim(get(gcf,'Children'));
%% no dist, near, far

% directly compare avg sd for no distractor, near distractor
% trials, final sacc only ; averaged over radial/tang...

% params_of_interest = {'f_sacc'}; 
% param_str = {'final sacc'};
% cond_colors = [0.7100    0.2128    0.4772;0 0 1;0 0 1;]; %red, blue
% cond_str = {'No distractor','Near dist','Far dist'};

% new avg sd in polar angle,deg, all bins 
% 
% %directly compare avg sd for no distractor, near distractor
% %trials, final sacc only ; averaged over radial/tang...
% 
% params_of_interest = {'f_sacc'}; 
% param_str = {'final sacc'};
% cond_colors = [0.7100    0.2128    0.4772;0 0 1;0 0 1;0 0 1; 0 0 1;]; %red, blue
% cond_str = {'No distractor','Bin 0 Dist','Bin 1 Dist','Bin 2 Dist','Bin 3 Dist'};
% 
% figure;
% 
% for pp = 1:length(params_of_interest) 
%     
%     subplot(1,length(params_of_interest),pp); hold on;
%     
%     % distractor cond x subj
%     thise = nan(5,length(subj)); %collect this error
%     
%     for ii = 1:length(all_errpol)
%        if ii == 1 || ii ==2
%            thise(ii,:) = squeeze(mean(all_errpol{ii}(:,pp,:,:),3))'; % mean over radial/tangential; distractor bin; for n_bins, then mean over these bins
%        else
%            thise([ii:5],:)= squeeze(mean(all_errpol{ii}(:,pp,:,:),3)); % mean over radial/tangential; distractor bin; for n_bins, then mean over these bins
%            
%        end
%     end
%     
%     plot(1:size(thise,1),thise,'-','Color',[0.3 0.3 0.3],'LineWidth',.5);
%    % plot(1:size(thise,1),thise,'-','Color',[0.5 0.5 0.5],'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1);
%     
%     for ii = 1:5
%         hold on;
%         tmpe = std(thise(ii,:))/sqrt(length(subj));
%         plot(ii*[1 1],mean(thise(ii,:))+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(ii,:));
%         plot(ii,mean(thise(ii,:)),'o','LineWidth',1.5,'Color',cond_colors(ii,:),'MarkerSize',10,'MarkerFaceColor',cond_colors(ii,:));
%     end
%     
%     xlim([0 6]);
%     end 
%     set(gca,'XTick',2,'XTickLabel',cond_str','XTickLabelRotation',45,'TickDir','out');
%     xlabel('Condition');
%     if pp == 1
%         ylabel('Precision (avg sd, polar \circ)');
%     end
%     title(param_str{pp});
%     set(gca,'XTick',1:5,'XTickLabel',cond_str','XTickLabelRotation',45,'TickDir','out');
% 
%     the_y = [thise(1,:)';thise(2,:)';];
%     the_iv =[ones(length(thise(1,:)'),1); 2*ones(length(thise(1,:)'),1);];
%     subj = [1 2 3 4 5 6 7]';
%     the_subj =[subj;subj;];
%     x = [the_y the_iv the_subj];
%     [fval,pval] = RMAOV1_gh(x,0.05); %one-way RM ANOVA
%     RMAOV1(x,0.05) 
%     
%    [p, h]= ttest(thise(1,:)',thise(2,:)')
%     text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',pval),'color','k','fontsize',15) %filled
% 
% %ylim([0 1.5])
% 
% set(gcf,'position',[ 549   724   499   571])
% match_ylim(get(gcf,'Children'));
% %% no dist, near, far
% 
% % directly compare avg sd for no distractor, near distractor
% % trials, final sacc only ; averaged over radial/tang...
% 
% params_of_interest = {'f_sacc'}; 
% param_str = {'final sacc'};
% cond_colors = [0.7100    0.2128    0.4772;0 0 1;0 0 1;]; %red, blue
% cond_str = {'No distractor','Near dist','Far dist'};
% 
% figure;
% 
% for pp = 1:length(params_of_interest) 
%     
%     subplot(1,length(params_of_interest),pp); hold on;
%     
%     % distractor cond x subj
%     thise = nan(3,length(subj)); %collect this error
%     
%     for ii = 1:length(all_errpol)
%        % thise(ii,:) = mean(mean(all_err{ii}(:,pp,:,:),3),1); % mean over radial/tangential; distractor bin; for n_bins, then mean over these bins
%        if ii == 1 || ii ==2
%            thise(ii,:) = squeeze(mean(all_errpol{ii}(:,pp,:,:),3))'; % mean over radial/tangential; distractor bin; for n_bins, then mean over these bins
%        else
%            thise(ii,:)= mean(squeeze(mean(all_errpol{ii}(:,pp,:,:),3))); % mean over radial/tangential; distractor bin; for n_bins, then mean over these bins
%            %all dist trials are lumped into a single index, then averaged 
%        end
%     end
%     
%     plot(1:size(thise,1),thise,'-','Color',[0.3 0.3 0.3],'LineWidth',.5);
%    % plot(1:size(thise,1),thise,'-','Color',[0.5 0.5 0.5],'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1);
%     
%     for ii = 1:3
%         hold on;
%         tmpe = std(thise(ii,:))/sqrt(length(subj));
%         plot(ii*[1 1],mean(thise(ii,:))+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(ii,:));
%         plot(ii,mean(thise(ii,:)),'o','LineWidth',1.5,'Color',cond_colors(ii,:),'MarkerSize',10,'MarkerFaceColor',cond_colors(ii,:));
%         
%         clear tmpe
%     end
%     
%     xlim([0 4]);
%     end 
%     xlabel('Condition');
%     if pp == 1
%         ylabel('Precision (avg sd, polar \circ)');
%     end
%     title(param_str{pp});
%     set(gca,'XTick',1:3,'XTickLabel',cond_str','XTickLabelRotation',45,'TickDir','out');
% 
%     the_y = [thise(1,:)';thise(2,:)';thise(3,:)'];
%     the_iv =[ones(length(thise(1,:)'),1); 2*ones(length(thise(1,:)'),1);3*ones(length(thise(1,:)'),1);];
%     subj = [1 2 3 4 5 6 7]';
%     the_subj =[subj;subj;subj];
%     x = [the_y the_iv the_subj];
%     [fval,pval] = RMAOV1_gh(x,0.05); %one-way RM ANOVA
%     RMAOV1(x,0.05) 
%   
%     text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',pval),'color','k','fontsize',15) %filled
% 
% 
% 
% set(gcf,'position',[ 549   724   499   571])
% match_ylim(get(gcf,'Children'));
%% new
%% Figure 1D : Precision, Distractor Absent vs. Distractor Present (all distractor bins collapsed) 

% % directly compare avg sd for no distractor, near distractor
% % trials, final sacc only ; averaged over radial/tang...
% 
% params_of_interest = {'f_sacc'}; 
% param_str = {'final sacc'};
% cond_colors = [0.7100    0.2128    0.4772;0 0 1;]; %red, blue
% cond_str = {'Distractor Absent','Distractor Present'};
% 
% figure;
% 
% for pp = 1:length(params_of_interest) 
%     
%     subplot(1,length(params_of_interest),pp); hold on;
%     
%     % distractor cond x subj
%     thise = nan(2,length(subj)); %collect this error
%     for ii = 1:2
%      if ii ==1
%         thise(ii,:) = squeeze(mean(all_errpol{1}(:,pp,:,:),3)); % mean over radial/tangential; distractor bin; for n_bins, then mean over these bin
%      elseif ii ==2
%         thise(ii,:) = squeeze(mean(all_errpol_alldist{1}(:,pp,:,:),3)); % this "all_err_alldist" container is used bc all_err separates near vs. far dists. Simplified w this {} where all dist trials are within a single cell 
%      end 
%     end
%  
%     plot(1:size(thise,1),thise,'-','Color',[0.3 0.3 0.3],'LineWidth',.5);
%     plot(1:size(thise,1),thise,'-','Color',[0.5 0.5 0.5],'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1);
%     for ii = 1:2
%         hold on;
%         tmpe = std(thise(ii,:))/sqrt(length(subj));
%         plot(ii*[1 1],mean(thise(ii,:))+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(ii,:));
%         plot(ii,mean(thise(ii,:)),'o','LineWidth',1.5,'Color',cond_colors(ii,:),'MarkerSize',10,'MarkerFaceColor',cond_colors(ii,:));
%         
%         clear tmpe
%     end
%     
%     xlim([0 3]);
%     
%     set(gca,'XTick',1:2,'XTickLabel',cond_str','XTickLabelRotation',45,'TickDir','out');
%     xlabel('Condition');
%     if pp == 1
%         ylabel('Precision (avg sd, polar \circ)');
%     end
%     title(param_str{pp});
%     set(gca,'XTick',1:2,'XTickLabel',cond_str','XTickLabelRotation',45,'TickDir','out');
% 
%     the_y = [thise(1,:)';thise(2,:)';];
%     the_iv =[ones(length(thise(1,:)'),1); 2*ones(length(thise(1,:)'),1);];
%     subj = [1 2 3 4 5 6 7]';
%     the_subj =[subj;subj;];
%     x = [the_y the_iv the_subj];
%     [fval,pval] = RMAOV1_gh(x,0.05); %one-way RM ANOVA
%     RMAOV1(x,0.05) 
%     
%     [p, h]= ttest(thise(1,:)',thise(2,:)')
%     text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',pval),'color','k','fontsize',15) %filled
% 
% %ylim([0 1.5])
% end
% set(gcf,'position',[ 549   724   499   571])
% match_ylim(get(gcf,'Children'));
%% REPORTED STATS ONLY - PRECISION, Distractor Absent, Abs(Bin) Distractor Present (4-level RM ANOVA) 
params_of_interest = {'f_sacc'}; %we're using the final saccade
param_str = {'final sacc'};
cond_colors = [0.7100    0.2128    0.4772;0 0 1;0 0 1;0 0 1;0 0 1;];
cond_str = {'No distractor','Bin 0 Dist','Bin 1 Dist','Bin 2 Dist','Bin 3 Dist'};
%(F(2,4) = 2.081, p = 0.1148
figure;

for pp = 1:length(params_of_interest)
    
    subplot(1,length(params_of_interest),pp); hold on;
    
    % distractor cond x subj
    thise = nan(5,length(subj)); %collect this error
    
    for ii = 1:length(all_err)
       % thise(ii,:) = mean(mean(all_err{ii}(:,pp,:,:),3),1); % mean over radial/tangential; distractor bin; for n_bins, then mean over these bins
       if ii == 1 || ii ==2
           thise(ii,:) = squeeze(mean(all_err{ii}(:,pp,:,:),3))'; % mean over radial/tangential; distractor bin; for n_bins, then mean over these bins
       else
           thise([ii:5],:)= squeeze(mean(all_err{ii}(:,pp,:,:),3)); % mean over radial/tangential; distractor bin; for n_bins, then mean over these bins
           
       end
    end
    
    plot(1:size(thise,1),thise,'-','Color',[0.3 0.3 0.3],'LineWidth',.5);
    plot(1:size(thise,1),thise,'-','Color',[0.5 0.5 0.5],'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1);
    
    for ii = 1:5
        hold on;
        tmpe = std(thise(ii,:))/sqrt(length(subj));
        plot(ii*[1 1],mean(thise(ii,:))+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(ii,:));
        plot(ii,mean(thise(ii,:)),'o','LineWidth',1.5,'Color',cond_colors(ii,:),'MarkerSize',10,'MarkerFaceColor',cond_colors(ii,:));
    end
    
    xlim([0 6]);
    
    set(gca,'XTick',1:5,'XTickLabel',cond_str','XTickLabelRotation',45,'TickDir','out');
    xlabel('Condition');
    if pp == 1
        ylabel('comp Precision (avg sd, \circ)');
    end
    title(param_str{pp});
    
    the_y = [thise(1,:)';thise(2,:)'; thise(3,:)'; thise(4,:)'; thise(5,:)';];
    the_iv =[ones(length(thise(1,:)'),1); 2*ones(length(thise(1,:)'),1); 3*ones(length(thise(1,:)'),1); 4*ones(length(thise(1,:)'),1); 5*ones(length(thise(1,:)'),1)];
    subj = [1 2 3 4 5 6 7]';
    the_subj =[subj;subj;subj;subj;subj];
    x = [the_y the_iv the_subj];
    RMAOV1(x,0.05) %one-way RM ANOVA
    [fval,pval] = RMAOV1_gh(x,0.05) %one-way RM ANOVA
    RMAOV1(x,0.05)

    text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',pval),'color','k','fontsize',15) %filled


end

match_ylim(get(gcf,'Children'));
set(gcf,'position',[ 549   724   499   571])

% %% RT, no dist to all_dist
% 
% 
% params_of_interest = {'i_sacc_rt'}; %we're using the final saccade
% param_str = {'RT'};
% cond_colors = [0.7100    0.2128    0.4772;0 0 1;];
% cond_str = {'Distractor Absent','Distractor Present'};
% 
% figure;
% 
% for pp = 1:length(params_of_interest)
%     
%     subplot(1,length(params_of_interest),pp); hold on;
%     
%     % distractor cond x subj
%     thise = nan(2,length(subj)); %collect this error
%     for ii = 1:2
%      if ii ==1
%         thise(ii,:) = squeeze(mean(all_rt{1}(:,pp,:,:),3)); % mean over radial/tangential; distractor bin; for n_bins, then mean over these bin
%      elseif ii ==2
%        thise(ii,:) = squeeze(mean(all_rt_alldist{1}(:,pp,:,:),3)); % mean over radial/tangential; distractor bin; for n_bins, then mean over these bin
%      end 
%     end
%     
%     
%     
%     plot(1:size(thise,1),thise,'-','Color',[0.3 0.3 0.3],'LineWidth',.5);
%     plot(1:size(thise,1),thise,'-','Color',[0.5 0.5 0.5],'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1);
%     for ii = 1:2
%         hold on;
%         tmpe = std(thise(ii,:))/sqrt(length(subj));
%         plot(ii*[1 1],mean(thise(ii,:))+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(ii,:));
%         plot(ii,mean(thise(ii,:)),'o','LineWidth',1.5,'Color',cond_colors(ii,:),'MarkerSize',10,'MarkerFaceColor',cond_colors(ii,:));
%     end
%     
%     xlim([0 3]);
%     
%     set(gca,'XTick',1:2,'XTickLabel',cond_str','XTickLabelRotation',45,'TickDir','out');
%     xlabel('Condition');
%     if pp == 1
%         ylabel('RT (seconds)');
%     end
%     title(param_str{pp});
%     
%     the_y = [thise(1,:)';thise(2,:)';];
%     the_iv =[ones(length(thise(1,:)'),1); 2*ones(length(thise(1,:)'),1);];
%     subj = [1 2 3 4 5 6 7]';
%     the_subj =[subj;subj;]
%     x = [the_y the_iv the_subj];
%     [fval,pval] = RMAOV1_gh(x,0.05) %one-way RM ANOVA
%     RMAOV1(x,0.05) 
%     text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',pval),'color','k','fontsize',15) %filled
% 
% 
% end
% set(gcf,'position',[ 549   724   499   571])
% match_ylim(get(gcf,'Children'));
% %% RT by bin 
% 
% params_of_interest = {'i_sacc_rt'}; %we're using the final saccade
% param_str = {'RT'};
% cond_colors = [0.7100    0.2128    0.4772;0 0 1;0 0 1;0 0 1;0 0 1;];
% cond_str = {'No distractor','Bin 0 Dist','Bin 1 Dist','Bin 2 Dist','Bin 3 Dist'};
% 
% figure;
% 
% for pp = 1:length(params_of_interest)
%     
%     subplot(1,length(params_of_interest),pp); hold on;
%     
%     % distractor cond x subj
%     thise = nan(5,length(subj)); %collect this error
%     for ii = 1:length(all_rt)
%        if ii == 1 || ii ==2
%            thise(ii,:) = squeeze(mean(all_rt{ii}(:,pp,:,:),3))'; % mean over radial/tangential; distractor bin; for n_bins, then mean over these bins
%            
%        else
%            thise([ii:5],:)= squeeze(mean(all_rt{ii}(:,pp,:,:),3)); % mean over radial/tangential; distractor bin; for n_bins, then mean over these bins
%            
%        end
%     end
%     
%     plot(1:size(thise,1),thise,'-','Color',[0.3 0.3 0.3],'LineWidth',.5);
%     plot(1:size(thise,1),thise,'-','Color',[0.5 0.5 0.5],'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1);
%     for ii = 1:5
%         hold on;
%         tmpe = std(thise(ii,:))/sqrt(length(subj));
%         plot(ii*[1 1],mean(thise(ii,:))+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(ii,:));
%         plot(ii,mean(thise(ii,:)),'o','LineWidth',1.5,'Color',cond_colors(ii,:),'MarkerSize',10,'MarkerFaceColor',cond_colors(ii,:));
%     end
%     
%     xlim([0 6]);
%     
%     set(gca,'XTick',1:5,'XTickLabel',cond_str','XTickLabelRotation',45,'TickDir','out');
%     xlabel('Condition');
%     if pp == 1
%         ylabel('RT');
%     end
%     title(param_str{pp});
%     
%     the_y = [thise(1,:)';thise(2,:)'; thise(3,:)'; thise(4,:)'; thise(5,:)';];
%     the_iv =[ones(length(thise(1,:)'),1); 2*ones(length(thise(1,:)'),1); 3*ones(length(thise(1,:)'),1); 4*ones(length(thise(1,:)'),1); 5*ones(length(thise(1,:)'),1)];
%     subj = [1 2 3 4 5 6 7]';
%     the_subj =[subj;subj;subj;subj;subj];
%     x = [the_y the_iv the_subj];
%    [fval,pval] = RMAOV1_gh(x,0.05) %one-way RM ANOVA
%     RMAOV1(x,0.05)
%     text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',pval),'color','k','fontsize',15) %filled
% 
% 
% end
% 
% match_ylim(get(gcf,'Children'));
% 
% %% Figure 1E : RT
% 
% cond_colors = [0.7100    0.2128    0.4772;0 0 1;0 0 1;0 0 1;0 0 1;]; %colors for illustrator 
% params_of_interest = {'i_sacc_rt'};
% param_str = {'RT'};
% 
% figure;
% 
% for pp = 1:length(params_of_interest)
%     
%     subplot(1,length(params_of_interest),pp); hold on;                                                                                    
%     
%     % nbins x nsubj
%     thisnod_rt = squeeze(all_rt{1}(1,pp,1,:))';
%     thisnear_dist = squeeze(all_rt{2}(:,pp,1,:))';
%     thisfar_dist = mean(squeeze(all_rt{3}(:,pp,1,:)));
%     thisrt= [thisnod_rt; thisnear_dist; thisfar_dist];
%     plot(1:length(flip_bins),thisrt,'-','color',[.2 .2 .2]);
%     
%     for ii = 1:size(thisrt,1)
%         tmpe = std(thisrt(ii,:))/sqrt(length(subj));
%         plot(ii*[1 1],mean(thisrt(ii,:))+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(ii,:));
%         plot(ii,mean(thisrt(ii,:)),'o','LineWidth',1.5,'color',cond_colors(ii,:),'markerfacecolor',cond_colors(ii,:),'MarkerSize',10);
%     end
%     title(param_str{pp});
%     xlabel('Distractor Condition');
%    
%     if pp == 1
%         ylabel('RT');
%         
%     end
%     set(gca,'XTick',[1 2 3],'Xticklabels',({'No Distractor', 'Near Distractor','Far Distractor'}), 'TickDir','out','XTickLabelRotation',45)
%     the_RT = [thisrt(1,:)'; thisrt(2,:)'; thisrt(3,:)']; %use for just dist
%     the_ivRT =[ones(length(thisrt(1,:)'),1); 2*ones(length(thisrt(2,:)'),1); 3*ones(length(thisrt(3,:)'),1);]; %use for just dist
%     the_subj =[subj;subj;subj;];
%     x = [the_RT the_ivRT the_subj];
%     RMAOV1(x,0.05)
% end
% match_ylim(get(gcf,'Children'));
% xlim([0 size(thisrt,1)+1]);
% 
%% Figure 1F : Bias 
params_of_interest = {'f_sacc'};
param_str = {'final sacc'};

figure;

for pp = 1:length(params_of_interest)
    
    subplot(1,length(params_of_interest),pp); hold on;
    
    % nbins x nsubj
    thism_near = squeeze(all_mu{2}(:,pp,2,:))'; %this is the only condition that is information for this analysis.
    plot(1,thism_near,'o','MarkerSize',5,'Markerfacecolor',[0.3 0.3 0.3], 'Color',[0.3 0.3 0.3]);
    plot(1,mean(thism_near,2),'o','MarkerSize',15,'Color',cond_colors(2,:),'MarkerFaceColor',cond_colors(2,:));
    hold on;
    
    tmpe = std(thism_near)/sqrt(length(subj));
    plot([1 1], mean(thism_near)+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(2,:));

    %errorbar(1, mean(thism_near), [std(thism_near)/sqrt(length(subj))],'linewidth',2,'Color',cond_colors(2,:))
    title(param_str{pp});
    
   
    ylabel('Memory Bias, toward distractor (\circ)');
end
xlabel('Near Distractor');

[h_bias_near p_bias_near,CI,STATS] = ttest(thism_near(1,:)') %t-test to determine if bias is different than zero
text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',p_bias_near),'color','k','fontsize',15) %filled

match_ylim(get(gcf,'Children'));
set(gcf,'position',[ 549   724   499   571])

% %% Figure 1F : Bias %% NEW%%
% params_of_interest = {'f_sacc'};
% param_str = {'final sacc'};
% 
% figure;
% 
% for pp = 1:length(params_of_interest)
%     
%     subplot(1,length(params_of_interest),pp); hold on;
%     
%     % nbins x nsubj
%     thism_near = squeeze(all_mupol{2}(:,pp,2,:))'; %this is the only condition that is information for this analysis.
%     plot(1,thism_near,'o','MarkerSize',5,'Markerfacecolor',[0.3 0.3 0.3], 'Color',[0.3 0.3 0.3]);
%     plot(1,mean(thism_near,2),'o','MarkerSize',15,'Color',cond_colors(2,:),'MarkerFaceColor',cond_colors(2,:));
%     hold on;
%     
%     tmpe = std(thism_near)/sqrt(length(subj));
%     plot([1 1], mean(thism_near)+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(2,:));
% 
%     %errorbar(1, mean(thism_near), [std(thism_near)/sqrt(length(subj))],'linewidth',2,'Color',cond_colors(2,:))
%     title(param_str{pp});
%     
%    
%     ylabel('Memory Bias, toward distractor (polar \circ)');
% end
% xlabel('Near Distractor');
% 
% [h_bias_near p_bias_near,CI,STATS] = ttest(thism_near(1,:)') %t-test to determine if bias is different than zero
% text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',p_bias_near),'color','k','fontsize',15) %filled
% 
% match_ylim(get(gcf,'Children'));
%set(gcf,'position',[ 549   724   499   571])
