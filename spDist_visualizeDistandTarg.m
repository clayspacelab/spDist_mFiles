% spDist_plotEyeData_visualizeDistandTarg.m


root = '/share/data/spDist/'

%subj = {'KD','CC','AY','MR','XL','EK','SF'};
subj = {'CC','AY','MR'};
%sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist2'},{'spDist1','spDist2'},{'spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}}; %two sessions removed
sess = {{'spDistLong1','spDistLong2'},{'spDistLong1','spDistLong2'},{'spDistLong1','spDistLong2'}}
%if nargin < 3
%WHICH_EXCL = [13 20 21 22]; % don't exclude trials w/ calibration failures for now...
WHICH_EXCL = []; % don't exclude trials w lg error, this is overridden due to geh inspection
%end


% for now, let's use cat_struct to load/concatenate all data...
all_subj = nan(1000*length(subj),1);
%u_subj = unique(cellfun(@(s) s(1:2),subj,'uniformoutput',0));

all_data = [];

startidx = 1;

for ss = 1:length(subj)
    for sessidx = 1:length(sess{ss})
        
        fn = sprintf('%s/spDist_behav_61520/%s_%s_scored.mat',root,subj{ss},sess{ss}{sessidx});
        fprintf('Loading scored eye data from %s\n',fn);
        this_scored = load(fn);
        
        this_data.s_all = this_scored.ii_sess;
        this_data.sess_all = sessidx;
        
        this_subj = ss;%find(strcmpi(u_subj,subj{ss}(1:2)));
        
        all_data = cat_struct(all_data,this_data);
        all_subj(startidx:(startidx-1+size(this_scored.ii_sess.trialinfo,1))) = this_subj;
        
        startidx = startidx+size(this_scored.ii_sess.trialinfo,1);
        
        clear this_subj this_data;
    end
end

% let's try this pattern for now
all_subj = all_subj(1:(startidx-1));
all_data.subj_all = all_subj;

% determine which trials to include
% first, narrow based on saccade preprocessing/scoring exclusions
% (wmChoose_extractSaccadeData1.m)
all_data.use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, WHICH_EXCL), all_data.s_all.excl_trial, 'UniformOutput',false));
all_data.use_trial(all_data.s_all.f_sacc_err>10) = 0;
all_data.use_trial(all_data.s_all.i_sacc_err>10) = 0;
all_data.use_trial(ismember(all_data.s_all.r_num, [8,12,14]) & all_subj==3) = 0;
% drop trials with very short (< 100 ms) or very long RT (> 1 s)
all_data.use_trial(all_data.s_all.i_sacc_rt<0.1 | all_data.s_all.i_sacc_rt>1.0) = 0;


%% collect error (std dev) for ccw & cw distractors
%  for distractor bin zero trials, compare distractor offset


which_bin = 1;

% within each:
% - row: bin
% - col: initial, final
% - page1:rad/tang (x/y)
% - page2:subj

all_mu = cell(1,1);

cond_colors = lines(2);

params_of_interest = {'f_sacc'};
param_str = {'final saccade'};

testmu_cw = nan(1,length(params_of_interest),2,length(subj));
tmp_mu = nan(2,length(params_of_interest),2,length(subj));
figure(1);
figure(2);

for ss = 1:length(subj)
    
    cwidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & all_data.s_all.trialinfo(:,6)== which_bin*-1 & all_data.s_all.trialinfo(:,10) < 0; %negative jitter = cw
    ccwidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & all_data.s_all.trialinfo(:,6)== which_bin & all_data.s_all.trialinfo(:,10) > 0; %positive jitter = ccw
    
 
    for pp = 1:length(params_of_interest)
        tmp_mu(1,pp,:,ss) = nanmean( all_data.s_all.(params_of_interest{pp})(cwidx,:),  1 );
        tmp_mu(2,pp,:,ss) = nanmean( all_data.s_all.(params_of_interest{pp})(ccwidx,:),  1 );
    end
    
    for ii =1:length(find(cwidx))
        tmp = find(cwidx);
        thistmp =tmp(ii);
        
        figure(1)
        subplot(1,length(subj),ss)
        hold on;
        plot([all_data.s_all.trialinfo(thistmp,2) all_data.s_all.trialinfo(thistmp,4)], [all_data.s_all.trialinfo(thistmp,3) all_data.s_all.trialinfo(thistmp,5)],'bo-');
        plot(all_data.s_all.trialinfo(thistmp,4),all_data.s_all.trialinfo(thistmp,5),'ro','markersize',20) % only plot the dist now at a diff size
        title(subj(ss))
    end
    %sgtitle('all clockwise distractor and targ positions')
    sgtitle(sprintf('bin %i distractor distractor and targ positions, neg jitter',which_bin))
    for ii =1:length(find(ccwidx));
        tmp = find(ccwidx);
        thattmp =tmp(ii);
        figure(2)
        subplot(1,length(subj),ss)
        hold on;
        plot([all_data.s_all.trialinfo(thattmp,2) all_data.s_all.trialinfo(thattmp,4)], [all_data.s_all.trialinfo(thattmp,3) all_data.s_all.trialinfo(thattmp,5)],'bo-');
        plot(all_data.s_all.trialinfo(thattmp,4),all_data.s_all.trialinfo(thattmp,5),'ro','markersize',10)
        title(subj(ss))
    end
    %sgtitle('all counterclockwise distractor and targ positions')
    
    sgtitle(sprintf('bin %i distractor and targ positions, pos jitter',which_bin))
end


%% do a test here to visualize how uniform bin 3 is 


which_bin = 3;

% within each:
% - row: bin
% - col: initial, final
% - page1:rad/tang (x/y)
% - page2:subj

all_mu = cell(1,1);

cond_colors = lines(3);

params_of_interest = {'f_sacc'};
param_str = {'final saccade'};

testmu_cw = nan(1,length(params_of_interest),2,length(subj));
tmp_mu = nan(2,length(params_of_interest),2,length(subj));
figure(3);
tmp_ccwidxstore = {};
tmp_cwidxstore = {};

% same plot
for ss = 1:length(subj)
    
    cwidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.s_all.trialinfo(:,6)== which_bin*-1; %negative jitter = cw
    ccwidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2  & all_data.s_all.trialinfo(:,6)== which_bin; %positive jitter = ccw
    
 
    for pp = 1:length(params_of_interest)
        tmp_mu(1,pp,:,ss) = nanmean( all_data.s_all.(params_of_interest{pp})(cwidx,:),  1 );
        tmp_mu(2,pp,:,ss) = nanmean( all_data.s_all.(params_of_interest{pp})(ccwidx,:),  1 );
    end
    
    for ii =1:length(find(cwidx))
        if ii==1
        sprintf('%i trials for subj %s', length(find(cwidx)), subj{ss})
        else
        end
        
        tmp = find(cwidx);
        tmp_cwidxstore{ss} = tmp;
        thistmp =tmp(ii);
        
        figure(3)
        hold on;
        %plot([all_data.s_all.trialinfo(thistmp,2) all_data.s_all.trialinfo(thistmp,4)], [all_data.s_all.trialinfo(thistmp,3) all_data.s_all.trialinfo(thistmp,5)],'bo-');
        plot(all_data.s_all.trialinfo(thistmp,4),all_data.s_all.trialinfo(thistmp,5),'ro','markersize',10) % only plot the dist now at a diff size
    end
    clear tmp
    for ii =1:length(find(ccwidx));
        if ii ==1
        sprintf('%i trials for subj %s', length(find(cwidx)), subj{ss})
        else
        end
        
        tmp = find(ccwidx);
        tmp_ccwidxstore{ss} = tmp;

        thattmp =tmp(ii)
        figure(3)
        hold on;
        %plot([all_data.s_all.trialinfo(thattmp,2) all_data.s_all.trialinfo(thattmp,4)], [all_data.s_all.trialinfo(thattmp,3) all_data.s_all.trialinfo(thattmp,5)],'bo-');
        plot(all_data.s_all.trialinfo(thattmp,4),all_data.s_all.trialinfo(thattmp,5),'ro','markersize',40)
    end
    
    sgtitle(sprintf('all bin %i distractor positions',which_bin))
end


%% do a test here to visualize how uniform bin 3 is 


which_bin = 3;

% within each:
% - row: bin
% - col: initial, final
% - page1:rad/tang (x/y)
% - page2:subj

all_mu = cell(1,1);

cond_colors = lines(3);


testmu_cw = nan(1,length(params_of_interest),2,length(subj));
tmp_mu = nan(2,length(params_of_interest),2,length(subj));
figure(3);
tmp_ccwidxstore = {};
tmp_cwidxstore = {};

% to confirm TARGET in RIGHT VF, aka trials to be collected for LEFT
% HEMISPHERE VOXELS 
for ss = 1:length(subj)
  thisidx = all_subj==ss & (all_data.s_all.trialinfo(:,6) == 3 |all_data.s_all.trialinfo(:,6) == -3) & (all_data.s_all.trialinfo(:,2)>0 & all_data.s_all.trialinfo(:,4)<0) ;

    for ii =1:length(find(thisidx))
        if ii==1
        sprintf('%i trials for subj %s', length(find(thisidx)),subj{ss})
        else
        end
        
        tmp = find(thisidx);
        thistmp =tmp(ii);
        
        figure(4)
        hold on;
        plot([all_data.s_all.trialinfo(thistmp,2) all_data.s_all.trialinfo(thistmp,4)], [all_data.s_all.trialinfo(thistmp,3) all_data.s_all.trialinfo(thistmp,5)],'bo-');
        plot(all_data.s_all.trialinfo(thistmp,4),all_data.s_all.trialinfo(thistmp,5),'ro','markersize',30)
    end
    if ss==3
       sgtitle('all bin 3 distractor positions, target in RVF')
    else
    end
end

clear thisidx ss ii
% to confirm TARGET in RIGHT VF, aka trials to be collected for LEFT
% HEMISPHERE VOXELS 
for ss = 1:length(subj)
  thisidx = all_subj==ss & (all_data.s_all.trialinfo(:,6) == 1 |all_data.s_all.trialinfo(:,6) == -1 |all_data.s_all.trialinfo(:,6) == 2 |all_data.s_all.trialinfo(:,6) == -2 | all_data.s_all.trialinfo(:,6) == 3 |all_data.s_all.trialinfo(:,6) == -3) & (all_data.s_all.trialinfo(:,2)<0 & all_data.s_all.trialinfo(:,4)>0) ;

    for ii =1:length(find(thisidx))
        if ii==1
        sprintf('%i trials for subj %s', length(find(thisidx)),subj{ss})
        else
        end
        
        tmp = find(thisidx);
        thistmp =tmp(ii);
        
        figure(8)
        hold on;
        plot([all_data.s_all.trialinfo(thistmp,2) all_data.s_all.trialinfo(thistmp,4)], [all_data.s_all.trialinfo(thistmp,3) all_data.s_all.trialinfo(thistmp,5)],'bo-');
        plot(all_data.s_all.trialinfo(thistmp,4),all_data.s_all.trialinfo(thistmp,5),'ro','markersize',30)
    end
     if ss==3
       sgtitle('all bin 3 distractor positions, target in LVF')
    else
    end
end

