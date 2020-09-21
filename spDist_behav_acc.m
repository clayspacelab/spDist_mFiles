%spDist_behav_acc

root = '/share/data/spDist/';

%load raw subject data 
subj = {'KD','CC','AY','MR','XL','EK','SF'};
sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist2'},{'spDist1','spDist2'},{'spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}}; %two sessions removed

WHICH_EXCL = [13 20 21 22];


% concatenate ALL subject data
all_subj = nan(1000*length(subj),1);

all_data = [];

startidx = 1;

for ss = 1:length(subj)
    for sessidx = 1:length(sess{ss})
        
        fn = sprintf('%sspDist_behav_61520/%s_%s_scored.mat',root,subj{ss},sess{ss}{sessidx}); 
        %fn = sprintf('%s/spDist_behav_82719/%s_%s_scored.mat',root,subj{ss},sess{ss}{sessidx});
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
%% what was acc?

acc_mu = mean(all_data.s_all.trialinfo(:,11)) %idx of task acc 