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
        fn = sprintf('%sspDist_behav_92220/%s_%s_scored.mat',root,subj{ss},sess{ss}{sessidx}); 
        %fn = sprintf('%sspDist_behav_61520/%s_%s_scored.mat',root,subj{ss},sess{ss}{sessidx}); 
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

acc_mu = mean(all_data.s_all.trialinfo(:,11));%idx of task acc 
coh_mu = mean(all_data.s_all.trialinfo(:,12))
%% what was acc by bin?

dist_bins =[-3 -2 -1 0 1 2 3];

bin_acc =nan(length(subj),length(dist_bins));

for ss =1:length(subj)
    for dd =1:length(dist_bins)
        
        thisidx =  all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & all_data.s_all.trialinfo(:,6)==dist_bins(dd); %negative CW jitter

        bin_acc(ss,dd) = mean(all_data.s_all.trialinfo(thisidx,8)); 
     %  bin_acc(ss,dd) = sum(all_data.s_all.trialinfo(thisidx,8))/length(all_data.s_all.trialinfo(thisidx,8));
    end
end


figure
plot([1:length(dist_bins)], bin_acc','-', 'color',[.5 .5 .5])
hold on;

for ii=1:length(dist_bins)
    sem = std(bin_acc(:,ii),1)/length(subj)
    plot(ii, mean(bin_acc(:,ii)),'ko','markersize',10)
    plot([ii ii], [mean(bin_acc(:,ii))+1.*sem mean(bin_acc(:,ii))-1.*sem], 'k-','linewidth',2)
    clear sem
end

the_acc = [bin_acc(:,1); bin_acc(:,2); bin_acc(:,3); bin_acc(:,4); bin_acc(:,5); bin_acc(:,6); bin_acc(:,7)]; %use for just dist
    the_ivacc =[ones(length(bin_acc(:,1)),1); 2*ones(length(bin_acc(:,2)),1); 3*ones(length(bin_acc(:,3)),1); 4*ones(length(bin_acc(:,4)),1); 5*ones(length(bin_acc(:,5)),1); 6*ones(length(bin_acc(:,6)),1); 7*ones(length(bin_acc(:,7)),1)]; %use for just dist
    subj =[ 1 2 3 4 5 6 7]';
    the_subj =[subj;subj;subj;subj; subj; subj; subj];
    x = [the_acc the_ivacc the_subj];
    RMAOV1(x,0.05)
set(gca,'xtick',[1:7], 'xticklabel', {'Bin -3', 'Bin -2', 'Bin -1', 'Bin 0', 'Bin 1', 'Bin 2', 'Bin 3'},'xlim',[0 8],'xticklabelrotation',45)
ylabel('% correct, distractor response')
title('Accuracy by bin')