% spDist_pilot_concatBehav.m
%
% combines all the behavioral data files into one for easy loading
%
% NOTE: for now, requires _scoreEyeData.m be run already - will load
% trialinfo from ii_sess file (_scored.mat) since that's where we compute
% most of the condition info, and no need to replicate...

function spDist_concatBehav(subj,sess)

task_dir = 'spDist'; % use this to construct all dirs

root = spDist_loadRoot;

root_raw = sprintf('%s/raw',root); % where raw behavioral data lives
root_target = sprintf('%s/%s_behav',root,task_dir); % where to save concat file (and where to look for _scored)

if nargin < 1 || isempty(subj)
    subj = {'AY','CC','KD','MR','XL'};
end


if nargin < 2 || isempty(sess)
    sess_template = {'spDist1','spDist2'};
    sess = cell(length(subj),1); for ss = 1:length(subj); sess{ss} = sess_template; end
    clear sess_template
end




for ss = 1:length(subj)
    
    for sess_idx = 1:length(sess{ss})
        
        this_subjID = sprintf('%s_%s',subj{ss},sess{ss}{sess_idx});
        
%	if strfind(sess{ss}{sess_idx},'Long')
%		expt_name = '';
%	else
        	expt_name = 'spDist_scanner';
%	end
        
        % look for all behavioral files (not saccades yet) 20171213T103006
        fnm_b = sprintf('%s/%s_behav/%s_r*_%s_*.mat',root_raw,this_subjID,this_subjID,expt_name); % because of datestr(now,30), there's a T in the behav files, but not in the saccade files
        tmp_f_b = dir(fnm_b);
        
        % get rid of "runTaskMap" & "iEye" files if they exist
        myfm_b = tmp_f_b( ~contains( {tmp_f_b.name}, {'runTaskMap','iEye'} ) );

        
        
        % initialize the variables we'll fill up
        t_all = [];
        tnum_all = [];
        c_all = []; % from trialinfo, maybe move it there...
        p_all.ecc = []; % is this saved in an obvious place? maybe load parameters.mat now...
        
        r_all = [];
        targ_ang_all = [];
        dist_ang_all = [];
        
        
        
%         t_map = [];   % timing
%         tnum_map = []; % trial number within each run
%         c_map = [];   % conditions/trial labels
%         p_map.radMean = []; % stimulus parameters
%         p_map.polarAngleOffset = [];
%         r_map = [];  % run #
%         coords_map = [];
        
        % (these should be same length - should probably check for that here)
        for ff = 1:length(myfm_b)
            
            fnb = sprintf('%s/%s_behav/%s',root_raw,this_subjID,myfm_b(ff).name);
            fprintf('loading %s...\n',fnb);
            bdata = load(fnb);
            bdata = bdata.p;
            
            ntrials = bdata.ntrials;
            
            r_all = [r_all; ff*ones(ntrials,1)];
            
            % compute timing - we want timing of delay onset, distractor
            % onset, response onset
            % col 1: onset of Delay 1
            % col 2: onset of distractor (for KD, sess1, not saved on
            % distractor trials!!!!)
            % col 3: onset of MGS cue
            t_tmp = nan(ntrials,3);
            t_tmp(:,1) = bdata.delay_start;%-bdata.raw_t.startExperimentTime;
            t_tmp(:,2) = bdata.dist_start(:,1);%raw_t.distractorStartTime-bdata.raw_t.startExperimentTime;
            t_tmp(:,3) = bdata.MGS_start;%raw_t.respStartTime-bdata.raw_t.startExperimentTime-0.25; % NOTE: 0.25 offset here is due to double-marking of events for some reason, go cue occurred 0.25 s before this...
            
            t_all = [t_all;t_tmp-bdata.expt_start];
            
            clear t_tmp;
            
            % we know all runs/trials have same ecc, so let's just get it
            % from the first one
            if ff == 1
                p_all.ecc = bdata.wm_ecc;
            end
            
            clear bdata;
        end
        
        % load the scored eyedata
        efn = sprintf('%s/%s_scored.mat',root_target,this_subjID);
        edata = load(efn);
        
        c_all = edata.ii_sess.trialinfo;
        
        clear edata;
        
        
        % compute angle for target, distractor using 2,3 and 4,5 columns
        targ_ang_all = atan2d(c_all(:,3),c_all(:,2));
        dist_ang_all = atan2d(c_all(:,5),c_all(:,4));
        
        % save everything
        fn2s = sprintf('%s/%s_behav.mat',root_target,this_subjID);
        fprintf('saving to %s...\n',fn2s);
        save(fn2s,'tnum_all','c_all','p_all','r_all','targ_ang_all','dist_ang_all','t_all');
        
        % clear those things..
        clear tnum_all c_all p_all r_all targ_ang_all dist_ang_all t_all;
    end
end

return
