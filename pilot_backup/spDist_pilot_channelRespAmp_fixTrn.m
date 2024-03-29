% spDist_pilot_channelRespAmp_fixTrn.m
% adapted from MGSMap_channelRespAmp_catSess_thruTime_fixTrn.m
%
% trains on all MGSMap runs found in trn_dir, using a fixed training time
% window, then reconstructs on all testing timepoints in tst_sess. voxel
% selection, if based on n_vox, is done using training data only
%
% computes reconstructions at each timepoint using a FIXED encoding model,
% estimated using specified timepoints
%
% input "sess" refers to TESTING sessions for now...
%
% TCS 8/23/2018
%
function spDist_pilot_channelRespAmp_fixTrn(subj,sess,ROIs,which_vox,trn_tpts)

tst_dir = 'spDist_pilot';
trn_dir = 'wmChoose';

trn_sess = 'MGSMap'; % files to load for training

root =  spDist_pilot_loadRoot;% '/deathstar/data/PrismaPilotScans/';
trn_root = sprintf('%s/../wmChoose_scanner/',root);

if nargin < 1
    subj = {'KD'};
end
if nargin < 2
    sess = {{'spDist_pilot1'}};
end

if nargin < 3
    ROIs = {'V1','V2','V3','V3AB','hV4','VO1','VO2','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS','iPCS'};
    
end


% analysis parameters:
n_chan = 8; % # of channels, evenly spaced around the screen
chan_centers = linspace(360/n_chan,360,n_chan);

% evaluate basis set at these
angs = linspace(-176,180,90);
if nargin < 4
    which_vox = 0.1; % top 1000 vox
end

if nargin < 5
    trn_tpts = [7:16]; % what we use to train model!
end

align_to = {'targ_ang_all','dist_ang_all'};

func_suffix = 'surf';
delay_tpts = -3:26; % 0.8 s TR ---- what we want to reconstruct


% loop over subj, ROIs and load each session, concatenate, and process
for ss = 1:length(subj)
    
    for vv = 1:length(ROIs)
        
        % load TESTING data from each session and concatenate
        data_tst = [];
        
        for sess_idx = 1:length(sess{ss})
            
            fn = sprintf('%s/%s_trialData/%s_%s_%s_%s_trialData.mat',root,tst_dir,subj{ss},sess{ss}{sess_idx},ROIs{vv},func_suffix);
            fprintf('loading TESTING data from %s...\n',fn);
            thisdata = load(fn);
            
            thisdata.sess = sess_idx*ones(size(thisdata.r_all));
            
            data_tst = cat_struct(data_tst,thisdata,{'rf','TR','which_TRs'}); % skip 'rf', these will be the same
            
        end
        
        
        % list training sessions for this ROI:
        
        tmp_fn_trn = dir( sprintf( '%s/%s_trialData/%s_%s*_%s_%s_trialData.mat', trn_root,trn_dir,subj{ss},trn_sess,ROIs{vv},func_suffix) );
        
        % loop over those and load (as above) - TRAINING data
        data_trn = [];
        fn_trn = cell(length(tmp_fn_trn),1);
        for ff = 1:length(tmp_fn_trn)
            fn_trn{ff} = sprintf('%s/%s',tmp_fn_trn(ff).folder,tmp_fn_trn(ff).name);
            fprintf('loading TRAINING data from %s...\n',fn_trn{ff});
            thisdata = load(fn_trn{ff});
            
            thisdata.sess = ff*ones(size(thisdata.r_map));
            
            data_trn = cat_struct(data_trn,thisdata,{'rf','TR','which_TRs'}); % skip 'rf', these will be the same
        end
        
        which_TRs_tst = data_tst.which_TRs;
        which_TRs_trn = data_trn.which_TRs;
        
        
        % because which_TRs doesn't necessarily start at 1...
        delay_idx = find(ismember(which_TRs_tst,delay_tpts));
        %IEM_trn_tpt_idx = find(ismbember(which_TRs,trn_tpts)); % tpts to average over when training IEM
        
        % save out recons rotated to align with target and with distractor
        % (nans for no-distractor trials)
        recons = cell(length(align_to),1);
        
        recons_raw = nan(size(data_tst.c_all,1),length(angs),length(delay_tpts));
        chan_resp = nan(size(data_tst.c_all,1),n_chan,length(delay_tpts));
        
        
        %w_all = cell(n_folds,1);
        %r2_all = nan(n_folds,length(delay_tpts));
        
        
        % loop over CV folds
        %for fold_idx = 1:n_folds
            
            % ~~~~~~~ first, estimate IEM ~~~~~~~~~~
            
            
            
            
            % pick CV sets
%             trn_runs = ones(length(unique(data.r_LORO)),1);
%             
%             
%             trn_runs(fold_idx) = 0;
%             trn_runs = trn_runs==1; % convert to logical
%             
%             
%             
%             tst_runs = trn_runs~=1;
            
            % LEFTOVER: when we used to use non-unitary #'s of runs for
            % trn/tst
            %trn_idx = ismember(data.r_LORO,ru(find(trn_runs)));
            %tst_idx = ismember(data.r_LORO,ru(find(tst_runs)));
            
            %tst_idx = data_tst.r_LORO==ru(fold_idx);
            %trn_idx = ~tst_idx;
            
            
            
            
            % voxel selection (from _genModelDecode2.m scripts)
            
            
            % first - filter voxels based on whether there's signal (if no signal,
            % there are 0's the entire dataset it seems)
            trndata = mean(data_trn.dt_mapz(:,:,ismember(which_TRs_trn,trn_tpts)),3);
            mystd = std(trndata,[],1);
            
            % if which_vox < 1, that means we're using the RF VE to
            % constrain voxel choice (which will be the same for all CV
            % folds) - choose voxels which have data (no nan/0 std dev)
            % and VE >= threshold specified
            
            if which_vox < 1
                these_vox = mystd~=0 & ~isnan(mystd) & data_trn.rf.ve>=which_vox;
                
                % otherwise, we're using the top N voxels sorted by
                % quadrant-wise F-score, or all voxels, whichever is smaller
            else
                
                
                %trndata = trndata(:,mystd~=0 & ~isnan(mystd));
                
                %extra_vox = sum(mystd==0 | isnan(mystd)); % in case we have to remove voxels, add this many to F-socre
                
                % voxel selection (training data only)
                allF = nan(size(trndata,2),1);
                allp = nan(size(trndata,2),1);
                thisG = data_tst.c_map(trn_idx,2);
                for voxidx = 1:size(trndata,2)
                    
                    thisX = trndata(:,voxidx);
                    
                    
                    [p,tab,stats] = anova1(thisX,thisG,'off');
                    allF(voxidx) = tab{2,5};
                    allp(voxidx) = p;
                    clear thisX p tab stats;
                end
                
                % get rid of NaN's, which can in principle happen when
                % zero std dev, etc.
                f_sorted = sort(allF(~isnan(allF)),'descend');
                if which_vox <= length(allF) % handle case of small ROI
                    f_thresh = f_sorted(which_vox);
                else
                    f_thresh = f_sorted(end);
                end
                
                % deal with rare chance that there are some identical F
                % scores, which would allow incorrect # of vox
                these_vox = allF>=f_thresh;
                if sum(these_vox)>which_vox
                    allidx = find(these_vox);
                    these_vox(allidx((which_vox+1):end)) = 0;
                end
                
            end
                        
            trndata = trndata(:,these_vox);
            
            % n_trials x n_vox
            %trn = mean(data.dt_map(:,:,delay_idx),3);
            trn = trndata;
            
            
            X = spDist_pilot_makeX1(data_trn.c_map(:,1),chan_centers);
            
            
            X = X/max(X(:));
            
            
            %% compute channel weights using only training set (wmMapSpace)
            
            fprintf('computing channel weights\n');
            w = X\trn;
            w_all = w;
            
            
            %             % compute predicted channel response given this W and X
            %             pred_resp_tst = X_tst * w_all{fold_idx};
            %             pred_resp_trn = X * w_all{fold_idx};
            %
            %             r2_all_trn(fold_idx,:) = diag(corr(pred_resp_trn,trn)).^2;
            
            
            % ~~~~~~~~ then, reconstruct each time point ~~~~~~~~~~
            
            for tpt_idx = 1:length(delay_tpts)
                fprintf('TPT: %i\n',tpt_idx);
                
                
                %badvox = (mystd==0);
                
                %mydata = trndata;
                tstdata = mean(data_tst.dt_allz(:,:,delay_idx(tpt_idx)),3);
                
                tstdata = tstdata(:,these_vox);
                tst = tstdata;
                
                
                
                % compute R2 - how well do the predicted voxel responses
                % (pred_resp) on this fold predict observed voxel
                % responses?
                %r2_all(fold_idx,tpt_idx) = 1 - (sum((tst(:) - pred_resp_tst(:)).^2) / sum((tst(:)-mean(tst(:))).^2));
                %                r2_all(fold_idx,:,tpt_idx) = diag(corr(tst,pred_resp_tst)).^2;
                
                
                %% use (optimized) design matrix to compute channel responses using testing data
                
                fprintf('computing channel responses\n');
                
                chan_resp(:,:,tpt_idx) = (inv(w*w.')*w*tst.').';
                
                
                %% compute reconstructions (weighted channel profiles)
                % one aligned to T1, another aligned to T2
                
                %             recons = cell(2,1);
                %             recons{1} = nan(size(data.c_task,1),length(angs));
                %             recons{2} = nan(size(data.c_task,1),length(angs));
                %
                
                
                clear tst tstdata;
                
            end % end of testing/reconstruction loop (within each fold)
            
            clear trn trn_idx tst_idx trn_runs tst_runs w trndata allF f_thresh thisG mystd;
            
       % end % end of CV loop
        
       
       %% after all channel resp computed, coregister each trial
       
       % basis set used for 'raw' reconstructions
       myb_orig = build_basis_polar_mat(angs,chan_centers);

       for aa = 1:length(align_to)
           recons{aa} = nan(size(data_tst.c_all,1),length(angs),length(delay_tpts));
           for tt = 1:size(chan_resp,1)
               
               % we have to build a basis set for each trial, each target
               
               % remove the polar angle of the aligned target
               rot_by = data_tst.(align_to{aa})(tt);%atan2d(data.xy_task(tt,2),data.xy_task(tt,1));
               
               if ~isnan(rot_by)
                   this_rfTh = chan_centers-rot_by; % rotate basis
                   myb = build_basis_polar_mat(angs,this_rfTh);
               end
               % myb is length(angs) x n_channels
               % we want to weight each channel (col) by this trials' channel
               % activation
               % (result shoudl be 1 x length(angs))
               
               
               for tpt_idx = 1:length(delay_tpts)
                   % on no-distractor trials, we won't have to rotate by
                   % anything, so those won't have a distractor-aligned
                   % reconstruction
                   if ~isnan(rot_by)
                       recons{aa}(tt,:,tpt_idx) = (myb * chan_resp(tt,:,tpt_idx).').';
                   end
                   
                   % only need to do this once (raw)
                   if aa == 1
                        recons_raw(tt,:,tpt_idx) = (myb_orig * chan_resp(tt,:,tpt_idx).').';
                   end
               end
               

               
               clear rot_by myb;
           end
           
           
       end
       
       % things we want to save
        
        c_all = data_tst.c_all;
        r_all = data_tst.r_all;
        p_all = data_tst.p_all;
        a_all = [data_tst.targ_ang_all data_tst.dist_ang_all];
        sess_all  = data_tst.sess;
        
        
        % save with VE marker when which_vox < 1, otherwise, number of
        % vox
        if which_vox < 1
            fn2s = sprintf('%s/%s_reconstructions/%s_%s_%s_%s_%ichan_VE%03.f_trn%ito%i_recon_thruTime1.mat',root,tst_dir,subj{ss},horzcat(sess{ss}{:}),ROIs{vv},func_suffix,n_chan,100*which_vox,trn_tpts(1),trn_tpts(end));
        else
            fn2s = sprintf('%s/%s_reconstructions/%s_%s_%s_%s_%ichan_%ivox_trn%ito%i_recon_thruTime1.mat'  ,root,tst_dir,subj{ss},horzcat(sess{ss}{:}),ROIs{vv},func_suffix,n_chan,    which_vox,trn_tpts(1),trn_tpts(end));
        end
        fprintf('saving to %s...\n',fn2s);
        
        save(fn2s,'c_all','r_all','p_all','n_chan','delay_tpts','angs','recons','recons_raw','chan_resp','w_all','which_vox','sess_all','these_vox','a_all','fn_trn');
        
        clear data c_all r_all p_all chan_resp w_all recons a_all;
        
        
        
    end
    
    
end
return