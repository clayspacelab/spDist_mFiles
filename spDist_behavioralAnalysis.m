function spDist_behavioralAnalysis(subj,sess) 
% based on spDist_plotEyeData.m
% in the actual task, cyan = distractor, magenta = no distractor 
% dependencies: misc_util, RMAOV1 

if nargin < 1 || isempty(subj)
    subj = {'AY','CC','EK','KD','MR','SF','XL'};
end

if nargin < 2 || isempty(subj)
sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}}; 
end
root = spDist_loadRoot;

%load raw subject data 

scatterplot_1BC = 0; % plot portions of figure that take quite a while? y/n 

WHICH_EXCL = [13 20 21 22]; % use all exclusion criteria
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


rng(spDist_randSeed);
%% organize data
distractor_bins = unique(all_data.s_all.trialinfo(all_data.s_all.trialinfo(:,1)~=1,6));
distractor_spacing = 360/length(distractor_bins);

%initialize storage cells
all_err = cell(3,1);
all_mu = cell(3,1);
all_rt = cell(3,1);


params_of_interest = {'f_sacc'};

% first, error for no-distractor trials
tmp_err = nan(1,length(params_of_interest),2,length(subj)); % condition x param x dimen of param (2) x subj
tmp_mu  = nan(1,length(params_of_interest),2,length(subj)); 
tmp_errpol = nan(1,length(params_of_interest),2,length(subj));
tmp_mupol = nan(1,length(params_of_interest),2,length(subj));

for ss = 1:length(subj)
    thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==1 & all_data.use_trial==1;
    for pp = 1:length(params_of_interest)
        % distractor bin x param x [radial; tangential] x subj
        tmp_err(1,pp,:,ss) = std( all_data.s_all.(params_of_interest{pp})(thisidx,:), [], 1 );
        tmp_mu(1,pp,:,ss)  = mean( all_data.s_all.(params_of_interest{pp})(thisidx,:), 1 );
        % new
        [dtheta_err, drad_err]= cart2pol(all_data.s_all.(params_of_interest{pp})(thisidx,1), all_data.s_all.(params_of_interest{pp})(thisidx,2));
        tmp_errpol(1,pp,:,ss) = std([drad_err rad2deg(dtheta_err)]);
        [dtheta_mu, drad_mu] = cart2pol(all_data.s_all.(params_of_interest{pp})(thisidx,1), all_data.s_all.(params_of_interest{pp})(thisidx,2));
        tmp_mupol(1,pp,:,ss) = mean([drad_mu rad2deg(dtheta_mu)]); %put rad as x FIRST, theta as Y second, in keeping with later indexing
    end
end
% no distractor std 
all_err{1} = tmp_err; % 1st cell of all err and all mu = no distractor
all_mu{1} = tmp_mu; %note: the first index of all_ is NO DIST
all_mupol{1} = tmp_mupol;
all_errpol{1} = tmp_errpol;

clear tmp_err tmp_mu tmp_errpol tmp_mupol;

% set-up separate, stand-alone structure for bias shuffled test 

tmp_mu = nan(1,length(params_of_interest),2,length(subj));
tmp_mupol = nan(1,length(params_of_interest),2,length(subj));
tic
for xx =1:1000
for ss = 1:length(subj)
    
    for bb =1 %placeholder, keep dimens same size as when we collect others 
       
        for pp = 1:length(params_of_interest)
            tmpd =  all_data.s_all.(params_of_interest{pp});
            tmpidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & all_data.s_all.trialinfo(:,6)==0; 
           
            
            orig_y =  tmpd(tmpidx,2);
            shuf_idx = randperm(length(orig_y));
            shuf_orig_y = orig_y(shuf_idx);
            
            n_tflip = round(length(orig_y)/2); %how many trials are there? flip half. make it a 1:number vector 
            flip_y =  -1.*shuf_orig_y(1:n_tflip);
            noflip_y = shuf_orig_y(n_tflip+1:end);
            tmpd(tmpidx,2) = [flip_y;noflip_y];

            [dtheta_mu, drad_mu] = cart2pol(tmpd(tmpidx,1), tmpd(tmpidx,2));
            tmp_mupol(bb,pp,:,ss,xx) = mean([drad_mu rad2deg(dtheta_mu)]); %put rad as x FIRST, theta as Y second, in keeping with later indexing 

            clear dtheta_err drad_err drad_mu dtheta_mu tmpd thisidx orig_y orig_y_flip
        end
    end
end
end
toc


shuf_mupol{1} = tmp_mupol; % this will only be used for bias testing 

% collect and store random zero bin trials 

tmp_errpol = nan(1,length(params_of_interest),2,length(subj));
tmp_mupol = nan(1,length(params_of_interest),2,length(subj));

for ss = 1:length(subj)
    
    for bb =1 % placeholder, keep dimens same size as when we collect others 
       
        for pp = 1:length(params_of_interest)
            tmpd =  all_data.s_all.(params_of_interest{pp});
            tmpidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & all_data.s_all.trialinfo(:,6)==0 & all_data.s_all.trialinfo(:,10) < 0; %negative CW jitter
            orig_y =  tmpd(tmpidx,2);
            orig_y_flip = orig_y*-1;
            tmpd(tmpidx,2) = orig_y_flip;  % here is where the actual y flip is inserted for CW bins, into tmpd ONLY 
            
            thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & all_data.s_all.trialinfo(:,6)==0; %now, collect all jitters.
 
            % distractor bin x param x [radial; tangential] x subj
            tmp_err(bb,pp,:,ss) = std( all_data.s_all.(params_of_interest{pp})(thisidx,:), [], 1 );
            tmp_mu(bb,pp,:,ss) = mean( tmpd(thisidx,:),  1 ); % use idx with all zero bin from flipped data 

            [dtheta_err, drad_err]= cart2pol(all_data.s_all.(params_of_interest{pp})(thisidx,1), all_data.s_all.(params_of_interest{pp})(thisidx,2));
            tmp_errpol(bb,pp,:,ss) = std([drad_err rad2deg(dtheta_err)]);
            [dtheta_mu, drad_mu] = cart2pol(tmpd(thisidx,1), tmpd(thisidx,2));
            tmp_mupol(bb,pp,:,ss) = mean([drad_mu rad2deg(dtheta_mu)]); %put rad as x FIRST, theta as Y second, in keeping with later indexing 

            clear dtheta_err drad_err drad_mu dtheta_mu tmpd thisidx orig_y orig_y_flip
        end
    end
end
%%%%%%%%



all_err{2} = tmp_err;
all_mu{2} = tmp_mu;
all_mupol{2} = tmp_mupol;
all_errpol{2} = tmp_errpol;

clear tmp_err tmp_mu tmp_mupol tmp_errpol;

%put other distractor bins in - we will collect on this basis, and .*-1 for
%like bins
bins = [1 2 3];
tmp_err = nan(length(bins),length(params_of_interest),2,length(subj));
tmp_mu = nan(length(bins),length(params_of_interest),2,length(subj));
tmp_errpol = nan(1,length(params_of_interest),2,length(subj));
tmp_mupol = nan(1,length(params_of_interest),2,length(subj));

for ss = 1:length(subj)
    for bb = 1:length(bins)        
        for pp = 1:length(params_of_interest)
            
            tmpd =  all_data.s_all.(params_of_interest{pp});
            tmpidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & all_data.s_all.trialinfo(:,6)==(bins(bb)*-1);
            orig_y =  tmpd(tmpidx,2);
            orig_y_flip = orig_y*-1;
            tmpd(tmpidx,2) = orig_y_flip;  % the y flip is inserted for CW bins.
            
            thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & (all_data.s_all.trialinfo(:,6)==bins(bb) | all_data.s_all.trialinfo(:,6)==bins(bb)*-1);
            

            % distractor bin x param x [radial; tangential] x subj
            
            
            tmp_err(bb,pp,:,ss) = std( all_data.s_all.(params_of_interest{pp})(thisidx,:), [], 1 );
            tmp_mu(bb,pp,:,ss) = mean( tmpd(thisidx,:),  1 ); %us thisidx, which collects both bin_n & bin_-n, only on the tmpd data which contains flipped y-vals on cw trials
            
            % new
            [dtheta_err, drad_err]= cart2pol(all_data.s_all.(params_of_interest{pp})(thisidx,1), all_data.s_all.(params_of_interest{pp})(thisidx,2));
            tmp_errpol(bb,pp,:,ss) = std([drad_err rad2deg(dtheta_err)]); %we're only going to be using the theta dimen here 
            [dtheta_mu, drad_mu] = cart2pol(tmpd(thisidx,1), tmpd(thisidx,2));
            tmp_mupol(bb,pp,:,ss) = mean([drad_mu rad2deg(dtheta_mu)]); %put rad as x FIRST, theta as Y second, in keeping with later (x,y) indexing 
            
            clear dtheta_err drad_err drad_mu dtheta_mu tmpd tmpd thisidx orig_y orig_y_flip
            
        end
    end
end

all_err{3} = tmp_err;
all_mu{3} = tmp_mu;
all_mupol{3} = tmp_mupol;
all_errpol{3} = tmp_errpol;
clear tmp_err tmp_mu tmp_errpol tmpmupol;

% create a new container for _alldist - for precision only 
tmp_err = nan(1,length(params_of_interest),2,length(subj)); 
tmp_errpol = nan(1,length(params_of_interest),2,length(subj));

for pp = 1:length(params_of_interest)
    for ss = 1:length(subj)
        thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1;
        
        % distractor bin x param x [radial; tangential] x subj
        tmp_err(1,pp,:,ss) = std( all_data.s_all.(params_of_interest{pp})(thisidx,:), [], 1 );
        [dtheta_err, drad_err]= cart2pol(all_data.s_all.(params_of_interest{pp})(thisidx,1), all_data.s_all.(params_of_interest{pp})(thisidx,2));
        tmp_errpol(1,pp,:,ss) = std([drad_err rad2deg(dtheta_err)]);
        
        clear dtheta_err drad_err
    end
end

all_err_alldist{1} = tmp_err;
all_errpol_alldist{1} = tmp_errpol;

clear tmp_err  tmp_errpol 

% first, RT for no-distractor trials

tmp_rt = nan(1,length(params_of_interest),1,length(subj));

params_of_interest = {'i_sacc_rt'};

for ss = 1:length(subj)
    
    thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==1 & all_data.use_trial==1;

    for pp = 1:length(params_of_interest)
        % distractor bin x param x rt x subj
        tmp_rt(1,pp,1,ss)  = nanmean( all_data.s_all.(params_of_interest{pp})(thisidx),1);

    end
    
end

all_rt{1} = tmp_rt;
clear tmp_rt

% RT for near dist 
tmp_rt = nan(1,length(params_of_interest),1,length(subj));

for ss = 1:length(subj)
    
    thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & all_data.s_all.trialinfo(:,6)==0; 

    for pp = 1:length(params_of_interest)
        % distractor bin x param x rt x subj
        tmp_rt(1,pp,1,ss)  = mean( all_data.s_all.(params_of_interest{pp})(thisidx),1); 
    
    end
    
end

all_rt{2} = tmp_rt;
bins = [1 2 3];
tmp_rt = nan(length(bins),length(params_of_interest),1,length(subj)); % # conditions x # params x dimen of param x subj

for ss = 1:length(subj)
    for bb = 1:length(bins)
        
        thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1 & (all_data.s_all.trialinfo(:,6)== bins(bb) | all_data.s_all.trialinfo(:,6)== bins(bb)*-1);
        
        for pp = 1:length(params_of_interest)
            % distractor bin x param x [radial; tangential] x subj
            tmp_rt(bb,pp,:,ss) = mean( all_data.s_all.(params_of_interest{pp})(thisidx,:),  1 );
            
        end
    end
end

all_rt{3} = tmp_rt;

% RT for _alldist

tmp_rt = nan(1,length(params_of_interest),1,length(subj));
params_of_interest = {'i_sacc_rt'};
param_str = {'RT'};

for ss = 1:length(subj)
    
    thisidx = all_subj==ss & all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1;

    for pp = 1:length(params_of_interest)
        % all distractor x param x rt x subj
        tmp_rt(1,pp,1,ss)  = mean( all_data.s_all.(params_of_interest{pp})(thisidx),1);

    end
    
end

all_rt_alldist{1} =tmp_rt; 

if scatterplot_1BC ==1 
%% Figure 1B : Example participant eye-trace 

dist_colors = [0.7100 0.2128 0.4772; 0 0 1;]; %1 is red, no distractor, 2 is blue (near distractor), 3 is green (far distractor)
figure;
hold on;
xdat_to_plot = [1 2 3 4 5 6 7 8 9 10]; % these are the annotations from our eye-data output, marking the distinct epochs of the task
dist_cond =[1 2];

for cc = 1:length(dist_cond)
    thisidx = find(all_data.use_trial==1 & all_data.subj_all ==2 & all_data.s_all.trialinfo(:,1) == dist_cond(cc)); %subj should be 'CC'
    
    for tt = 1:length(thisidx)
        figure(1)
        %transform the cartesian X,Y coords we have for the trace and put
        %it in polar
        [tmpth,tmpr] = cart2pol(all_data.s_all.X{thisidx(tt)},all_data.s_all.Y{thisidx(tt)});
        
        % change th, keeping r the same, based on th of all_data.s_all.targ
        [adjth,~] = cart2pol(all_data.s_all.targ(thisidx(tt),1),all_data.s_all.targ(thisidx(tt),2));
        
        [aligned_x,~] = pol2cart(tmpth-adjth,tmpr); % we're not using y here, just x
        
        aligned_x = aligned_x(ismember(all_data.s_all.XDAT{thisidx(tt)},xdat_to_plot));
        
        if  cc ==1 %no distractor
            
            this_t = (1:length(aligned_x))/500; % time in seconds (how many samples were taken / frequency of recording device)
            subplot(2,1,1)
            plot(this_t,aligned_x,'-','LineWidth',0.3,'Color', dist_colors(cc,:));
            ylim([-18 18]);
            xlim([0 15])
            xticks([0 1.5 6 13.5])
            xticklabels({'-1.5','0','4.5','12'})
            set(gca,'YTick',[-12 0 12],'TickDir','out');
            hold on;
            set(gcf,'Renderer','painters')
        else
            aligned_x = aligned_x.*-1; % for visualization purposes, flip the trace of the distractor condition to be mirrored over the x-axis
            this_t = (1:length(aligned_x))/500; % time in seconds (how many samples were taken / frequency of recording device)
            subplot(2,1,1)
            plot(this_t,aligned_x,'-','LineWidth',0.3,'Color',[dist_colors(2,:)  0.5000]) %dist_colors(cc,:));
            ylim([-18 18]);
            xlim([0 15])
            xticks([0 1.5 6 13.5])
            xticklabels({'-1.5','0','4.5','12'})
            set(gca,'YTick',[-12 0 12],'TickDir','out');
            hold on;
            set(gcf,'Renderer','painters')
            
        end
        
        if tt ==1
            %for this subject, this trial, this epoch, turn XDAT marker on
            epoch_1_condc = all_data.s_all.XDAT{thisidx(tt)} ==1 ; %  distractor cue, store the entire epoch in a variable
            epoch_2_targ = all_data.s_all.XDAT{thisidx(tt)} ==2; % targets
            epoch_3_d1 = all_data.s_all.XDAT{thisidx(tt)}==3; % delay 1
            epoch_4_dist = all_data.s_all.XDAT{thisidx(tt)} ==4; % distractor
            epoch_5_d2 = all_data.s_all.XDAT{thisidx(tt)} ==5; % delay 2
            epoch_6_go = all_data.s_all.XDAT{thisidx(tt)}==6; % response
            epoch_7_fb = all_data.s_all.XDAT{thisidx(tt)}==7; %feedback
            epoch_8_iti = all_data.s_all.XDAT{thisidx(tt)}==8; % waiting for subj to respond
            
            %epoch 9 = distractor feed-back (no eye-movement)
            %epoch 10 = ITI
            subplot(2,1,2)
            hold on;
            plot(this_t,epoch_1_condc+6,'k-','LineWidth',3);
            plot(this_t,epoch_2_targ+4,'k-','LineWidth',3);
            plot(this_t,epoch_4_dist+2,'k-','LineWidth',3);
            plot(this_t,epoch_6_go,'k-','LineWidth',3);
            xlim([0 15])
            xticks([0 1.5 6 13.5])
            xticklabels({'-1.5','0','4.5','12'})
            yticks([1 3 5 7])
            yticklabels({'Response','Distractor', 'Target','Condition Cue'})
            xlabel('Time relative to delay onset (s)');
            ylabel('Trial Epoch')
            set(gcf,'Renderer','painters')
        else
        end
        clear aligned_x this_t tmpth tmpr
    end
end

set(gcf,'Renderer','painters')
set(gcf,'position',[  1000         626        1147         712])
xlabel('Time relative to delay onset (s)');
ylabel('Eye position (DVA)');

%% Figure 1C : Aligned Saccadic Endpoints; Distractor Absent, Distractor Present 

dist_colors = spDist_condColors; % 1= RED = no distractor; 2 =blue = DISTRACTOR

figure;
params_of_interest= {'f_sacc'}; % what fields do we want to plot?

dist_bins = [1 2]; %1 = no dist, 2 = dist 

for pp = 1:length(params_of_interest)
    
    for dd =1:length(dist_bins)
        
        if dist_bins(dd) == 1
            thisidx = all_data.s_all.trialinfo(:,1)==1 & all_data.use_trial==1;
            x = all_data.s_all.(params_of_interest{pp})(thisidx,1);
            y = all_data.s_all.(params_of_interest{pp})(thisidx,2);
            y_save = nanmean(y);
            
        elseif dist_bins(dd) == 2
            
            thisidx = all_data.s_all.trialinfo(:,1)==2 & all_data.use_trial==1;
            
            x = all_data.s_all.(params_of_interest{pp})(thisidx,1);
            y = all_data.s_all.(params_of_interest{pp})(thisidx,2);
            y_save = nanmean(y);
       
        end
        
        if dd ==1
            subplot(1,2,1)
            title('Distractor Absent')
        elseif dd==2
            subplot(1,2,2)
             title('Distractor Present')
        end
        hold on;
        scatter(x,y,30,dist_colors(dd,:),'filled','MarkerFaceAlpha',.2)
        scatter(mean(x),mean(y),30,'k','filled')
        hold on;
        plot([3 15], [0 0],'--','linewidth',1,'color',[.1 .1 .1])
        clear tmpy
        clear thisidx
        ylim([-6 6])
        xlim([3 15])
        xticks([3:1:15])
        yticks([-6:1:6])
        xticklabels([3:15])
        yticklabels([-6:6])
        axis equal
    end
end

set(gcf,'Renderer','painters')
match_ylim(get(gcf,'Children'))
else
%% Figure 1D : Precision, Distractor Absent vs. Distractor Present (all distractor bins collapsed), polar deg 

% directly compare avg sd for no distractor, near distractor, far
% distractor 
% trials, final sacc only ; averaged over radial/tang...

params_of_interest = {'f_sacc'}; 
param_str = {'final sacc'};
cond_colors = spDist_condColors;  %red, blue
cond_str = {'Distractor Absent','Distractor Present'};

figure;

for pp = 1:length(params_of_interest) 
    
    subplot(1,length(params_of_interest),pp); hold on;
    
    % distractor cond x subj
    thise = nan(2,length(subj)); %collect this error
    for ii = 1:2
     if ii ==1
        thise(ii,:) = squeeze(all_errpol{1}(:,pp,2,:)); % DO NOT mean over 3rd dimension here, look only at the 'y' component 
     elseif ii ==2
        thise(ii,:) = squeeze(all_errpol_alldist{1}(:,pp,2,:)); % this "all_err_alldist" container is used bc all_err separates near vs. far dists. Simplified w this {} where all dist trials are within a single cell 
     end 
    end
 
    plot(1:size(thise,1),thise,'-','Color',[0.3 0.3 0.3],'LineWidth',.5);
    
    for ii = 1:2
        hold on;
        tmpe = std(thise(ii,:))/sqrt(length(subj));
        plot(ii*[1 1],mean(thise(ii,:))+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(ii,:));
        plot(ii,mean(thise(ii,:)),'o','LineWidth',1.5,'Color',cond_colors(ii,:),'MarkerSize',10,'MarkerFaceColor',cond_colors(ii,:));
        
        clear tmpe
    end
    
    xlim([0 3]);
    
    set(gca,'XTick',1:2,'XTickLabel',cond_str','XTickLabelRotation',45,'TickDir','out');
    xlabel('Condition');
    if pp == 1
        ylabel('Error, Polar Angle (\circ)'); % 
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
    
    [h,p,ci,stats] = ttest(thise(1,:)',thise(2,:)') % compare no dist, dist 
    text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',pval),'color','k','fontsize',15) %filled

%ylim([0 1.5])
end
set(gcf,'position',[ 549   724   499   571])
match_ylim(get(gcf,'Children'));

%% RT : No distractor, all distractor bins collapsed (seconds) 


params_of_interest = {'i_sacc_rt'}; %we're using the final saccade
param_str = {'RT'};
cond_colors = spDist_condColors;
cond_str = {'Distractor Absent','Distractor Present'};

figure;

for pp = 1:length(params_of_interest)
    
    subplot(1,length(params_of_interest),pp); hold on;
    
    % distractor cond x subj
    thisrt = nan(2,length(subj)); %collect this error
    for ii = 1:2
     if ii ==1
        thisrt(ii,:) = squeeze(mean(all_rt{1}(:,pp,:,:),3)); % mean over radial/tangential; distractor bin; for n_bins, then mean over these bin
     elseif ii ==2
        thisrt(ii,:) = squeeze(mean(all_rt_alldist{1}(:,pp,:,:),3)); % mean over radial/tangential; distractor bin; for n_bins, then mean over these bin
     end 
    end

    plot(1:size(thisrt,1),thisrt,'-','Color',[0.3 0.3 0.3],'LineWidth',.5);
    
    for ii = 1:2
        hold on;
        tmpe = std(thisrt(ii,:))/sqrt(length(subj));
        plot(ii*[1 1],mean(thisrt(ii,:))+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(ii,:));
        plot(ii,mean(thisrt(ii,:)),'o','LineWidth',1.5,'Color',cond_colors(ii,:),'MarkerSize',10,'MarkerFaceColor',cond_colors(ii,:));
    end
    
    xlim([0 3]);
    
    set(gca,'XTick',1:2,'XTickLabel',cond_str','XTickLabelRotation',45,'TickDir','out');
    xlabel('Condition');
    if pp == 1
        ylabel('RT (seconds)');
    end
    title(param_str{pp});
    
    the_y = [thisrt(1,:)';thisrt(2,:)';];
    the_iv =[ones(length(thisrt(1,:)'),1); 2*ones(length(thisrt(1,:)'),1);];
    subj = [1 2 3 4 5 6 7]';
    the_subj = [subj;subj;];
    x = [the_y the_iv the_subj];
    [fval,pval] = RMAOV1_gh(x,0.05) %one-way RM ANOVA
    RMAOV1(x,0.05) 
    
    [h_rt,p_rt,ci_rt,stats_rt] = ttest(thisrt(1,:)',thisrt(2,:)')
    text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',pval),'color','k','fontsize',15) %filled


end
set(gcf,'position',[ 549   724   499   571])
match_ylim(get(gcf,'Children'));

 %% Figure 1F : Near Distractor trials: Bias, polar degrees 
 
params_of_interest = {'f_sacc'};
param_str = {'final sacc'};

figure;

for pp = 1:length(params_of_interest)
    
    subplot(1,length(params_of_interest),pp); hold on;
    
    % nbins x nsubj
    thism_near = squeeze(all_mupol{2}(:,pp,2,:))'; %this is the only condition that is information for this analysis.
    plot(1+0.15,thism_near,'o','MarkerSize',5,'Markerfacecolor',[0.3 0.3 0.3], 'Color',[0.3 0.3 0.3]);
    plot(1,mean(thism_near,2),'o','MarkerSize',15,'Color',cond_colors(2,:),'MarkerFaceColor',cond_colors(2,:));
    hold on;
    
    tmpe = std(thism_near)/sqrt(length(subj));
    plot([1 1], mean(thism_near)+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(2,:));

    %errorbar(1, mean(thism_near), [std(thism_near)/sqrt(length(subj))],'linewidth',2,'Color',cond_colors(2,:))
    title(param_str{pp});
    
   
    ylabel('Bias, Polar Angle (\circ), toward distractor');
end
xlabel('Near Distractor')
xlim([0 2])

[h_bias_near p_bias_near,CI,STATS] = ttest(thism_near(1,:)') %t-test to determine if bias is different than zero
realT = STATS.tstat;
text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',p_bias_near),'color','k','fontsize',15) %filled

match_ylim(get(gcf,'Children'));
set(gcf,'position',[ 549   724   499   571])

%%%%
% permutation bias testing ; take randomly flipped bias (mu), subject to
% ttest, iter x 
params_of_interest = {'f_sacc'};
param_str = {'final sacc'};

figure;
for xx=1:1000
for pp = 1:length(params_of_interest)
    
    subplot(1,length(params_of_interest),pp); hold on;
    
    % nbins x nsubj
    thism_near = squeeze(shuf_mupol{1}(:,pp,2,:,xx))'; %this is the only condition that is information for this analysis.
    plot(xx+0.15,thism_near,'o','MarkerSize',5,'Markerfacecolor',[0.3 0.3 0.3], 'Color',[0.3 0.3 0.3]);
    plot(xx,mean(thism_near,2),'o','MarkerSize',15,'Color',cond_colors(2,:),'MarkerFaceColor',cond_colors(2,:));
    hold on;
    
    tmpe = std(thism_near)/sqrt(length(subj));
    plot([xx xx], mean(thism_near)+tmpe*[-1 1],'-','LineWidth',1.5,'Color',cond_colors(2,:));

    %errorbar(1, mean(thism_near), [std(thism_near)/sqrt(length(subj))],'linewidth',2,'Color',cond_colors(2,:))
    title('Shuffled Bias');
    
   [h_bias p_bias,CI,STATS] = ttest(thism_near(1,:)');
   shufT(xx) = STATS.tstat;
   
    ylabel('Bias, Polar Angle (\circ), toward distractor');
end
end
xlabel('Near Distractor')
xlim([0 xx+1])
perm_pval = 2 * min(mean(shufT<=realT), mean(shufT>=realT));

text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',perm_pval),'color','k','fontsize',15) %filled

match_ylim(get(gcf,'Children'));
set(gcf,'position',[ 549   724   499   571])

%%%%
%%%%%%%%%%%%%%%%%% SUPPLEMENATARY FIGURE 1 A  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% REPORTED STATS ONLY - PRECISION, Distractor Absent, Distractor Present Abs(Bin)
% STATS ARE ONLY ON Distractor Present (4-level RM ANOVA) : polar deg
% updated 110620 with only y-component 
params_of_interest = {'f_sacc'}; %we're using the final saccade
param_str = {'final sacc'};
cond_colors = [0.7100    0.2128    0.4772;0 0 1;0 0 1;0 0 1;0 0 1;];
cond_str = {'No distractor','Bin 0 Dist','Bin 1 Dist','Bin 2 Dist','Bin 3 Dist'};
figure;

for pp = 1:length(params_of_interest)
    
    subplot(1,length(params_of_interest),pp); hold on;
    
    % distractor cond x subj
    thise = nan(5,length(subj)); %collect this error
    
    for ii = 1:length(all_errpol)
       if ii == 1 || ii ==2
           thise(ii,:) = squeeze(all_errpol{ii}(:,pp,2,:))'; % mean over radial/tangential; distractor bin; for n_bins, then squeeze over these bins to make 1x n_subj vect
       else
           thise([ii:5],:)= squeeze(all_errpol{ii}(:,pp,2,:)) ; % mean over radial/tangential; distractor bin; for n_bins,  then squeeze over these bins to make 1x n_subj vect
           
       end
    end
    
    plot(1:size(thise,1),thise,'-','Color',[0.3 0.3 0.3],'LineWidth',.5);
    
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
        ylabel('Error, Polar Angle (\circ)');
    end
    title(param_str{pp});
    

    the_y = [thise(2,:)'; thise(3,:)'; thise(4,:)'; thise(5,:)';]; % do not use no distractor condition here 
    the_iv =[ones(length(thise(1,:)'),1); 2*ones(length(thise(1,:)'),1); 3*ones(length(thise(1,:)'),1); 4*ones(length(thise(1,:)'),1);];
    subj = [1 2 3 4 5 6 7]';
    the_subj =[subj;subj;subj;subj;];
    x = [the_y the_iv the_subj];
    RMAOV1(x,0.05) %one-way RM ANOVA
    [fval,pval] = RMAOV1_gh(x,0.05) %one-way RM ANOVA
    RMAOV1(x,0.05)


    text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',pval),'color','k','fontsize',15) %filled


end

match_ylim(get(gcf,'Children'));
set(gcf,'position',[ 549   724   499   571])


%%%%%%%%%%%%%%%%%% SUPPLEMENATARY FIGURE 1 B  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% REPORTED STATS ONLY RT :Distractor Absent, Distractor Present abs(bins)

params_of_interest = {'i_sacc_rt'}; %we're using the final saccade
param_str = {'RT'};
cond_colors = [0.7100    0.2128    0.4772;0 0 1;0 0 1;0 0 1;0 0 1;];
cond_str = {'No distractor','Bin 0 Dist','Bin 1 Dist','Bin 2 Dist','Bin 3 Dist'};

figure;

for pp = 1:length(params_of_interest)
    
    subplot(1,length(params_of_interest),pp); hold on;
    
    % distractor cond x subj
    thise = nan(5,length(subj)); %collect this error
    for ii = 1:length(all_rt)
       if ii == 1 || ii ==2
           thise(ii,:) = squeeze(all_rt{ii}(:,pp,:,:))'; 
           
       else
           thise([ii:5],:)= squeeze(all_rt{ii}(:,pp,:,:)); 
           
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
        ylabel('RT');
    end
    title(param_str{pp});
    
    the_y = [thise(2,:)'; thise(3,:)'; thise(4,:)'; thise(5,:)';];
    the_iv =[ones(length(thise(1,:)'),1); 2*ones(length(thise(1,:)'),1); 3*ones(length(thise(1,:)'),1); 4*ones(length(thise(1,:)'),1);];
    subj = [1 2 3 4 5 6 7]';
    the_subj =[subj;subj;subj;subj;];
    x = [the_y the_iv the_subj];
    [fval,pval] = RMAOV1_gh(x,0.05) %one-way RM ANOVA
    RMAOV1(x,0.05)
    text(max(xlim)-(.2*max(xlim)),max(ylim)-(.1*max(ylim)),sprintf('p = %.3f',pval),'color','k','fontsize',15) %filled
    set(gcf,'position',[549   724   499   571])

end

end 

fprintf('behave done')
%c/p to illustrator, scale by %60
return 