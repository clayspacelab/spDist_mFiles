% spDist_pilot_plotReconstructions_thruTime_gh.m
% adapted from MGSMap_plotReconstructions_cv_thruTime1.m
%
% for plotting WM reconstructions during trials with/without spatial
% distractor, plotting distractor-aligned reconstructions, and sorting
% trials based on relative distractor position

% for plotting cross-validated WM reconstructions from mapping task,
% computed using MGSMap_channelRespAmp* scripts
%
% TODO: extend to compareReconstruction, which can load multiple sets of
% sessions per subj, and will compare across sessions (and/or across sets
% of tpts, etc... - only one set of comparisons at a time?)

function spDist_plotReconstructions_thruTime_Figure4(subj,sess,ROIs)

root = spDist_loadRoot;

task_dir = 'spDist';

if nargin < 1 || isempty(subj)
    subj = {'CC','KD','AY','MR','XL','SF','EK'};
end

if nargin < 2 || isempty(sess)
    % each subj gets one cell, with strings for each sess
    % TODO: automate...
    sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
end

if nargin < 3 || isempty(ROIs)
    %ROIs =
    %{'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'}; ORIG
    ROIs = {'V1','V2','V3','V3AB','IPS0','IPS1','IPS2','IPS3','sPCS'};
    %ROIs = {'V1V2V3','V3AB','IPS0IPS1','IPS2IPS3','sPCS'}
end


func_suffix = 'surf';

cat_mode = 1; % if 1, look for catSess1Ses...SessN_ files, otherwise, just look for each session in turn

nchan = 8;
which_vox = 0.1; % if > 1 , look for WHICH_VOXvox string; otherwise, look for VE<100*WHICH_VOX>

smooth_by = 1; % if this is 1, use regular files, otherwise, load smooth_by files

myTR = 0.75;

t_range_to_plot = [-inf 12]; % plot b/w these (s)

trn_tpts = 7:15; % if blank, load files w/ no _trn%ito%i, otherwise,

plot_tpts = 7:15; % for a plot where we average reconstructions over a fixed time window

dist_time = 4.5; % onset at 4 s


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

delay_tpt_range = [2 5; 5 9;9 10.5; 10.5 12];



% for fidelity timecourses
tmpcolors = lines(7);

ROI_colors = [repmat(tmpcolors(1,:),3,1); % V1, V2, V3
    tmpcolors(4,:);             % V3AB
    tmpcolors(1,:);             % hV4
    repmat(tmpcolors(3,:),1,1); % VO1
    repmat(tmpcolors(1,:),2,1); % LO1/2
    
    repmat(tmpcolors(2,:),2,1); % TO1-2
    
    repmat(tmpcolors(5,:),2,1); % IPS0-1
    repmat(tmpcolors(6,:),2,1); % IPS2-3
    tmpcolors(7,:);             % sPCS
    ];             % iPCS (color 1...)


%% load data
startidx = 1;
for ss = 1:length(subj)
    
    for vv = 1:length(ROIs)
        
        if cat_mode == 1
            % just one file to load
            fn = sprintf('%s/%s_reconstructions/%s_%s_%s_%s_%ichan%s%s%s_recon_thruTime1.mat',root,task_dir,subj{ss},horzcat(sess{ss}{:}),ROIs{vv},func_suffix,nchan,vox_str,smooth_str,trn_str);
            
            fprintf('loading %s...\n',fn);
            data = load(fn);
            
              
            if vv == 1 && ss == 1
                % initialize variables...
                
                
                nblankt = length(ROIs)*size(data.recons{1},1);
                all_recons = cell(size(data.recons));
                for aa = 1:length(data.recons)
                    all_recons{aa} = nan(nblankt,size(data.recons{aa},2),size(data.recons{aa},3));
                end
                
                all_recons_nodist = nan(nblankt,size(data.recons_nodist,2),size(data.recons_nodist,3));
                
                all_conds = nan(nblankt,size(data.c_all,2));
                all_angs = nan(nblankt,size(data.a_all,2));
                
                all_fidelity = nan(nblankt,size(data.recons{1},3),length(data.recons)); % timecourse of fidelity for each alignment condition
                all_fidelity_nodist = nan(nblankt,size(data.recons_nodist,3));
                
                all_subj = nan(nblankt,1);
                all_ROIs = nan(nblankt,1);
                all_sess = nan(nblankt,1);
                
                
                angs = data.angs;
                tpts = data.delay_tpts;
                
                % ugh have to do this in a multi-D array...
                %all_r2 = nan(length(ROIs),length(tpts),length(subj));
                
            end
            
            
            
            thisidx = startidx:(startidx+size(data.c_all,1)-1);
            
            for aa = 1:length(all_recons)
                all_recons{aa}(thisidx,:,:) = data.recons{aa};
                all_fidelity(thisidx,:,aa) = squeeze(mean(cosd(angs) .* data.recons{aa},2));
            end
            
            all_recons_nodist(thisidx,:,:) = data.recons_nodist;
            all_fidelity_nodist(thisidx,:) =  squeeze(mean(cosd(angs) .* data.recons_nodist,2));
            
            all_conds(thisidx,:) = data.c_all;
            all_angs(thisidx,:) = data.a_all;
            
            
            all_subj(thisidx) = ss;
            
            
            all_ROIs(thisidx) = vv;
            
            all_sess(thisidx) = data.sess_all;
            
            
            startidx = thisidx(end)+1;
            
            clear data;
            
        else
            % NOT SUPPORTED YET!!!!
            
            for sess_idx = 1:length(sess{ss})
                % build fn
                fn = sprintf('%swmChoose_reconstructions/%s_%s_%s_%s_%ichan%s%s%s_recon_cv_thruTime1.mat',root,subj{ss},sess{ss}{sess_idx},ROIs{vv},func_suffix,nchan,vox_str,smooth_str,trn_str);
                
                fprintf('loading %s...\n',fn);
                data = load(fn);
                
                
                if vv == 1 && ss == 1
                    % initialize variables...
                    
                    
                    nblankt = length(ROIs)*numel(sess)*size(data.recons,1);
                    
                    all_recons = nan(nblankt,size(data.recons,2),size(data.recons,3));
                    all_conds = nan(nblankt,size(data.c_map,2));
                    
                    all_fidelity = nan(nblankt,size(data.recons,3)); % timecoruse of fidelity
                    
                    
                    all_subj = nan(nblankt,1);
                    all_ROIs = nan(nblankt,1);
                    all_sess = nan(nblankt,1);
                    
                    angs = data.angs;
                    tpts = data.delay_tpts;
                    
                    all_r2 = nan(length(ROIs),length(tpts),length(subj));
                    
                end
                
                % set up our variable used to compute R2
                if sess_idx == 1
                    tmp_r2 = nan(length(tpts),length(sess{ss})); % average acorss sessions...
                end
                
                thisidx = startidx:(startidx+size(data.c_map,1)-1);
                
                
                all_recons(thisidx,:,:) = data.recons;
                all_fidelity(thisidx,:) = squeeze(mean(cosd(angs) .* data.recons,2));
                
                all_conds(thisidx,:) = data.c_map;
                
                
                
                all_subj(thisidx) = ss;
                
                
                all_ROIs(thisidx) = vv;
                
                all_sess(thisidx) = sess_idx;
                
                tmp_r2(:,sess_idx) = squeeze(mean(mean(data.r2_all,1),2));
                
                startidx = thisidx(end)+1;
                
                clear data;
                
            end
            
            %all_r2(vv,:,ss) = mean(tmp_r2,2);
            %clear tmp_r2;
            
            
        end
    end
    
end


%% which tpts are we plotting throughout?
tpts_to_plot = (tpts*myTR) >= t_range_to_plot(1) & (tpts*myTR) <= t_range_to_plot(2);

delay_tpts = cell(size(delay_tpt_range,1),1);
for dd = 1:size(delay_tpt_range,1)
    delay_tpts{dd} = (tpts*myTR) >= delay_tpt_range(dd,1) & (tpts*myTR) < delay_tpt_range(dd,2);
end

% last tpt of delay_tpt{4} is at ind 19, and last ind of tpts_to_plot is at
% 20. delay_tpts{3} = 16,17 ; delay_tpts{4} =18,19 % this is bc for this
% delay tpts plot, TRs doesnt exactly align with chosen subplot timing 


%% align like distractor bins (and flip/average) 1D FOR ALL BINS, first near, then far
% goal here is to align cw/ccw distractor bins and flip one set to match
% - for 0-bin, need to determine which trials are CW/CCW and flip
%   accordingly
% new code 
% flip negative angles

roi_str = {'V1','V2','V3','V3AB','IPS0','IPS1','IPS2','IPS3','sPCS'};
cond_colors =lines(length(roi_str));
tpts_to_plot = (tpts*myTR) >= t_range_to_plot(1) & (tpts*myTR) <= t_range_to_plot(2);
% look for all trials where <> is < 0, flipLR the reconstruction

tmprel =  all_angs(:,2) - all_angs(:,1);
this_rel = mod((tmprel+180), 360)-180;
% sign of this matches all_conds(:,6) (relative distractor angle bin)

flipidx = this_rel<0;

all_recons_flipped = all_recons{1}; 

all_recons_flipped(flipidx,:) = fliplr(all_recons_flipped(flipidx,:));

% collect the last 4 TRs of the delay period 
tpts_to_plot(1:16)=0; %this aligns with 12 sec end from tpt_to_plot derivation, but doesnt with delay_tpts. 
to_plot =[1 2];
store_b = nan(length(to_plot),length(ROIs),length(subj));

figure
for ff = 1 %:length(to_plot)
    for vv = 1:length(ROIs)
        subplot(length(to_plot),length(ROIs),(ff-1)*length(ROIs)+vv);hold on;
        bias = nan(size(all_recons{1},3),size(all_recons_flipped,2),length(subj));
        thisd = nan(size(all_recons{1},3),size(all_recons_flipped,2),length(subj));
        thisb = nan(1,length(subj));
        for ss = 1:length(subj)
            if ff ==1
            thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==2 & all_conds(:,6)==0; 
            thisd(:,:,ss) = squeeze(mean(all_recons_flipped(thisidx,:,:),1)).';
            bias(:,:,ss)  = atan2d(squeeze(mean(all_recons_flipped(thisidx,:,:),1)).'.*sind(angs),squeeze(mean(all_recons_flipped(thisidx,:,:),1)).'.*cosd(angs));
            thisb(ss) = atan2d(sum(mean(mean(all_recons_flipped(thisidx,:,tpts_to_plot),1),3).*sind(angs)),sum(mean(mean(all_recons_flipped(thisidx,:,tpts_to_plot),1),3).*cosd(angs)));
            store_b(ff,vv,ss) = thisb(ss);
            else
            thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==2 & all_conds(:,6)~=0; 
            thisd(:,:,ss) = squeeze(mean(all_recons_flipped(thisidx,:,:),1)).';
            bias(:,:,ss)  = atan2d(squeeze(mean(all_recons_flipped(thisidx,:,:),1)).'.*sind(angs),squeeze(mean(all_recons_flipped(thisidx,:,:),1)).'.*cosd(angs));
            thisb(ss) = atan2d(sum(mean(mean(all_recons_flipped(thisidx,:,tpts_to_plot),1),3).*sind(angs)),sum(mean(mean(all_recons_flipped(thisidx,:,tpts_to_plot),1),3).*cosd(angs)));
            store_b(ff,vv,ss) = thisb(ss);
            end
        end
       hold on
 
    t = mean(thisd((tpts_to_plot),:,:),1);
    my_sem = nanstd(t,[],3) /sqrt(length(subj));
    plot(linspace(-180,180,90).*[1 1]' ,mean(mean(thisd(tpts_to_plot,:,:),3)),'LineWidth',.5,'color',cond_colors(1,:))
    hold on;
    plot(linspace(-180,180,90),mean(mean(thisd(tpts_to_plot,:,:),3))+1.*my_sem,'--','LineWidth',.5,'color',cond_colors(1,:))
    hold on;
    plot(linspace(-180,180,90),mean(mean(thisd(tpts_to_plot,:,:),3))-1.*my_sem,'--','LineWidth',.5,'color',cond_colors(1,:))

    plot([min(xlim) max(xlim)], [0 0], '--', 'color', [.2 .2 .2])
    x = linspace(-180,180,90);
    [y_max, y_max_ind] = max(mean(mean(thisd(tpts_to_plot,:,:),3)));
    x_ofymax = x(y_max_ind);
    line([x_ofymax x_ofymax], [0 2], 'color','k','linewidth',.5,'linestyle','-')
    line([0 0], [-.8 max(ylim)], 'color',[.2 .2 .2],'linewidth',0.5,'linestyle','--')

    if ff == 1
        title(ROIs{vv});
    end
   
    set(gca,'XTick',-180:40:180,'XTickLabelRotation',45,'TickDir','out')
    
    xlabel('Polar angle (\circ)');
    
    if ff ==1 
    ylabel('Near Distractor');
    else
    ylabel('Far Distractor');
    end
    
    ylim([-.8 2])
    
    set(gca,'XTickLabel',{'-180','-140','-100','-60','-20','20','60','100','140','180'});
    set(gca,'YTickLabel')
    xlim([-180 180]);
   

    end
clear thisidx
clear thisb
clear thisd
end


%store_b dimens
% 1: near or far distractor
% 2: ROI
% 3: subj

for tt =1:size(store_b,1)
for rr =1:length(ROIs)
[H(tt,rr) P(tt,rr)] =  ttest(store_b(tt,rr,:));

end 


end

sprintf('P-values = %s',P) 
