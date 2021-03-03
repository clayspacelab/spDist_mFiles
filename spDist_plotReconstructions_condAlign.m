% spDist_pilot_plotReconstructions_thruTime.m
% adapted from MGSMap_plotReconstructions_cv_thruTime1.m
%
% for plotting WM reconstructions during trials with/without spatial
% distractor, plotting distractor-aligned reconstructions, and sorting
% trials based on relative distractor position
%
% for plotting cross-validated WM reconstructions from mapping task,
% computed using MGSMap_channelRespAmp* scripts
%
% TODO: extend to compareReconstruction, which can load multiple sets of
% sessions per subj, and will compare across sessions (and/or across sets
% of tpts, etc... - only one set of comparisons at a time?)
%
% TODO: automate colorlims across figures

function spDist_plotReconstructions_condAlign(subj,sess,ROIs)

root = spDist_loadRoot;

task_dir = 'spDist';

if nargin < 1 || isempty(subj)
   subj = {'AY','CC','EK','KD','MR','SF','XL'}; %geh updated subj to alpha
 
end

if nargin < 2 || isempty(sess)
    % each subj gets one cell, with strings for each sess

    sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
end

if nargin < 3 || isempty(ROIs)
  ROIs = {'V1V2V3','V3AB','hV4','LO1','IPS0IPS1','IPS2IPS3','sPCS'}; %update ROIs 103020

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



%% load data
startidx = 1;
for ss = 1:length(subj)
    
    for vv = 1:length(ROIs)
        
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
                
                all_fidelity = nan(nblankt,size(data.recons{1},3),length(data.recons)); % timecoruse of fidelity for each alignment condition
                all_fidelity_nodist = nan(nblankt,size(data.recons_nodist,3));
                
                all_subj = nan(nblankt,1);
                all_ROIs = nan(nblankt,1);
                all_sess = nan(nblankt,1);
                
                
                angs = data.angs;
                tpts = data.delay_tpts;
                
          
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
            
      
    end
    
end


%% which tpts are we plotting throughout?
% TCS: lte to lessthan
tpts_to_plot = (tpts*myTR) >= t_range_to_plot(1) & (tpts*myTR) < t_range_to_plot(2);

delay_tpts = cell(size(delay_tpt_range,1),1);
for dd = 1:size(delay_tpt_range,1)
    delay_tpts{dd} = (tpts*myTR) >= delay_tpt_range(dd,1) & (tpts*myTR) < delay_tpt_range(dd,2);
end


%% plot each condition (target-locked) (average over subj)
% length(cu) x n_rois

% store all axes across figures
recon_ax = [];


cu = 1;
cond_str = {'No distractor trials'};

figure('name','Figure4A');
for cc = 1:length(cu)
    for vv = 1:length(ROIs)
        
        subplot(length(cu),length(ROIs),(cc-1)*length(ROIs)+vv);hold on;
        
        
        thisd = nan(size(all_recons{1},3),size(all_recons{1},2),length(subj));
        for ss = 1:length(subj)
            thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==cu(cc);
            thisd(:,:,ss) = squeeze(mean(all_recons{1}(thisidx,:,:),1)).';
        end
        % TCS: center each row over middle of TR
        imagesc(angs,tpts(tpts_to_plot)*myTR + myTR/2,mean(thisd(tpts_to_plot,:,:),3));
        colormap viridis
        if cc == 1
            title(ROIs{vv});
        end
        axis ij tight
        set(gca,'XTick',-180:90:180);
        if vv == 1
            xlabel('Polar angle (\circ)');
            ylabel(sprintf('%s - time (s)',cond_str{cc}));
            set(gca,'XTickLabel',{'-180','','0','','180'});
            
        else
            set(gca,'YTick',[],'XTickLabel',[],'YTickLabel',[]);
        end
        xlim([-180 180]);
        caxis([-1.4808 1.7960])
    end
end

set(get(gcf,'Children'),'TickDir','out','Box','off','TickLength',[0.015 0.015],'YTick',[0 4.5 12]);
set(gcf,'Position',[102         405        1353         174]);
recon_ax = [recon_ax; get(gcf,'Children')];


%% plot each condition (target-locked) (average over subj)
% length(cu) x n_rois

cu = 2;
cond_str = {'Distractor trials'};

figure('name','Figure4B');
for cc = 1:length(cu)
    for vv = 1:length(ROIs)
        
        subplot(length(cu),length(ROIs),(cc-1)*length(ROIs)+vv);hold on;
        
        
        thisd = nan(size(all_recons{1},3),size(all_recons{1},2),length(subj));
        for ss = 1:length(subj)
            thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==cu(cc);
            thisd(:,:,ss) = squeeze(mean(all_recons{1}(thisidx,:,:),1)).';
        end
        % TCS: center each row over middle of TR
        imagesc(angs,tpts(tpts_to_plot)*myTR + myTR/2,mean(thisd(tpts_to_plot,:,:),3));
        colormap viridis
        if cc == 1
            title(ROIs{vv});
        end
        axis ij tight
        set(gca,'XTick',-180:90:180);
        if vv == 1
            xlabel('Polar angle (\circ)');
            ylabel(sprintf('%s - time (s)',cond_str{cc}));
            set(gca,'XTickLabel',{'-180','','0','','180'});
            
        else
            set(gca,'YTick',[],'XTickLabel',[],'YTickLabel',[]);
        end
        xlim([-180 180]);
        caxis([-1.4808 1.7960])
    end
end

set(get(gcf,'Children'),'TickDir','out','Box','off','TickLength',[0.015 0.015],'YTick',[0 4.5 12]);
set(gcf,'Position',[102         405        1353         174]);
recon_ax = [recon_ax;get(gcf,'Children')];

%% plot distractor-locked (for distractor)
cond_str ={'Distractor locked'};

figure('name','Figure4C');
for vv = 1:length(ROIs)
    
    subplot(1,length(ROIs),vv);hold on;
    
    
    thisd = nan(size(all_recons{1},3),size(all_recons{1},2),length(subj));
    for ss = 1:length(subj)
        thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==2;
        thisd(:,:,ss) = squeeze(mean(all_recons{2}(thisidx,:,:),1)).';
    end
    % TCS: as above...
    imagesc(angs,tpts(tpts_to_plot)*myTR + myTR/2,mean(thisd(tpts_to_plot,:,:),3));
    colormap viridis;
    
    
    title(ROIs{vv});
    
    axis ij tight
    set(gca,'XTick',-180:90:180);
    if vv == 1
        xlabel('Polar angle (\circ)');
        ylabel(sprintf('%s - time (s)',cond_str{cc}));
        set(gca,'XTickLabel',{'-180','','0','','180'});
        
    else
        set(gca,'YTick',[],'XTickLabel',[],'YTickLabel',[]);
    end
    xlim([-180 180]);
    caxis([-1.4808 1.7960])
end

set(get(gcf,'Children'),'TickDir','out','Box','off','TickLength',[0.015 0.015],'YTick',[0 4.5 12]);
set(gcf,'Position',[102         405        1353         174]);
recon_ax = [recon_ax;get(gcf,'Children')];

%% plot distractor-removed target representation (all positions)
cond_str ={'Target, dist removed'};
figure;

for vv = 1:length(ROIs)
    
    subplot(1,length(ROIs),vv);hold on;
    
    
    thisd = nan(size(all_recons_nodist,3),size(all_recons_nodist,2),length(subj));
    for ss = 1:length(subj)
        thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==2;
        thisd(:,:,ss) = squeeze(mean(all_recons_nodist(thisidx,:,:),1)).';
    end
    
    % TCS: as above
    imagesc(angs,tpts(tpts_to_plot)*myTR + myTR/2,mean(thisd(tpts_to_plot,:,:),3));
    colormap viridis;
    
    
    title(ROIs{vv});
    
    axis ij tight
    set(gca,'XTick',-180:90:180);
    if vv == 1
        xlabel('Polar angle (\circ)');
        ylabel(sprintf('%s - time (s)',cond_str{cc}));
        set(gca,'XTickLabel',{'-180','','0','','180'});
        
    else
        set(gca,'YTick',[],'XTickLabel',[],'YTickLabel',[]);
    end
    xlim([-180 180]);
    caxis([-1.4808 1.7960])
end

set(get(gcf,'Children'),'TickDir','out','Box','off','TickLength',[0.015 0.015],'YTick',[0 4.5 12]);
set(gcf,'Position',[102         405        1353         174]);

recon_ax = [recon_ax;get(gcf,'Children')];

match_clim(recon_ax);


%% plot target fidelity on distractor-/+ trials and distractor fidelity
% row 1: target fidelity
% row 2: distractor fidelity
% row 3: target fidelity (after removing distractor)

% target: without and with distractor; distractor
fidelity_colors = lines(7); fidelity_colors = fidelity_colors(4:6,:);

t_markers = [0 4.5 12]; % onset of delay, distractor, response
mh1 = nan(length(ROIs),length(t_markers));
mh2 = nan(length(ROIs),length(t_markers));
mh3 = nan(length(ROIs),length(t_markers));

mu_fidelity = nan(length(ROIs),size(all_fidelity,2),4); % ROIs x tpts x targ w/ and w/out distractor; distractor; with-distractor after removing distractor...

figure;
% first, plot the target fidelity
for vv = 1:length(ROIs)
    
    subplot(3,length(ROIs),vv); hold on;
    
    
    
    for cc = 1:length(cu)
        
        thisd = nan(length(subj),size(all_fidelity,2));
        
        for ss = 1:length(subj)
            
            thisidx = all_conds(:,1)==cu(cc) & all_subj==ss & all_ROIs==vv;
            thisd(ss,:) = mean(all_fidelity(thisidx,:,1));
            
        end
        
        mu_fidelity(vv,:,cc) = mean(thisd,1);
        
        thise = std(thisd,[],1)/sqrt(length(subj));
        
        % TCS: updated x axis
        % plot mean
        plot(myTR*tpts + myTR/2,mean(thisd,1),'-','LineWidth',1.5,'Color',fidelity_colors(cc,:));
        
        % plot error bars
        plot((myTR*tpts.*[1;1]).' + myTR/2,(mean(thisd,1)+[-1;1].*thise).','--','LineWidth',1,'Color',fidelity_colors(cc,:));
        
        
        yline(0);
        % TODO: plot std error across subj
        
        title(ROIs{vv});
        if vv == 1
            ylabel('Target fidelity');
        else
            set(gca,'YTickLabel',[]);
        end
        
        set(gca,'XTick',[0:6:24],'TickDir','out','XTickLabel',[]);
        
        clear thisd thise;
    end
    
    mh1(vv,:) = plot(t_markers.*[1;1],[0 .1],'-','Color',[0.7 0.7 0.7],'LineWidth',0.75);
    
    
    % ---------- SECOND ROW ---------------------
    subplot(3,length(ROIs),vv+length(ROIs)); hold on;
    
    
    
    thisd = nan(length(subj),size(all_fidelity,2));
    
    for ss = 1:length(subj)
        
        thisidx = all_conds(:,1)==2 & all_subj==ss & all_ROIs==vv;
        thisd(ss,:) = mean(all_fidelity(thisidx,:,2));
        
    end
    
    mu_fidelity(vv,:,1+length(cu)) = mean(thisd,1);
    
    thise = std(thisd,[],1)/sqrt(length(subj));

    
    
    plot(myTR*tpts + myTR/2,mean(thisd,1),'-','LineWidth',1.5,'Color',fidelity_colors(3,:));
    plot((myTR*tpts.*[1;1]).' + myTR/2,(mean(thisd,1)+[-1;1].*thise).','--','LineWidth',1,'Color',fidelity_colors(3,:));
    
    yline(0);
    
    if vv == 1
        ylabel('Distractor fidelity');
        xlabel('Time (s)');
    else
        set(gca,'YTickLabel',[]);
    end
    
    set(gca,'XTick',[0:6:24],'TickDir','out');
    
    mh2(vv,:) = plot(t_markers.*[1;1],[0 .1],'-','Color',[0.7 0.7 0.7],'LineWidth',0.75);
    
    
    clear thisd;
    % ---------- THIRD ROW ------------ 
    
    subplot(3,length(ROIs),vv+2*length(ROIs)); hold on;
    
    thisd = nan(length(subj),size(all_fidelity,2));
    
    for ss = 1:length(subj)
        
        thisidx = all_conds(:,1)==2 & all_subj==ss & all_ROIs==vv;
        thisd(ss,:) = mean(all_fidelity_nodist(thisidx,:));
        
    end
    
    mu_fidelity(vv,:,2+length(cu)) = mean(thisd,1);
    
    thise = std(thisd,[],1)/sqrt(length(subj));

    
    
    plot(myTR*tpts + myTR/2,mean(thisd,1),'-','LineWidth',1.5,'Color',fidelity_colors(2,:));
    plot((myTR*tpts.*[1;1]).' + myTR/2,(mean(thisd,1)+[-1;1].*thise).','--','LineWidth',1,'Color',fidelity_colors(2,:));
    
    yline(0);
    
    if vv == 1
        ylabel('WM fidelity (minus distractor)');
        xlabel('Time (s)');
    else
        set(gca,'YTickLabel',[]);
    end
    
    set(gca,'XTick',[0:6:24],'TickDir','out');
    
    mh3(vv,:) = plot(t_markers.*[1;1],[0 .1],'-','Color',[0.7 0.7 0.7],'LineWidth',0.75);
    
    
    clear thisd;
    
    
    
end

%myy = cell2mat(get(get(gcf,'Children'),'YLim'));

myy = match_ylim(get(gcf,'Children'));
set(mh1,'YData',[min(myy(:,1)) max(myy(:,2))]);
set(mh2,'YData',[min(myy(:,1)) max(myy(:,2))]);
set(mh3,'YData',[min(myy(:,1)) max(myy(:,2))]);

set(gcf,'Position',[185         745        1843         470]);



return