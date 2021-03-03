% spDist_exampleTrial.m
%
% loads and plots data from an example trial
%
% (for Fig. 2)
%
% good trials:
% CC, V3AB: 30, 237, 297 (with distractor); no distractor: 134, 135

close all;

subj_sess = 'CC_spDist1spDist2';
ROI  = 'V3AB';
datasuff = 'surf_8chan_VE010_trn7to15_recon_thruTime1.mat';


thisdata = load(sprintf('%s/spDist_reconstructions/%s_%s_%s',spDist_loadRoot,subj_sess,ROI,datasuff));

myangs = thisdata.angs+180;

nrows = 8;
ncols = 20;

thisrow = 1;
thiscol = 1;

figure;
for tt = 1:size(thisdata.c_all,1)
    
    
    trialnum = tt;
    subplot(nrows,ncols,(thisrow-1)*ncols+thiscol);
    hold on;
    imagesc(myangs,0.75*thisdata.delay_tpts,squeeze(thisdata.recons_raw(trialnum,:,:)).');
    
    plot(180+thisdata.a_all(trialnum,1),0,'cv','MarkerFaceColor','c');
    plot(180+thisdata.a_all(trialnum,2),4.5,'rv','MarkerFaceColor','r');
    
    xlim([0 360]);
    ylim([-2.25 12]);
    axis ij;
    
    set(gca,'TickDir','out','XTick',[0 180 360],'YTick',[0 5 10]);
    
    colormap viridis;
    
    if thiscol==ncols
        thiscol = 1;
        if thisrow == nrows
            thisrow = 1;
            figure;
        else
            thisrow = thisrow+1;
        end
    else
        thiscol = thiscol+1;
    end
    
end


%% plot specific trials

extrial_dist   = 297;
extrial_nodist = 175;

figure('name','Figure3B');

subplot(1,2,1);
hold on;
imagesc(myangs,0.75*thisdata.delay_tpts,squeeze(thisdata.recons_raw(extrial_nodist,:,:)).');

plot(180+thisdata.a_all(extrial_nodist,1),0,  'cv','MarkerFaceColor','c');

xlim([0 360]);
ylim([-2.25 12]);
axis ij;

set(gca,'TickDir','out','XTick',[0 180 360],'YTick',[0 5 10]);
title('No distractor');
colormap viridis;


subplot(1,2,2);
hold on;
imagesc(myangs,0.75*thisdata.delay_tpts,squeeze(thisdata.recons_raw(extrial_dist,:,:)).');

plot(180+thisdata.a_all(extrial_dist,1),0,  'cv','MarkerFaceColor','c');
plot(180+thisdata.a_all(extrial_dist,2),4.5,'rv','MarkerFaceColor','r');

xlim([0 360]);
ylim([-2.25 12]);
axis ij;

set(gca,'TickDir','out','XTick',[0 180 360],'YTick',[0 5 10],'YTickLabel',[]);
title('Distractor');
colormap viridis;

match_clim(get(gcf,'Children'));