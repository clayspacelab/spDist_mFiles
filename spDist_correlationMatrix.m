function spDist_correlationMatrix(subj,sess,ROIs)


root = '/share/data/spDist/';

task_dir = 'spDist';

if nargin < 1 || isempty(subj)
    subj = {'CC','AY','MR'};
       %subj = {'CC','AY','MR','EK','KD','AY','XL'};
end

if nargin < 2 || isempty(sess)
    % each subj gets one cell, with strings for each sess
    % TODO: automate...
    %sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
    sess = {{'spDistLong1','spDistLong2'},{'spDistLong1','spDistLong2'},{'spDistLong1','spDistLong2'}};
    %sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}}
end

if nargin < 3 || isempty(ROIs)
    ROIs = {'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'};
end


func_suffix = 'surf';

cat_mode = 1; % if 1, look for catSess1Ses...SessN_ files, otherwise, just look for each session in turn

nchan = 8;
which_vox = 0.1; % if > 1 , look for WHICH_VOXvox string; otherwise, look for VE<100*WHICH_VOX>

smooth_by = 1; % if this is 1, use regular files, otherwise, load smooth_by files

myTR = 0.75;


if ismember(sess{1},{'spDistLong1','spDistLong2'})
    t_range_to_plot = [-inf 16.5];
    delay_tpt_range = [1 4; 6 9; 12 15];
    total_delay_range =[0 16.5];
else
    t_range_to_plot = [-inf 12];
    delay_tpt_range = [1 3;4.5 5.5; 9 12];
    total_delay_range =[0 12];
end


trn_tpts = 7:15; % if blank, load files w/ no _trn%ito%i, otherwise,


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

delay_tpts = cell(size(delay_tpt_range,1),1);
for dd = 1:size(delay_tpt_range,1)
    delay_tpts{dd} = (tpts*myTR) >= delay_tpt_range(dd,1) & (tpts*myTR) < delay_tpt_range(dd,2);
end

total_delay_tpts = cell(size(total_delay_range,1),1);
for dd = 1:size(total_delay_range,1)
    total_delay_tpts{dd} = (tpts*myTR) >= total_delay_range(dd,1) & (tpts*myTR) < total_delay_range(dd,2);
end
%%
%% distractor chopped delay - plot by ROI - avg over epoch, per subject
% CORRELATION MATRIX
% for every subj, epoch, roi correlate the vector of avg fidelity for every
% i-j ROI. store this correlation in the i'th,j'th element of the matrix
epoch ={'pre','during','post'};
cmap = cbrewer('div','RdBu',100);
colormap(cmap)
thiscorr = [];
thisp=[];
for ss =1:length(subj)
    figure(ss);
    thisii =[];
    
    for dd = 1:length(delay_tpts)
       
        for ii= 1:length(ROIs)
           this_ii_idx = all_conds(:,1)==2 & all_subj==ss & all_ROIs==ii;
           thisii(ii,:,dd) = mean(all_fidelity_nodist(this_ii_idx,delay_tpts{dd}),2);
        
 %   subplot(length(subj),length(ROIs),vv+(ss-1)*length(ROIs));hold on;

        subplot(length(delay_tpts),length(ROIs),(dd-1)*length(ROIs)+ii); hold on;
        histogram(thisii(ii,:,dd))
        set(gcf,'Position',[861   349   778   643]);
        title(sprintf(' %s',ROIs{ii}))
        end 
    end
        sgtitle(sprintf('Mean fidelity on each trial, subj %s',subj{ss}))  

end
%% distractor chopped delay - plot by ROI - avg over epoch, per subject
% CORRELATION MATRIX
% for every subj, epoch, roi correlate the vector of avg fidelity for every
% i-j ROI. store this correlation in the i'th,j'th element of the matrix
epoch ={'pre','during','post'};
cmap = cbrewer('div','RdBu',100);
colormap(cmap)
thiscorr = [];
thisp=[];
for ss =1:length(subj)
    figure(ss);
    thisjj =[];
    thisii =[];
    
    for dd = 1:length(delay_tpts)
       
        for ii= 1:length(ROIs)
            this_ii_idx = all_conds(:,1)==2 & all_subj==ss & all_ROIs==ii;
            for jj = 1:length(ROIs)
                this_j_idx = all_conds(:,1)==2 & all_subj==ss & all_ROIs==jj;
                thisjj(jj,:,dd) = mean(all_fidelity_nodist(this_j_idx,delay_tpts{dd}),2)';
                thisii(ii,:,dd) = mean(all_fidelity_nodist(this_ii_idx,delay_tpts{dd}),2)';
                [r, p] = corrcoef(thisjj(jj,:,dd),thisii(ii,:,dd));
                thiscorr(ii,jj,dd,ss) = r(1,2);
                thisp(ii,jj,dd,ss) = p(1,2);
            end
        end
        
        subplot(1,length(delay_tpts),dd); hold on;
        imagesc(thiscorr(:,:,dd,ss),[-1 1])
        colormap(cmap)
        axis ij tight
        axis square
        set(gca,'YTick',1:length(ROIs),'Yticklabel',{'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'});
        set(gca,'XTick',1:length(ROIs),'Xticklabel',{'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'},'XTicklabelrotation',45);
        set(gcf,'Position',[861   349   778   643]);

        title(sprintf(' %s - dist',epoch{dd}))
        colorbar

    end
        sgtitle(sprintf('Distractor-removed target fidelity, avg, subj %s',subj{ss}))  
if ss == length(subj)
        figure
        for as = 1:3
        subplot(1,length(delay_tpts),as); hold on;
        imagesc(mean(thiscorr(:,:,as,:),4),[-1 1])
        colormap(cmap)
        axis ij tight
        axis square
        set(gca,'YTick',1:length(ROIs),'Yticklabel',{'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'});
        set(gca,'XTick',1:length(ROIs),'Xticklabel',{'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'},'XTicklabelrotation',45);
        set(gcf,'Position',[861   349   778   643]);
        title(sprintf(' %s - dist',epoch{as}))
        colorbar

        end
else 
end
end
%  for tt = 1:3
%         figure; subplot(1,length(delay_tpts),tt);
%         hold on; 
%         plot([1:size(thiscorr,1)],thiscorr(:,:,as,tt),'o')
%  end 
 fprintf('end avg sec')
%% use max
% %% distractor chopped delay - plot by ROI - avg over epoch, per subject
% CORRELATION MATRIX
% for every subj, epoch, roi correlate the vector of avg fidelity for every
% i-j ROI. store this correlation in the i'th,j'th element of the matrix
epoch ={'pre','during','post'};
cmap = cbrewer('div','RdBu',100);
colormap(cmap)
thiscorr = [];
for ss =1:length(subj)
    figure;
    thisjj =[];
    thisii =[];
    
    for dd = 1:length(delay_tpts)
       
        for ii= 1:length(ROIs)
            this_ii_idx = all_conds(:,1)==2 & all_subj==ss & all_ROIs==ii;
            [thisii(ii,:,dd), ~] = max(all_fidelity_nodist(this_ii_idx,delay_tpts{dd}),[],2);

            for jj = 1:length(ROIs)
                this_jj_idx = all_conds(:,1)==2 & all_subj==ss & all_ROIs==jj;
                [thisjj(jj,:,dd),~] = max(all_fidelity_nodist(this_jj_idx,delay_tpts{dd}),[],2);
                r = corrcoef(thisjj(jj,:,dd),thisii(ii,:,dd));
                thiscorr(ii,jj,dd,ss) = r(1,2);
                
            end
        end
     
        subplot(1,length(delay_tpts),dd); hold on;
        imagesc(thiscorr(:,:,dd,ss),[-1 1])
        colormap(cmap)
        axis ij tight
        axis square
        set(gca,'YTick',1:length(ROIs),'Yticklabel',{'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'});
        set(gca,'XTick',1:length(ROIs),'Xticklabel',{'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'});
        title(sprintf(' %s - dist',epoch{dd}))
        colorbar


    end
        sgtitle(sprintf('Distractor-removed target fidelity, max, subj %s',subj{ss}))  
if ss == length(subj)
        figure
        for as = 1:length(delay_tpts)
        subplot(1,length(delay_tpts),as); hold on;
        imagesc(mean(thiscorr( :,:,as,:),4),[-1 1])
        colormap(cmap)
        axis ij tight
        axis square
        set(gca,'YTick',1:length(ROIs),'Yticklabel',{'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'});
        set(gca,'XTick',1:length(ROIs),'Xticklabel',{'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'});
        title(sprintf(' %s - dist',epoch{as}))
        colorbar


        end
else 
end
end
%% total delay - plot by ROI - avg over epoch, per subject
% CORRELATION MATRIX
% for every subj, epoch, roi correlate the vector of avg fidelity for every
% i-j ROI. store this correlation in the i'th,j'th element of the matrix
total_epoch ={'total delay'};
cmap = cbrewer('div','RdBu',100);
colormap(cmap)
thiscorr = [];
for ss =1:length(subj)
    figure;
    thisjj =[];
    thisii =[];
    
    for dd = 1:length(total_delay_tpts)
       
        for ii= 1:length(ROIs)
            this_ii_idx = all_conds(:,1)==2 & all_subj==ss & all_ROIs==ii;
            for jj = 1:length(ROIs)
                this_j_idx = all_conds(:,1)==2 & all_subj==ss & all_ROIs==jj;
                thisjj(jj,:,dd) = mean(all_fidelity_nodist(this_j_idx,total_delay_tpts{dd}),2)';
                thisii(ii,:,dd) = mean(all_fidelity_nodist(this_ii_idx,total_delay_tpts{dd}),2)';
                r = corrcoef(thisjj(jj,:,dd),thisii(ii,:,dd));
                thiscorr(ii,jj,dd,ss) = r(1,2);
                
            end
        end
        
        subplot(1,length(total_delay_tpts),dd); hold on;
        imagesc(thiscorr(:,:,dd,ss),[-1 1])
        colormap(cmap)
        axis ij tight
        axis square
        set(gca,'YTick',1:length(ROIs),'Yticklabel',{'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'});
        set(gca,'XTick',1:length(ROIs),'Xticklabel',{'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'});
        title(sprintf(' %s - dist',epoch{dd}))
        colorbar


    end
        sgtitle(sprintf('Distractor-removed target fidelity, avg, subj %s',subj{ss}))  
if ss == length(subj)
        figure
        for as = 1:length(total_delay_tpts)
        subplot(1,length(total_delay_tpts),as); hold on;
        imagesc(mean(thiscorr(:,:,as,:),4),[-1 1])
        colormap(cmap)
        axis ij tight
        axis square
        set(gca,'YTick',1:length(ROIs),'Yticklabel',{'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'});
        set(gca,'XTick',1:length(ROIs),'Xticklabel',{'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'});
        title(sprintf(' %s - dist,all subj',total_epoch{as}))
        colorbar


        end
else 
end
end
%% use max
% %% distractor chopped delay - plot by ROI - avg over epoch, per subject
% CORRELATION MATRIX
% for every subj, epoch, roi correlate the vector of avg fidelity for every
% i-j ROI. store this correlation in the i'th,j'th element of the matrix
total_epoch ={'total delay'};
cmap = cbrewer('div','RdBu',100);
colormap(cmap)
thiscorr = [];
for ss =1:length(subj)
    figure;
    thisjj =[];
    thisii =[];
    
    for dd = 1:length(total_delay_tpts)
       
        for ii= 1:length(ROIs)
            this_ii_idx = all_conds(:,1)==2 & all_subj==ss & all_ROIs==ii;
            [thisii(ii,:,dd), ~] = max(all_fidelity_nodist(this_ii_idx,total_delay_tpts{dd}),[],2);

            for jj = 1:length(ROIs)
                this_jj_idx = all_conds(:,1)==2 & all_subj==ss & all_ROIs==jj;
                [thisjj(jj,:,dd),~] = max(all_fidelity_nodist(this_jj_idx,total_delay_tpts{dd}),[],2);
                r = corrcoef(thisjj(jj,:,dd),thisii(ii,:,dd));
                thiscorr(ii,jj,dd,ss) = r(1,2);
                
            end
        end
     
        subplot(1,length(total_delay_tpts),dd); hold on;
        imagesc(thiscorr(:,:,dd,ss),[-1 1])
        colormap(cmap)
        axis ij tight
        axis square
        set(gca,'YTick',1:length(ROIs),'Yticklabel',{'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'});
        set(gca,'XTick',1:length(ROIs),'Xticklabel',{'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'});
        title(sprintf(' %s - dist',total_epoch{dd}))
        colorbar


    end
        sgtitle(sprintf('Distractor-removed target fidelity, max, subj %s',subj{ss}))  
        if ss == length(subj)
        figure
        for as = 1:length(total_delay_tpts)
        subplot(1,length(total_delay_tpts),as); hold on;
        imagesc(mean(thiscorr( :,:,as,:),4),[-1 1])
        colormap(cmap)
        axis ij tight
        axis square
        set(gca,'YTick',1:length(ROIs),'Yticklabel',{'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'});
        set(gca,'XTick',1:length(ROIs),'Xticklabel',{'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'});
        title(sprintf(' %s - dist',total_epoch{as}))
        colorbar


        end
else 
end
end
end


