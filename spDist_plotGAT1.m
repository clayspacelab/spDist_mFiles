% spDist_plotGAT1.m
%
% for plotting cross-validated WM GAT matrices (via IEM) on distractor trials,
% computed using spDist_channelRespAmp_GATdist.m

% 9/21/2020 UPDATED: COLORMAP, ROIs, GCF POS - geh
root = spDist_loadRoot;
root = '/share/data/spDist/';

subj = {'AY','CC','EK','KD','MR','XL','SF'};
%sess= {{'spDistLong1','spDistLong2'},{'spDistLong1','spDistLong2'},{'spDistLong1','spDistLong2'}};
sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
%ROIs = {'V1','V2','V3','V3AB','hV4','VO1','VO2','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS','iPCS'};
ROIs = {'V1','V2','V3','V3AB','hV4','LO1','IPS0','IPS1','IPS2','IPS3','sPCS'};

func_suffix = 'surf';

nchan = 8;
which_vox = 0.1;

if which_vox < 1
    vox_str = sprintf('_VE%03.f',100*which_vox);
else
    vox_str = sprintf('_%ivox',which_vox);
end

GAT_str = {'Target','Distractor'};

myTR = 0.75;

%% load data
startidx = 1;
for ss = 1:length(subj)
    
    
    
    for vv = 1:length(ROIs)
        
            
            % just one file to load
            fn = sprintf('%sspDist_reconstructions/%s_%s_%s_%s_%ichan%s_GATdist.mat',root,subj{ss},horzcat(sess{ss}{:}),ROIs{vv},func_suffix,nchan,vox_str);
            
            fprintf('loading %s...\n',fn);
            data = load(fn);
            
            
            if vv == 1 && ss == 1
                % initialize variables...
                
                
                nblankt = length(ROIs)*size(data.recons,1);
                
                %all_recons = nan(nblankt,size(data.recons,2),size(data.recons,3));
                all_conds = nan(nblankt,size(data.c_all,2));
                
                % each stimulus type (target/distractor)
                all_fidelity{1} = nan(size(data.recons,1),size(data.recons,2),nblankt); % timecoruse of fidelity
                all_fidelity{2} = nan(size(data.recons,1),size(data.recons,2),nblankt); % timecoruse of fidelity
                
                all_subj = nan(nblankt,1);
                all_ROIs = nan(nblankt,1);
                all_sess = nan(nblankt,1);
                
                
                angs = data.angs;
                tpts = data.delay_tpts;
                
                % ugh have to do this in a multi-D array...
                %all_r2 = nan(length(ROIs),length(tpts),length(subj));
                
            end
            
            
            
            thisidx_all = startidx:(startidx+size(data.c_all,1)-1);
            
            this_fidelity = cellfun(@(x) squeeze(mean(cosd(angs) .* x,2)),data.recons,'UniformOutput',false);
            
            for tt1 = 1:size(this_fidelity,1)
                for tt2 = 1:size(this_fidelity,2)
                    all_fidelity{1}(tt1,tt2,thisidx_all) = this_fidelity{tt1,tt2,1};
                    all_fidelity{2}(tt1,tt2,thisidx_all) = this_fidelity{tt1,tt2,2};
                end
            end
            
            %all_recons(thisidx_map,:,:) = data.recons;
            %all_fidelity(thisidx_map,:) = squeeze(mean(cosd(angs) .* data.recons,2));
            
            %all_r2(vv,:,ss) = squeeze(mean(mean(data.r2_all,1),2)); % average over folds (dim1) and vox (dim 2)
            
            all_conds(thisidx_all,:) = data.c_all;
            
            
            all_subj(thisidx_all) = ss;
            
            
            all_ROIs(thisidx_all) = vv;
            
            all_sess(thisidx_all) = data.sess_all;
            
            
            startidx = thisidx_all(end)+1;
            
            clear data;
            
          
        
    end
    
end


%% plot (average over subj)
% 1x n_rois
figure;
for aa = 1:length(all_fidelity)
    for vv = 1:length(ROIs)
        
        subplot(length(all_fidelity),length(ROIs),vv+(aa-1)*length(ROIs));hold on;
        
        
        thisd = nan(size(all_fidelity{aa},1), size(all_fidelity{aa},2), length(subj));
        for ss = 1:length(subj)
            thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==2;
            thisd(:,:,ss) = squeeze(mean(all_fidelity{aa}(:,:,thisidx),3));
        end
        imagesc(tpts*myTR,tpts*myTR,mean(thisd,3));
        colormap magma
        cbh = colorbar('h','location','eastoutside')
        title(ROIs{vv});
        axis tight square
        set(gca,'XTick',0:4:24,'YTick',0:4:24,'TickDir','out');
        if vv == 1
            xlabel('Test time (s)');
            ylabel('Train time (s)');
        else
            set(gca,'XTickLabel',[],'YTickLabel',[]);
        end
        
    end
end
match_clim(get(gcf,'Children'));
set(cbh,'XTick',[-.2:.1:.3])
%set(gcf,'Position',[32         778        1910         122]);
set(gcf,'position', [ 23         245        2386         453])




%% plot (individual subj as rows)
% NOTE: for now, only matched clim within figure...
for aa = 1:length(all_fidelity)
    figure;
    for vv = 1:length(ROIs)
        
        for ss = 1:length(subj)
            subplot(length(subj),length(ROIs),vv+(ss-1)*length(ROIs));hold on;
            
            
            thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==2;
            thisd = squeeze(mean(all_fidelity{aa}(:,:,thisidx),3));
            
            imagesc(tpts*myTR,tpts*myTR,thisd);
            
            
            if ss == 1
                title(ROIs{vv});
            end
            axis tight square
            if vv == 1
                xlabel('Test time (s)');
                ylabel(sprintf('%s - Train time (s)',subj{ss}));
            end
        end
        
    end
    set(gcf,'Name',GAT_str{aa},'NumberTitle','off');
    match_clim(get(gcf,'Children'));
end


