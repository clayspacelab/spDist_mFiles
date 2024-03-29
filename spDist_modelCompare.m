function spDist_modelCompare(subj,sess,ROIs)
% NOTE: LORO data must already be generated to run this script. see
% spDist_channelRespAmp_GAT.m to generate (can take ~24 hours)
root = spDist_loadRoot;

task_dir = 'spDist';

if nargin < 1 || isempty(subj)
    subj = {'AY','CC','EK','KD','MR','SF','XL'};
end

if nargin < 2 || isempty(sess)
    
    sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
    
    
end

if nargin < 3 || isempty(ROIs)
    ROIs = {'V1','V2','V3','V3AB','hV4','LO1','IPS0','IPS1','IPS2','IPS3','sPCS'};
    %ROIs = {'V1V2V3','V3AB','hV4','LO1','IPS0IPS1','IPS2IPS3','sPCS'};
end


save_stats = 0; 

func_suffix = 'surf';

nchan = 8;

cat_mode = 1; % if 1, look for catSess1Ses...SessN_ files, otherwise, just look for each session in turn

which_vox = 0.1; % if > 1 , look for WHICH_VOXvox string; otherwise, look for VE<100*WHICH_VOX>

smooth_by = 1; % if this is 1, use regular files, otherwise, load smooth_by files

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

rng(spDist_randSeed);

tmpcolors = spDist_condColors;
dist_color = tmpcolors(2,:); clear tmpcolors;
LORO_color = [255,165,0]./255; % orange

model_colors = [dist_color;LORO_color];

% load both original (1: _thruTime) and GAT files trained/tested per epoch (2)
n_files= [1 2]; % how many data files do we care about? % TCS: be careful! you reset this later and use the opposite number scheme

%% load data


for yy = 1:length(n_files)
    startidx = 1;
    for ss = 1:length(subj)
        
        for vv = 1:length(ROIs)
            
            % 'original' files, model estimated w/ held-out data
            if yy ==1
                fn = sprintf('%s/%s_reconstructions/%s_%s_%s_%s_%ichan%s%s%s_recon_thruTime1.mat',root,task_dir,subj{ss},horzcat(sess{ss}{:}),ROIs{vv},func_suffix,nchan,vox_str,smooth_str,trn_str);
                
                fprintf('loading %s...\n',fn);
                data = load(fn);
                
                if vv == 1 && ss == 1 && yy ==1
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
             
            % LORO/GAT model    
            else
                
                fn = sprintf('%sspDist_reconstructions/%s_%s_%s_%s_%ichan%s_GATdist_epochTPTS_gh.mat',root,subj{ss},horzcat(sess{ss}{:}),ROIs{vv},func_suffix,nchan,vox_str);
                
                fprintf('loading %s...\n',fn);
                data = load(fn);
                
                if vv == 1 && ss == 1 && yy ==2
                    
                    % initialize variables...
                    
                    
                    nblankt = length(ROIs)*size(data.recons{1},1);
                    all_recons_gat = cell(size(data.recons));
                    all_fidelity_gat = cell(size(data.recons));
                    

                    for pp =1:size(data.recons,3)
                        for aa =1:size(data.recons,1)
                            for ee = 1:size(data.recons,2)
                                all_recons_gat{aa,ee,pp} = nan(nblankt,size(data.recons{aa},2));
                                all_fidelity_gat{aa,ee,pp} = nan(nblankt,1);
                            end
                        end
                    end
                    
                    all_conds_gat = nan(nblankt,size(data.c_all,2));
                    all_angs_gat = nan(nblankt,size(data.a_all,2));
                    
                    all_subj_gat = nan(nblankt,1);
                    all_ROIs_gat = nan(nblankt,1);
                    all_sess_gat = nan(nblankt,1);
                    all_fn_gat = nan(nblankt,1);
                    
                    angs_gat = data.angs;
                    tpts_gat = data.delay_tpts;
                end
                
                
                % TCS: what are each of these dimensions??
                % 1: training tpt idx
                % 2: testing tpt idx
                % 3: align_to (1 = WM, 2 = distractor)
                
                
                thisidx = startidx:(startidx+size(data.c_all,1)-1);
                for pp =1:size(data.recons,3)
                    for aa =1:size(data.recons,1)
                        for ee = 1:size(data.recons,2)
                            all_recons_gat{aa,ee,pp}(thisidx,:) = data.recons{aa,ee,pp};
                            all_fidelity_gat{aa,ee,pp}(thisidx,:) = squeeze(mean(cosd(angs) .* data.recons{aa,ee,pp},2));
                        end
                    end
                end
                
                
                all_conds_gat(thisidx,:) = data.c_all;
                all_angs_gat(thisidx,:) = data.a_all;
                
                
                all_subj_gat(thisidx) = ss;
                
                
                all_ROIs_gat(thisidx) = vv;
                
                all_sess_gat(thisidx) = data.sess_all;
                
                % TCS: not sure what this is here for? I don't think we use
                % this variable?
                all_fn_gat(thisidx) = yy;
                
                startidx = thisidx(end)+1;
                
                clear data;
                
                
                
            end
                        
            
        end
        
        
    end
end
%% Figure 6C
%plot only like trn/tst combinations

trn_epoch =[1 2 3];
tst_epoch =[1 2 3];
gat_align ={'trn/tst:target/target'};

cond_colors = [ 0 0 1; 0 0 1; 0 0 1];


for pg=1:length(gat_align) % can be length 1 - target aligned recon or length 2, targ aling and dist align. for now, need only target
    figure ('Name','trnlines;Figure6C')
    
    for vv = 1:length(ROIs)
        
        h=[];
        for aa =1:length(trn_epoch)


            % TCS: easier way to do this...
            ee = aa;
            
            thisd = nan(length(subj),size(all_recons_gat{1},2));
            for ss = 1:length(subj)
                
                subplot(size(all_recons_gat,1),length(ROIs),vv+(ee-1)*length(ROIs));hold on;
                
                thisidx = all_subj_gat==ss & all_ROIs_gat==vv & all_conds_gat(:,1)==2;
                thisd(ss,:) = mean(all_recons_gat{aa,ee,pg}(thisidx,:)); %aa = TRN idx, rows; ee = TST idx, col; blue = 1, red =2 , yellow =3
                
            end
            
            % TCS: important to use "angs_gat" here - that's the actual
            % angles used, which aren't precisely in line w/ the linspace
            % command
            % fixed std to use dim 1, not mode 1 (mode 0 is default)
            my_sem = std(thisd,[],1)/sqrt((length(subj)));
            h(aa) = plot(angs_gat,mean(thisd,1),'-','LineWidth',1,'color',LORO_color);
            
            % TCS: updated the colors below to be consistent
            hold on;
            
            plot(angs_gat,mean(thisd,1)+1.*my_sem,'-','LineWidth',.5,'color',LORO_color)
            plot(angs_gat,mean(thisd,1)-1.*my_sem,'-','LineWidth',.5,'color',LORO_color)
            
            
            btwn_fill = [mean(thisd,1)+1.*my_sem fliplr((mean(thisd,1)-1.*my_sem))];
            fill([angs_gat fliplr(angs_gat)],btwn_fill,LORO_color,'linestyle','none','facealpha',0.3);
            
            
            clear thisidx
            clear thisd
            
            
            if ee== 1 && aa==1
                title(ROIs{vv});
            else
            end
            
            if ee==1 && aa==1 && vv ==1
                ylabel('PRE')
                set(gca,'XTick',-180:90:180,'Xticklabel',{'','',''},'TickDir','out');
                set(gca,'YTick',0:0.5:1.5,'Yticklabel',{'','','',''},'TickDir','out')
                
            elseif ee==2 && aa==2 && vv ==1
                ylabel('DIST')
                set(gca,'YTick',0:0.5:1.5,'Yticklabel',{'','','',''},'TickDir','out')
                set(gca,'XTick',-180:90:180,'Xticklabel',{'','',''},'TickDir','out');
            elseif ee==3 && aa==3 && vv ==1
                ylabel('POST')
                set(gca,'XTick',-180:90:180,'Xticklabel',{'-180','-90','0','90','180'},'TickDir','out');
                set(gca,'YTick',0:0.5:1.5,'Yticklabel',{'0','0.5','1.0','1.5'},'TickDir','out')
                ylabel('TRN/TST Matched TPTS, Recon')
                xlabel('Polar angle (\circ)');
            else
                set(gca,'XTick',-180:90:180,'Xticklabel',{'','',''},'TickDir','out');
                set(gca,'YTick',0:0.5:1.5,'Yticklabel',{'','','',''},'TickDir','out')
                
            end
            
        end
        
    end
    
end




set(gcf,'Position',[-132         503        2651         495])
match_ylim(get(gcf,'Children'));
match_xlim(get(gcf,'Children'));
%legend(h, {'PRE','DIST','POST'}) % because we're only showing matched data
%this isn't necessary





%% Figure 6C - plot no longer in use, but keep for stats 06292021
%plot matched trn/tst GAT fidelty data, on same plot, plot independently trained fidelity data
% plot fidelity plotting TRAIN as dv
% which tpts are we plotting throughout?

delay_tpt_range = [3.75 5.25; 8.25 9.75; 10.5 12]; %change middle from 809.5 to 7.5-9 on oct 27 geh
myTR = 0.75;

delay_tpts = cell(size(delay_tpt_range,1),1);
for dd = 1:size(delay_tpt_range,1)
    delay_tpts{dd} = (tpts*myTR) >= delay_tpt_range(dd,1) & (tpts*myTR) < delay_tpt_range(dd,2);
end



plot_indiv = 1; % if 1, plot individual subj lines in line plot, otherwise, don't.
hoffset = 0.15; % +/- this much


% TCS: we need to collect data for an ANOVA, so an "X" with each datapoint,
% and a "labels" variable with [ROI model_type epoch subj]

data_ANOVA   = nan(length(n_files)*length(subj)*size(all_fidelity_gat,1)*length(ROIs),1);
labels_ANOVA = nan(length(n_files)*length(subj)*size(all_fidelity_gat,1)*length(ROIs),4); % [ROI model_type epoch subj]

idx = 1;


gat_align ={'Target aligned'};
for n_files =1:2 % looping over separatee data types
    
    if n_files ==1
        for pg=1:length(gat_align)
            modelcomp =  figure('name','modelcomp');
            
            h=[];
            
            for vv = 1:length(ROIs)
                thisd = []; %recon data
                for ee =1:length(tst_epoch)
                    
                    
                    for aa =1:length(trn_epoch)
                        
                        for ss = 1:length(subj)
                            
                            subplot(1,length(ROIs),vv);hold on;
                            thisidx = all_subj_gat==ss & all_ROIs_gat==vv & all_conds_gat(:,1)==2;
                            thisd(aa,ee,ss) = mean(all_fidelity_gat{aa,ee,pg}(thisidx,:)); %aa = TRN idx, rows; ee = tst idx, blue = 1, red =2 , yellow =3
                            
                            % TCS hack...
                            if aa == ee
                                data_ANOVA(idx) = thisd(aa,ee,ss);
                                labels_ANOVA(idx,:) = [vv 2 ee ss];
                                idx = idx+1;
                            end
                            
                            clear thisidx
                        end
                        
                    end
                end
                hold on;
                h(n_files) = plot([1 2 3], [mean(thisd(1,1,:),3) mean(thisd(2,2,:),3) mean(thisd(3,3,:),3)],'b-','linewidth',.5); % CHANGING THIS FROM -- black to - black, may be confusing w previous figs!!
                
                for ii=1:3
                    my_sem = std(thisd(ii,ii,:),[],3)/sqrt((length(subj)));
                    
                    hold on;
                    plot([ii ii],[mean(thisd(ii,ii,:),3)+1.*my_sem mean(thisd(ii,ii,:),3)-1.*my_sem],'b-','linewidth',.5)
                    
                    clear my_sem
                end
                
                
                                if plot_indiv == 1
                                    plot([1 2 3]+hoffset,[squeeze(thisd(1,1,:)) squeeze(thisd(2,2,:)) squeeze(thisd(3,3,:))],'-','LineWidth',0.25,'Color',[0.4 0.4 0.4]);
                                    hold on;
                                    plot([1 2 3]+hoffset,[squeeze(thisd(1,1,:)) squeeze(thisd(2,2,:)) squeeze(thisd(3,3,:))],'.','markersize',10,'Color','b');
                
                                end
                
                match_ylim(get(gcf,'Children'));
                
                if vv ==1
                    title(ROIs{vv})
                    xlim([0.5 3.4])
                    ylabel('WM Target Fidelity')
                    set(gca,'Xtick',[0 1 2 3 4],'Xticklabel',{'','PRE','DIST','POST',''},'TickDir','out');
                    set(gca,'YTick',0:0.2:.6,'Yticklabel',{'0','.2','.4','.6'},'TickDir','out')
                else
                    set(gca,'Xtick',[0 1 2 3 4],'Xticklabel',{'','','','',''},'XTickLabelRotation',45,'TickDir','out');
                    set(gca,'YTick',0:0.2:.6,'Yticklabel',{'','','',''},'TickDir','out')
                end
            end
        end
        
        
    else
        
        for pg=1:length(gat_align)
            
            
            hh=[];
            
            for vv = 1:length(ROIs)
                thisdata = []; %recon data
                
                for dd =1:length(delay_tpts)
                    
                    for ss = 1:length(subj)
                        
                        subplot(1,length(ROIs),vv);hold on;
                        
                        thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==2;
                        thisdata(dd,ss) = mean(mean( all_fidelity(thisidx,delay_tpts{dd},1) ,2),1); % ensure this is identical to the tpts used for GAT
                        
                        data_ANOVA(idx) = thisdata(dd,ss);
                        labels_ANOVA(idx,:) = [vv 1 dd ss];
                        idx = idx+1;
                        
                        clear thisidx
                    end
                    
                end
                hold on;
                h(n_files) = plot([1 2 3], [mean(thisdata(1,:),2) mean(thisdata(2,:),2) mean(thisdata(3,:),2)],'-.','color', [0 0 1],'linewidth',.5); %collect n_files plot handle for legend use
                
                for ii=1:3
                    my_sem = std(thisdata(ii,:),[],2)/sqrt((length(subj)));
                    
                    hold on;
                    plot([ii ii],[mean(thisdata(ii,:),2)+1.*my_sem mean(thisdata(ii,:),2)-1.*my_sem],'-','color', [0 0 1],'linewidth',.5)
                    clear my_sem
                end
                
                
                xlim([0.5 3.5])
                %ylim([-0.05 .6])
                if vv ==1
                    
                    ylabel('WM target Fidelity')
                    set(gca,'Xtick',[0 1 2 3 4],'Xticklabel',{'','PRE','DIST','POST',''},'XTickLabelRotation',45,'TickDir','out');
                    
                else
                    
                    set(gca,'Xtick',[0 1 2 3 4],'Xticklabel',{'','','','',''},'TickDir','out');
                end
                
                
                clear thisdata
                title(ROIs{vv});
            end
            
            
            
        end
        
    end
end

set(gcf,'Renderer','painters')
set(gcf,'Position',[-132         503        2651         495])
legend(h, 'LORO  model', 'Ind Model')
% sgtitle('Model Comparison')%% do this, over average delay epochs, for each condition


%% try a new version of this now that the data is extracted (data_ANOVA)
% cols: ROI
% subplots: lines are models, x is epoch
%
% (under assumption we're ALWAYS looking at matched training/testing epoch)

model_cond = [1 2]; 


figure('name','Figure6D');
for vv = 1:length(ROIs)
    subplot(1,length(ROIs),vv); hold on;

    % data for each ROI: n_epoch x n_model_cond x n_subj
    thisd = nan(length(delay_tpts),length(model_cond),length(subj));
    
    for cc = 1:length(model_cond) % [1 2]
        
        for ee = 1:length(delay_tpts)
            
            for ss = 1:length(subj)
                % should just be one trial...
                thisd (ee,cc,ss) = data_ANOVA(labels_ANOVA(:,1)==vv & labels_ANOVA(:,2)==cc & labels_ANOVA(:,3)==ee & labels_ANOVA(:,4)==ss);
            end
            
            
        end
    
        thism = squeeze(mean(thisd(:,cc,:),3));
        thise = squeeze(std(thisd(:,cc,:),[],3))/sqrt(length(subj));
        plot(1:length(delay_tpts),thism,'o-','Color',model_colors(cc,:),'LineWidth',1.5,'MarkerSize',5,'MarkerFaceColor',model_colors(cc,:));
        %plot(1:length(delay_tpts),thism+[-1;1].*thise,'-','LineColor',model_colors(cc,:),'LineWidth',1.5);
      % this one  plot([1;1].*(1:length(delay_tpts)),(thism+[-1 1].*thise).','-','Color',model_colors(cc,:),'LineWidth',1.5)
    end
    
    for ee = 1:length(delay_tpts)
        plot(ee+[-1 1]*hoffset,squeeze(thisd(ee,:,:)).','-','Color',[0.4 0.4 0.4],'LineWidth',0.5);
        plot(ee-hoffset       ,squeeze(thisd(ee,1,:)),'.',  'Color',dist_color,'MarkerSize',6); % 255,165,0
        plot(ee+hoffset       ,squeeze(thisd(ee,2,:)),'.',  'Color',LORO_color,'MarkerSize',6);
    end
    
    title(ROIs{vv},'FontWeight','normal','FontAngle','italic'); % TODO: italics
    
    set(gca,'TickDir','out','XTick',1:length(delay_tpts),'YTick',0:0.3:0.9,'XLim',[0.5 length(delay_tpts)+0.5]);
    
    if vv == 1
        set(gca,'XTickLabel',{'PRE','DIST','POST'});
        ylabel('WM target fidelity');
    else
        set(gca,'XTickLabel',[],'YTickLabel',[]);
    end
    
    hold off;
end

match_ylim(get(gcf,'Children'));


%% TCS version of stats

niter = 1000;

% first - 3-way ANOVA against shuffled null


realF_ANOVA33 = RMAOV33_gh([data_ANOVA labels_ANOVA]);
permF_ANOVA33 = nan(niter,length(realF_ANOVA33));

fprintf('Stats - 3-way RMAOV against shuffled null\n');
for ii = 1:niter
    
    this_shuf_data = nan(size(data_ANOVA));
    % shuffle datapoints within each subj
    for ss = 1:length(subj)
        tmpdata = data_ANOVA(labels_ANOVA(:,end)==ss);
        shufidx = randperm(length(tmpdata));
        this_shuf_data(labels_ANOVA(:,end)==ss) = tmpdata(shufidx);
        clear tmpdata shufidx;
    end
    
    
    permF_ANOVA33(ii,:) = RMAOV33_gh([this_shuf_data labels_ANOVA]);
    clear this_shuf_data;
end

p_ANOVA33 = mean(permF_ANOVA33>=realF_ANOVA33,1);
IV_labels = {'ROI','model','epoch'};

% IV1, IV2, IV3, IV1xIV2, IV1xIV3, IV2xIV3, IV1xIV2xIV3
fprintf('ME ROI:      p = %0.03f\n',p_ANOVA33(1));
fprintf('ME model:    p = %0.03f\n',p_ANOVA33(2));
fprintf('ME epoch:    p = %0.03f\n',p_ANOVA33(3));
fprintf('ROI x model: p = %0.03f\n',p_ANOVA33(4));
fprintf('ROI x epoch: p = %0.03f\n',p_ANOVA33(5));
fprintf('mod x epoch: p = %0.03f\n',p_ANOVA33(6));
fprintf('3-way inter: p = %0.03f\n',p_ANOVA33(7));

% loop over ROIs, 2-way ANOVA against shuffled null
realF_ANOVA2 = nan(length(ROIs),3);
permF_ANOVA2 = nan(length(ROIs),3,niter);

p_ANOVA2 = nan(length(ROIs),3);

for vv = 1:length(ROIs)
    
    fprintf('Stats - 2-way RMAOV against shuffled null: %s\n',ROIs{vv});
    
    ROI_data = data_ANOVA(labels_ANOVA(:,1)==vv);
    ROI_labels = labels_ANOVA(labels_ANOVA(:,1)==vv,2:end);
    
    % real 2-way ANOVA
    realF_ANOVA2(vv,:) = RMAOV2_gh([ROI_data ROI_labels]);
    %tic;
    for ii = 1:niter
        
        ROI_data_shuf = nan(size(ROI_data));
        for ss = 1:length(subj)
            
            tmpdata = ROI_data(ROI_labels(:,end)==ss);
            shufidx = randperm(length(tmpdata));
            ROI_data_shuf(ROI_labels(:,end)==ss) = tmpdata(shufidx);
            clear tmpdata shufidx;
            
        end
        
        permF_ANOVA2(vv,:,ii) = RMAOV2_gh([ROI_data_shuf ROI_labels]);
        
        clear ROI_data_shuf;
        
    end
    
    p_ANOVA2(vv,:) = mean(permF_ANOVA2(vv,:,:)>=realF_ANOVA2(vv,:),3);
    
    fprintf('ME model: p = %0.03f\n',p_ANOVA2(vv,1));
    fprintf('ME epoch: p = %0.03f\n',p_ANOVA2(vv,2));
    fprintf('interact: p = %0.03f\n',p_ANOVA2(vv,3));
    
    %toc;
end

% FDR thresh
fdr_thresh_ANOVA2 = nan(1,size(p_ANOVA2,2));
for ii = 1:size(p_ANOVA2,2)
    fdr_thresh_ANOVA2(ii) = fdr(p_ANOVA2(:,ii),0.05);
end

fprintf('\nFDR thresh - model: p = %0.03f, epoch: p = %0.03f, interaction: p = %0.03f\n',fdr_thresh_ANOVA2(1),fdr_thresh_ANOVA2(2),fdr_thresh_ANOVA2(3));



if save_stats == 1
    fn2s = sprintf('%s/spDist_stats/n%i_GATEpochStats_shuf_%iIter_%s.mat',root,length(subj),niter,datestr(now,30));
    fprintf('Saving stats to %s\n',fn2s);
    save(fn2s,'realF_ANOVA33','permF_ANOVA33','p_ANOVA33','IV_labels','realF_ANOVA2','permF_ANOVA2','p_ANOVA2','fdr_thresh_ANOVA2','subj','ROIs','delay_tpt_range');
end

% geh adding stats markers 01/18/21

stats_markers ={'M','E','X'};
stats_pos = [-1.30 -.03; -1.0 -.04; -.75 -.03]; % in plotted units, x,y, relative to top right

signif_color = [0 0 0];
trend_color  = [0.5 0.5 0.5];

for vv = 1:length(ROIs)
    for ii =1:3 % effect #; IV1, IV2, IV1 x IV2, MODEL, EPOCH, MODEL x EPOCH
        figure(3); hold on;
        subplot(1,length(ROIs),vv)
        
        tmpxlim = get(gca,'XLim');
        tmpylim = get(gca,'YLim');
        
        topright = [tmpxlim(2) tmpylim(2)];
        
        if p_ANOVA2(vv,ii) <= fdr_thresh_ANOVA2(ii)
            text(topright(1)+stats_pos(ii,1),topright(2)+stats_pos(ii,2),stats_markers{ii},'FontSize',18,'Color',signif_color);
        elseif p_ANOVA2(vv,ii) <= 0.05
            text(topright(1)+stats_pos(ii,1),topright(2)+stats_pos(ii,2),stats_markers{ii},'FontSize',18,'Color',trend_color);
        end
        
        
    end
    
end




if 0


%% 2-way & 1-way true ANOVA w/_thruTime1 data
% get fidelity , over average delay epochs, for each condition, permutation test
%%%% this takes ~24 minutes to run!!!!!

thisfide_store = [];
the_y_store= [];

for dd = 1:length(delay_tpts)
    for vv = 1:length(ROIs)
        
        for ss = 1:length(subj)
            thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==2;
            thisfide_store(dd,vv,ss) = mean(mean(all_fidelity(thisidx,delay_tpts{dd},1),2));  %targ aligned
            the_y_store = [the_y_store;   thisfide_store(dd,vv,ss)  dd vv  ss];
        end
        
    end
end



% sort y_store into individual columns for transparency
y = the_y_store(:,1);
epoch_var = the_y_store(:,2);
roi_var =the_y_store(:,3);
subj_var = the_y_store(:,4);


P_truth_1= cell(length(ROIs),1);
T_truth_1 = cell(length(ROIs),1);
F_store_truth_1 = nan(length(ROIs),1);

P_truth_2= cell(1,3);
T_truth_2 = cell(1,3);
F_store_truth_2 = nan(1,3);


%perform the true ANOVA

thisy = y;
thisepoch = epoch_var;
thisroi= roi_var;
thissubj=subj_var;

[P_truth_2,T_truth_2,~] = anovan(thisy,{thissubj,thisepoch,thisroi},'model','full','random',1,'varnames',{'subj','epoch','roi'},'display','off');
[F_store_truth_2] = RMAOV2_gh([thisy,thisepoch,thisroi,thissubj],0.05);

% perform 1-way anova on a per ROI basis, with epoch as factor
for vv = 1:length(ROIs)
    
    thisroiidx = roi_var ==vv;
    thisy = y(thisroiidx);
    thisepoch = epoch_var(thisroiidx);
    thissubj=subj_var(thisroiidx);
    [P_truth_1{vv},T_truth_1{vv},~] = anovan(thisy,{thissubj,thisepoch},'random',1,'varnames',{'subj','epoch'},'display','off');
    [F_store_truth_1(vv,:)] = RMAOV1_gh([thisy,thisepoch,thissubj],0.05);
    
    clear thisroiidx thisepoch thiscond thissubj thisy
    
    
end





%%  2-way permuation ANOVA w/_thruTime1 data

thisy = y;
thisepoch = epoch_var;
thisroi= roi_var;
thissubj=subj_var;
iter = 1000;
fprintf(sprintf('computing %i 1 & 2-way ANOVAs on _thruTime1',iter))

T_iter_1 =cell(iter,1);
P_iter_1 =cell(iter,1);
F_store_iter_1 =nan(iter,1);

T_iter_2 =cell(iter,1);
P_iter_2 =cell(iter,1);
F_store_iter_2 =nan(iter,3);

tic

for xx=1:iter
    y_shuf=y;
    for sbj=1:length(subj)
        thisidx =  thissubj==sbj;
        tmp_y = y(thisidx);
        shuff_idx = randperm(size(tmp_y,1))';
        thisidx_val =find(thisidx);
        y_shuf(thisidx_val,:) = tmp_y(shuff_idx,:);
    end
    
    thisepoch = epoch_var;
    thissubj=subj_var;
    thisroi =roi_var;
    
    [P_iter_2{xx},T_iter_2{xx},~] = anovan(y_shuf,{thissubj,thisepoch,thisroi},'model','full','random',1,'varnames',{'subj','epoch','roi'},'display','off');
    [F_store_iter_2(xx,:)] = RMAOV2_gh([y_shuf,thisepoch,thisroi,thissubj],0.05);
    
    clear thisidx thisy tmpy shuff_idx
    
    
end
toc

%plot F distributions
%which F vals do we want? idx is as follows {iter,1}{var_str,6} %1 is
%fixed because we store iters vertically(1 column). 6 is fixed bc this is the Fval col
var_str ={'Source';'subj';'epoch';'roi';'subj*epoch';...
    'subj*roi';'epoch*roi';...
    'subj*epoch*roi';...
    'Error';'Total'}


col_idx = 6; % F val col

%what factor do we care about? concurs w var_str(vs)
vs = [3 4 7]; %3 epoch 4 roi 7epoch*cond
extract_store=[];

figure
for which_test =1:length(vs)
    extract_vals=[];
    for ii=1:iter
        extract_vals = [extract_vals; T_iter_2{ii}{vs(which_test),col_idx}];
    end
    hold on;
    
    subplot(1,length(vs),which_test)
    histogram(extract_vals)
    line([T_truth_2{vs(which_test),col_idx} T_truth_2{vs(which_test),col_idx}], [min(ylim) max(ylim)], 'color','r','linewidth',2,'linestyle','--')
    title(var_str{vs(which_test)},'Interpreter','none');
    exact_p = sum(extract_vals >= T_truth_2{vs(which_test),col_idx})/iter;
    text(max(xlim)-(.5*max(xlim)),max(ylim)-(.25*max(ylim)),sprintf('%i',exact_p),'FontSize',9)
    exact_store(:,which_test) = [extract_store; exact_p];
    clear exact_p
    if which_test ==1
        xlabel('T-stat')
        ylabel('Frequency of T-stat')
    else
    end
    
end

ta = table(exact_store(:,1),exact_store(:,2),exact_store(:,3));
ta.Properties.VariableNames={'Epoch','ROI','EpochROI'}
%% follow-up 1-way test within each ROI

thisfide_store = [];
the_y_store= [];

T_iter_1 =cell(length(ROIs),iter);
P_iter_1 =cell(length(ROIs),iter);
F_store_iter_1 =nan(length(ROIs),iter);


for dd = 1:length(delay_tpts)
    for vv = 1:length(ROIs)
        
        for ss = 1:length(subj)
            thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==2;
            thisfide_store(dd,vv,ss) = mean(mean(all_fidelity(thisidx,delay_tpts{dd},1),2));  %targ aligned
            the_y_store = [the_y_store;   thisfide_store(dd,vv,ss)  dd vv  ss];
        end
        
        
        
        
    end
end

y = the_y_store(:,1);
epoch_var = the_y_store(:,2);
roi_var =the_y_store(:,3);
subj_var = the_y_store(:,4);

for xx=1:iter
    for vv = 1:length(ROIs)
        y_shuf =y;
        thisroiidx = roi_var ==vv;
        for sbj = 1:length(subj)
            thisidx = roi_var ==vv & subj_var==sbj;
            tmp_y = y(thisidx);
            shuff_idx = randperm(size(tmp_y,1))';
            thisidx_val = find(thisidx);
            y_shuf(thisidx_val,:) = tmp_y(shuff_idx,:);
        end
        thisepoch = epoch_var(thisroiidx);
        
        thissubj=subj_var(thisroiidx);
        
        [P_iter_1{vv,xx},T_iter_1{vv,xx},~] = anovan(y_shuf(thisroiidx),{thissubj,thisepoch},'random',1,'varnames',{'subj','epoch'},'display','off');
        [F_store_iter_1(vv,xx,:)] = RMAOV1_gh([y_shuf(thisroiidx),thisepoch,thissubj],0.05);
        
        clear thisroiidx thisepoch thiscond thissubj thisy tmpy shuffidx
        
        
    end
end
toc


%plot F distributions
%which F vals do we want? idx is as follows {iter,1}{var_str,6} %1 is
%fixed because we store iters vertically(1 column). 6 is fixed bc this is the Fval col
iv_str ={'Source';'subj';'epoch';'Error';'Total'};


col_idx = 6; % F val col

%what factor do we care about? concurs w var_str(vs)
vs = [3]; %3 epoch
extract_store=[];
% here, i want to see the results of the two anova functions, separately
%%%%%%%%%%%%%%%%%%%%%%%%% do for anovan

which_effect = [3]; %3 epoch
exact_store_tmp=[];
figure('name','1-way perm;anovan')
for vv =1:length(ROIs)
    
    for ww =1:length(which_effect)
        extract_vals=[];
        for ii=1:iter
            extract_vals = [extract_vals; T_iter_1{vv,ii}{which_effect(ww),col_idx}];
        end
        hold on;
        
        subplot(1,length(ROIs),vv)
        histogram(extract_vals)
        line([T_truth_1{vv}{which_effect(ww),col_idx} T_truth_1{vv}{which_effect(ww),col_idx}], [min(ylim) max(ylim)], 'color','r','linewidth',2,'linestyle','--')
        title(iv_str{which_effect(ww)},'Interpreter','none');
        exact_p = sum(extract_vals >= T_truth_1{vv}{which_effect(ww),col_idx})/iter;
        text(max(xlim)-(.5*max(xlim)),max(ylim)-(.25*max(ylim)),sprintf('%i',exact_p),'FontSize',9)
        exact_store_1(vv,ww) = [exact_store_tmp; exact_p];
        clear exact_p
        
        if ww ==1
            xlabel('T-stat')
            ylabel('Frequency of T-stat')
        else
        end
        
    end
end

ta_anovan1 = table(ROIs',exact_store_1(:,1)); % outcome of 1-way shuffle derived p-vals for _thruTime, compare w/ exact_store_rma
ta_anovan1.Properties.VariableNames={'ROI','Epoch'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% do for RMAOV1

which_effect =1; % col of main effect
exact_store_tmp=[];
figure('name','1-way perm;RMAOV1')

for vv =1:length(ROIs)
    
    for ww =1:length(which_effect)
        extract_vals=[];
        for ii=1:iter
            extract_vals = [extract_vals; F_store_iter_1(vv,ii)];
        end
        hold on;
        
        subplot(1,length(ROIs),vv)
        histogram(extract_vals)
        line([F_store_truth_1(vv,which_effect(ww)) F_store_truth_1(vv,which_effect(ww))], [min(ylim) max(ylim)], 'color','r','linewidth',2,'linestyle','--')
        title(iv_str{which_effect(ww)},'Interpreter','none');
        exact_p = sum(extract_vals >= F_store_truth_1(vv,which_effect(ww)))/iter;
        text(max(xlim)-(.5*max(xlim)),max(ylim)-(.25*max(ylim)),sprintf('%i',exact_p),'FontSize',9)
        exact_store_rma(vv,ww) = [exact_store_tmp; exact_p];
        clear exact_p
        
        if ww ==1
            xlabel('T-stat')
            ylabel('Frequency of T-stat')
        else
        end
        
    end
end

ta_rma1 = table(ROIs',exact_store_rma(:,1));
ta_rma1.Properties.VariableNames={'ROI','Epoch'};
%% plot sig markers on the subplots - based on anovan output


[p_fdr_perm, p_masked_perm] = fdr(exact_store_1,0.05);


sig_mrkr ={'o'};
y_mod = [.1];
for vv = 1:length(ROIs)
    for ee =1
        figure(modelcomp); hold on;
        subplot(1,length(ROIs),vv)
        
        
        if exact_store_1(vv,ee) > p_fdr_perm(ee) && exact_store_1(vv,ee) <= 0.05
            text(max(xlim)-(.2*max(xlim)),max(ylim)-(y_mod(ee)*max(ylim)),sprintf('%s',sig_mrkr{ee}),'color',[.5 .5 .5],'fontsize',15) %unfilled grey for IND
        elseif  exact_store_1(vv,ee) <= p_fdr_perm(ee)
            text(max(xlim)-(.2*max(xlim)),max(ylim)-(y_mod(ee)*max(ylim)),sprintf('%s',sig_mrkr{ee}),'color','k','fontsize',15) %filled black for IND
        else
        end
        
        
    end
    
end



fprintf('1 & 2-way ANOVA permutations on _thru data are complete. ')


%% 1& 2-way true ANOVA, roi & roi x epoch, w/_GAT data
% get fidelity , over average delay epochs, for each condition, permutation test


thisfide_store = [];
the_y_store= [];
aa = [];

for vv = 1:length(ROIs)
    for ee = 1:3
        
        for ss = 1:length(subj)
            
            thisidx = all_subj_gat==ss & all_ROIs_gat==vv & all_conds_gat(:,1)==2;
            
            if ee==1
                thisfide_store(vv,ee,ss) = mean(all_fidelity_gat{1,ee,1}(thisidx,:));
            elseif ee==2
                thisfide_store(vv,ee,ss) = mean(all_fidelity_gat{2,ee,1}(thisidx,:));
            elseif ee==3
                thisfide_store(vv,ee,ss) = mean(all_fidelity_gat{3,ee,1}(thisidx,:));
            else
            end
            
            the_y_store = [the_y_store;   thisfide_store(vv,ee,ss)  vv ee ss];
        end
        
        
        
        
    end
end


% sort y_store into individual columns for transparency
y = the_y_store(:,1);
roi_var = the_y_store(:,2);
trntst_var = the_y_store(:,3);
subj_var = the_y_store(:,4);

P_truth_gat_1= cell(length(ROIs),1);
T_truth_gat_1= cell(length(ROIs),1);
F_store_truth_gat_1 = nan(length(ROIs),1);

P_truth_gat_2= cell(1,3);
T_truth_gat_2= cell(1,3);
F_store_truth_gat_2 = nan(1,3);


%perform the "real" ANOVA

thisy = y;
thisepoch = trntst_var;
thisroi= roi_var;
thissubj=subj_var;
[P_truth_gat_2,T_truth_gat_2,~] = anovan(thisy,{thissubj,thisepoch,thisroi},'model','full','random',1,'varnames',{'subj','epoch','roi'},'display','off');
[F_store_truth_gat_2] = RMAOV2_gh([thisy,thisepoch,thisroi,thissubj],0.05);

clear thisroiidx

% 1-way within each ROI
for vv = 1:length(ROIs)
    
    thisroiidx = roi_var ==vv;
    thisy = y(thisroiidx);
    thisepoch = trntst_var(thisroiidx);
    thissubj=subj_var(thisroiidx);
    [P_truth_gat_1{vv},T_truth_gat_1{vv},~] = anovan(thisy,{thissubj,thisepoch},'random',1,'varnames',{'subj','epoch'},'display','off');
    [F_store_truth_gat_1(vv,:)] = RMAOV1_gh([thisy,thisepoch,thissubj],0.05);
    
    clear thisroiidx thisepoch thiscond thissubj thisy
    
    
end




fprintf(sprintf('computing %i 1 & 2-way ANOVAs on _GAT data',iter))

T_iter_gat_2 =cell(iter,1);
P_iter_gat_2 =cell(iter,1);
F_store_iter_gat_2 =nan(iter,3);

T_iter_gat_1 =cell(length(ROIs),iter);
P_iter_gat_1 =cell(length(ROIs),iter);
F_store_iter_gat_1 =nan(length(ROIs),iter);
% two-way anova for each ROI
tic

thisy = y;
thisepoch = trntst_var;
thisroi= roi_var;
thissubj=subj_var;


for xx=1:iter
    y_shuf=y;
    for sbj=1:length(subj)
        thisidx =  thissubj==sbj;
        tmp_y = y(thisidx);
        shuff_idx = randperm(size(tmp_y,1))';
        thisidx_val =find(thisidx);
        y_shuf(thisidx_val,:) = tmp_y(shuff_idx,:);
    end
    
    thisepoch = trntst_var;
    thissubj=subj_var;
    thisroi =roi_var;
    
    [P_iter_gat_2{xx},T_iter_gat_2{xx},~] = anovan(y_shuf,{thissubj,thisepoch,thisroi},'model','full','random',1,'varnames',{'subj','epoch','roi'},'display','off');
    [F_store_iter_gat_2(xx,:)] = RMAOV2_gh([y_shuf,thisepoch,thisroi,thissubj],0.05);
    
    
    clear thisidx thisy tmpy shuff_idx
    
    
end
toc

%plot F distributions
%which F vals do we want? idx is as follows {iter,1}{var_str,6} %1 is
%fixed because we store iters vertically(1 column). 6 is fixed bc this is the Fval col
var_str ={'Source';'subj';'epoch';'roi';'subj*epoch';...
    'subj*roi';'epoch*roi';...
    'subj*epoch*roi';...
    'Error';'Total'}


col_idx = 6; % F val col

%what factor do we care about? concurs w var_str(vs)
vs = [3 4 7]; %3 epoch 4 roi 7epoch*roi
extract_store=[];

figure
for which_test =1:length(vs)
    extract_vals=[];
    for ii=1:iter
        extract_vals = [extract_vals; T_iter_gat_2{ii}{vs(which_test),col_idx}];
    end
    hold on;
    
    subplot(1,length(vs),which_test)
    histogram(extract_vals)
    line([T_truth_gat_2{vs(which_test),col_idx} T_truth_gat_2{vs(which_test),col_idx}], [min(ylim) max(ylim)], 'color','r','linewidth',2,'linestyle','--')
    title(var_str{vs(which_test)},'Interpreter','none');
    exact_p = sum(extract_vals >= T_truth_gat_2{vs(which_test),col_idx})/iter;
    text(max(xlim)-(.5*max(xlim)),max(ylim)-(.25*max(ylim)),sprintf('%i',exact_p),'FontSize',9)
    exact_store(:,which_test) = [extract_store; exact_p];
    clear exact_p
    if which_test ==1
        xlabel('T-stat')
        ylabel('Frequency of T-stat')
    else
    end
    
end

ta_gat_2 = table(exact_store(:,1),exact_store(:,2),exact_store(:,3));
ta_gat_2.Properties.VariableNames={'Epoch','ROI','EpochROI'};

% 1-way perm w _GAT
y = the_y_store(:,1);
roi_var = the_y_store(:,2);
epoch_var = the_y_store(:,3);
subj_var = the_y_store(:,4);

for xx=1:iter
    for vv = 1:length(ROIs)
        y_shuf =y;
        thisroiidx = roi_var ==vv;
        for sbj = 1:length(subj)
            thisidx = roi_var ==vv & subj_var==sbj;
            tmp_y = y(thisidx);
            shuff_idx = randperm(size(tmp_y,1))';
            thisidx_val = find(thisidx);
            y_shuf(thisidx_val,:) = tmp_y(shuff_idx,:);
        end
        thisepoch = epoch_var(thisroiidx);
        
        thissubj=subj_var(thisroiidx);
        
        [P_iter_gat_1{vv,xx},T_iter_gat_1{vv,xx},~] = anovan(y_shuf(thisroiidx),{thissubj,thisepoch},'random',1,'varnames',{'subj','epoch'},'display','off');
        [F_store_iter_gat_1(vv,xx,:)] = RMAOV1_gh([y_shuf(thisroiidx),thisepoch,thissubj],0.05);
        
        clear thisroiidx thisepoch thiscond thissubj thisy tmpy shuffidx
        
        
    end
end
toc


%plot F distributions
%which F vals do we want? idx is as follows {iter,1}{var_str,6} %1 is
%fixed because we store iters vertically(1 column). 6 is fixed bc this is the Fval col
iv_str ={'Source';'subj';'epoch';'Error';'Total'};


col_idx = 6; % F val col

%what factor do we care about? concurs w var_str(vs)
vs = [3]; %3 epoch
extract_store=[];

%%%%%%%%%%%%%%%%%%%%%%%%% do for anovan
which_effect = [3];
exact_store_tmp=[];
figure('name','1-way perm;anovan GAT')
for vv =1:length(ROIs)
    
    for ww =1:length(which_effect)
        extract_vals=[];
        for ii=1:iter
            extract_vals = [extract_vals; T_iter_gat_1{vv,ii}{which_effect(ww),col_idx}];
        end
        hold on;
        
        subplot(1,length(ROIs),vv)
        histogram(extract_vals)
        line([T_truth_gat_1{vv}{which_effect(ww),col_idx} T_truth_gat_1{vv}{which_effect(ww),col_idx}], [min(ylim) max(ylim)], 'color','r','linewidth',2,'linestyle','--')
        title(iv_str{which_effect(ww)},'Interpreter','none');
        exact_p = sum(extract_vals >= T_truth_gat_1{vv}{which_effect(ww),col_idx})/iter;
        text(max(xlim)-(.5*max(xlim)),max(ylim)-(.25*max(ylim)),sprintf('%i',exact_p),'FontSize',9)
        exact_store_gat_1(vv,ww) = [exact_store_tmp; exact_p];
        clear exact_p
        
        if ww ==1
            xlabel('T-stat')
            ylabel('Frequency of T-stat')
        else
        end
        
    end
end

ta_gat_anovan1 = table(ROIs',exact_store_gat_1(:,1));
ta_gat_anovan1.Properties.VariableNames= {'ROI','Epoch'}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% do for RMAOV1
which_effect =1;
exact_store_tmp=[];
figure('name','1-way perm;RMAOV1 GAT')
for vv =1:length(ROIs)
    
    for ww =1:length(which_effect)
        extract_vals=[];
        for ii=1:iter
            extract_vals = [extract_vals; F_store_iter_gat_1(vv,ii)];
        end
        hold on;
        
        subplot(1,length(ROIs),vv)
        histogram(extract_vals)
        line([F_store_truth_gat_1(vv,which_effect(ww)) F_store_truth_gat_1(vv,which_effect(ww))], [min(ylim) max(ylim)], 'color','r','linewidth',2,'linestyle','--')
        title(iv_str{which_effect(ww)},'Interpreter','none');
        exact_p = sum(extract_vals >= F_store_truth_gat_1(vv,which_effect(ww)))/iter;
        text(max(xlim)-(.5*max(xlim)),max(ylim)-(.25*max(ylim)),sprintf('%i',exact_p),'FontSize',9)
        exact_store_gat_rma(vv,ww) = [exact_store_tmp; exact_p];
        clear exact_p
        
        if ww ==1
            xlabel('T-stat')
            ylabel('Frequency of T-stat')
        else
        end
        
    end
end

ta_gat_rma1 = table(ROIs',exact_store_gat_rma(:,1));
ta_gat_rma1.Properties.VariableNames={'ROI','Epoch'};



[p_fdr_perm_gat, p_masked_perm_gat] = fdr(exact_store_gat_1,0.05);


sig_mrkr ={'+'};
y_mod = [.15];
for vv = 1:length(ROIs)
    for ee =1 % only one effect, 'epoch'
        figure(modelcomp); hold on;
        subplot(1,length(ROIs),vv)
        
        
        if exact_store_gat_1(vv,ee) > p_fdr_perm_gat(ee) && exact_store_gat_1(vv,ee) <= 0.05
            text(max(xlim)-(.2*max(xlim)),max(ylim)-(y_mod(ee)*max(ylim)),sprintf('%s',sig_mrkr{ee}),'color',[.5 .5 .5],'fontsize',15) %unfilled grey for GAT
        elseif  exact_store_gat_1(vv,ee) <= p_fdr_perm_gat(ee)
            text(max(xlim)-(.2*max(xlim)),max(ylim)-(y_mod(ee)*max(ylim)),sprintf('%s',sig_mrkr{ee}),'color','k','fontsize',15) %filled black for GAT
        else
        end
        
        
    end
    
end


fprintf('the end')
end

end

