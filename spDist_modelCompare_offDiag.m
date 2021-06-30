function spDist_modelCompare_offDiag(subj,sess,ROIs)
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
   % ROIs = {'V1','V2','V3','V3AB','hV4','LO1','IPS0','IPS1','IPS2','IPS3','sPCS'};
  ROIs = {'V1V2V3','V3AB','hV4','LO1','IPS0IPS1','IPS2IPS3','sPCS'};
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
n_files= [1 2 3]; % how many data files do we care about? % TCS: be careful! you reset this later and use the opposite number scheme
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
            elseif yy == 2
                
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
            elseif  yy ==3    
                fn = sprintf('%sspDist_reconstructions/%s_%s_%s_%s_%ichan%s_GATdist_epochTPTS_gh_shuf1000.mat',root,subj{ss},horzcat(sess{ss}{:}),ROIs{vv},func_suffix,nchan,vox_str);
                
                fprintf('loading %s...\n',fn);
                data = load(fn);
                
                if vv == 1 && ss == 1 && yy ==3
                    
                    % initialize variables...
                    
                    
                    nblankt = length(ROIs)*size(data.recons{1},1);
                    all_recons_gat_shuf = cell(size(data.recons));
                    all_fidelity_gat_shuf = cell(size(data.recons));
                    

                    for pp =1:size(data.recons,3)
                        for aa =1:size(data.recons,1)
                            for ee = 1:size(data.recons,2)
                                all_recons_gat_shuf{aa,ee,pp} = nan(nblankt,size(data.recons{aa},2));
                                all_fidelity_gat_shuf{aa,ee,pp} = nan(nblankt,1);
                            end
                        end
                    end
                    
                    all_conds_gat_shuf = nan(nblankt,size(data.c_all,2));
                    all_angs_gat_shuf = nan(nblankt,size(data.a_all,2));
                    
                    all_subj_gat_shuf = nan(nblankt,1);
                    all_ROIs_gat_shuf = nan(nblankt,1);
                    all_sess_gat_shuf = nan(nblankt,1);
                    all_fn_gat_shuf = nan(nblankt,1);
                    
                    angs_gat_shuf = data.angs;
                    tpts_gat_shuf = data.delay_tpts;
                end
                
                
                % TCS: what are each of these dimensions??
                % 1: training tpt idx
                % 2: testing tpt idx
                % 3: align_to (1 = WM, 2 = distractor)
                
                
                thisidx = startidx:(startidx+size(data.c_all,1)-1);
                for pp =1:size(data.recons,3)
                    for aa =1:size(data.recons,1)
                        for ee = 1:size(data.recons,2)
                            all_recons_gat_shuf{aa,ee,pp}(thisidx,:) = data.recons{aa,ee,pp};
                            all_fidelity_gat_shuf{aa,ee,pp}(thisidx,:) = squeeze(mean(cosd(angs) .* data.recons{aa,ee,pp},2));
                        end
                    end
                end
                
                
                all_conds_gat_shuf(thisidx,:) = data.c_all;
                all_angs_gat_shuf(thisidx,:) = data.a_all;
                
                
                all_subj_gat_shuf(thisidx) = ss;
                
                
                all_ROIs_gat_shuf(thisidx) = vv;
                
                all_sess_gat_shuf(thisidx) = data.sess_all;
                
                % TCS: not sure what this is here for? I don't think we use
                % this variable?
                % wasnt sure if i would need it
                all_fn_gat_shuf(thisidx) = yy;
                
                startidx = thisidx(end)+1;
                
                clear data;
                
                
            end
                        
            
        end
        
        
    end
end
%% Figure 6B - ORIGINAL (MAINTAINING FOR REFERENCE)
%plot only like trn/tst combinations
% 
% trn_epoch =[1 2 3];
% tst_epoch =[1 2 3];
% gat_align ={'trn/tst:target/target'};
% 
% cond_colors = [ 0 0 1; 0 0 1; 0 0 1];
% 
% new_col = [LORO_color+.1; LORO_color+.2;LORO_color+.3;];
% for pg=1:length(gat_align) % can be length 1 - target aligned recon or length 2, targ aling and dist align. for now, need only target
%     figure ('Name','trnlines;Figure6C')
%     
%     for vv = 1:length(ROIs)
%         
%         h=[];
%         for aa =1:length(trn_epoch)
% 
% 
%             % TCS: easier way to do this...
%             ee = aa;
%             
%             thisd = nan(length(subj),size(all_recons_gat{1},2));
%             for ss = 1:length(subj)
%                 
%                 subplot(size(all_recons_gat,1),length(ROIs),vv+(ee-1)*length(ROIs));hold on;
%                 
%                 thisidx = all_subj_gat==ss & all_ROIs_gat==vv & all_conds_gat(:,1)==2;
%                 thisd(ss,:) = mean(all_recons_gat{aa,ee,pg}(thisidx,:)); %aa = TRN idx, rows; ee = TST idx, col; blue = 1, red =2 , yellow =3
%                 
%             end
%             
%             % TCS: important to use "angs_gat" here - that's the actual
%             % angles used, which aren't precisely in line w/ the linspace
%             % command
%             % fixed std to use dim 1, not mode 1 (mode 0 is default)
%             my_sem = std(thisd,[],1)/sqrt((length(subj)));
%             h(aa) = plot(angs_gat,mean(thisd,1),'-','LineWidth',1,'color',LORO_color);
%             
%             % TCS: updated the colors below to be consistent
%             hold on;
%             
%             plot(angs_gat,mean(thisd,1)+1.*my_sem,'-','LineWidth',.5,'color',LORO_color)
%             plot(angs_gat,mean(thisd,1)-1.*my_sem,'-','LineWidth',.5,'color',LORO_color)
%             
%             
%             btwn_fill = [mean(thisd,1)+1.*my_sem fliplr((mean(thisd,1)-1.*my_sem))];
%             fill([angs_gat fliplr(angs_gat)],btwn_fill,LORO_color,'linestyle','none','facealpha',0.3);
%             
%             
%             clear thisidx
%             clear thisd
%             
%             
%             if ee== 1 && aa==1
%                 title(ROIs{vv});
%             else
%             end
%             
%             if ee==1 && aa==1 && vv ==1
%                 ylabel('PRE')
%                 set(gca,'XTick',-180:90:180,'Xticklabel',{'','',''},'TickDir','out');
%                 set(gca,'YTick',0:0.5:1.5,'Yticklabel',{'','','',''},'TickDir','out')
%                 
%             elseif ee==2 && aa==2 && vv ==1
%                 ylabel('DIST')
%                 set(gca,'YTick',0:0.5:1.5,'Yticklabel',{'','','',''},'TickDir','out')
%                 set(gca,'XTick',-180:90:180,'Xticklabel',{'','',''},'TickDir','out');
%             elseif ee==3 && aa==3 && vv ==1
%                 ylabel('POST')
%                 set(gca,'XTick',-180:90:180,'Xticklabel',{'-180','-90','0','90','180'},'TickDir','out');
%                 set(gca,'YTick',0:0.5:1.5,'Yticklabel',{'0','0.5','1.0','1.5'},'TickDir','out')
%                 ylabel('TRN/TST Matched TPTS, Recon')
%                 xlabel('Polar angle (\circ)');
%             else
%                 set(gca,'XTick',-180:90:180,'Xticklabel',{'','',''},'TickDir','out');
%                 set(gca,'YTick',0:0.5:1.5,'Yticklabel',{'','','',''},'TickDir','out')
%                 
%             end
%             
%         end
%         
%     end
%     
% end
% 
% 
% 
% 
% set(gcf,'Position',[-132         503        2651         495])
% match_ylim(get(gcf,'Children'));
% match_xlim(get(gcf,'Children'));
% %legend(h, {'PRE','DIST','POST'}) % because we're only showing matched data
% %this isn't necessary


%% Figure 6B - REVISIONS VERSION
%plot only like trn/tst combinations
% 
% trn_epoch =[1 2 3];
% tst_epoch =[1 2 3];
% gat_align ={'trn/tst:target/target'};
% 
% cond_colors = [ 0 0 1; 0 0 1; 0 0 1];
% 
% new_col = [1 1 0; .91 .41 .17; 1 0 0;];
% for pg=1:length(gat_align) % can be length 1 - target aligned recon or length 2, targ aling and dist align. for now, need only target
%     figure ('Name','trnlines;Figure6C')
%     
%     for vv = 1:length(ROIs)
%         
%         h=[];
%         for aa =1:length(trn_epoch)
% 
% 
%             % TCS: easier way to do this...
%             %ee = aa;
%             for ee =1:length(tst_epoch)
%             thisd = nan(length(subj),size(all_recons_gat{1},2));
%             for ss = 1:length(subj)
%                 
%                 subplot(size(all_recons_gat,1),length(ROIs),vv+(aa-1)*length(ROIs));hold on;
%                 
%                 thisidx = all_subj_gat==ss & all_ROIs_gat==vv & all_conds_gat(:,1)==2;
%                 thisd(ss,:) = mean(all_recons_gat{aa,ee,pg}(thisidx,:)); %aa = TRN idx, rows; ee = TST idx, col; blue = 1, red =2 , yellow =3
%                 
%             end
%             
%             % TCS: important to use "angs_gat" here - that's the actual
%             % angles used, which aren't precisely in line w/ the linspace
%             % command
%             % fixed std to use dim 1, not mode 1 (mode 0 is default)
%             my_sem = std(thisd,[],1)/sqrt((length(subj)));
%             h(ee) = plot(angs_gat,mean(thisd,1),'-','LineWidth',1,'color',new_col(ee,:));
%             
%             % TCS: updated the colors below to be consistent
%             hold on;
%             
%             plot(angs_gat,mean(thisd,1)+1.*my_sem,'-','LineWidth',.5,'color',new_col(ee,:))
%             plot(angs_gat,mean(thisd,1)-1.*my_sem,'-','LineWidth',.5,'color',new_col(ee,:))
%             
%             
%             btwn_fill = [mean(thisd,1)+1.*my_sem fliplr((mean(thisd,1)-1.*my_sem))];
%             fill([angs_gat fliplr(angs_gat)],btwn_fill,new_col(ee,:),'linestyle','none','facealpha',0.3);
%             
%             
%             clear thisidx
%             clear thisd
%             
%             
%             if ee== 1 && aa==1
%                 title(ROIs{vv});
%             else
%             end
%             
%             if ee==1 && aa==1 && vv ==1
%                 ylabel('Train PRE')
%                 set(gca,'XTick',-180:90:180,'Xticklabel',{'','',''},'TickDir','out');
%                 set(gca,'YTick',0:0.5:1.5,'Yticklabel',{'','','',''},'TickDir','out')
%                 
%             elseif ee==2 && aa==2 && vv ==1
%                 ylabel('Train DIST')
%                 set(gca,'YTick',0:0.5:1.5,'Yticklabel',{'','','',''},'TickDir','out')
%                 set(gca,'XTick',-180:90:180,'Xticklabel',{'','',''},'TickDir','out');
%             elseif ee==3 && aa==3 && vv ==1
%                 ylabel('Train POST')
%                 set(gca,'XTick',-180:90:180,'Xticklabel',{'-180','-90','0','90','180'},'TickDir','out');
%                 set(gca,'YTick',0:0.5:1.5,'Yticklabel',{'0','0.5','1.0','1.5'},'TickDir','out')
%                % ylabel('TRN/TST Matched TPTS, Recon')
%                 xlabel('Polar angle (\circ)');
%             else
%                 set(gca,'XTick',-180:90:180,'Xticklabel',{'','',''},'TickDir','out');
%                 set(gca,'YTick',0:0.5:1.5,'Yticklabel',{'','','',''},'TickDir','out')
%                 
%             end
%             end
%         end
%         
%     end
%     
% end
% 
% 
% 
% 
% set(gcf,'Position',[-132         503        2651         495])
% match_ylim(get(gcf,'Children'));
% match_xlim(get(gcf,'Children'));
% %legend(h, {'TEST PRE','DIST','POST'}) % because we're only showing matched data
% %this isn't necessary
% 
% 

%% REVISION: Figure 6E 

% plot matched trn/tst GAT fidelty data

% 
delay_tpt_range = [3.75 5.25; 8.25 9.75; 10.5 12]; 
myTR = 0.75;

delay_tpts = cell(size(delay_tpt_range,1),1);
for dd = 1:size(delay_tpt_range,1)
    delay_tpts{dd} = (tpts*myTR) >= delay_tpt_range(dd,1) & (tpts*myTR) < delay_tpt_range(dd,2);
end



plot_indiv = 1; % if 1, plot individual subj lines in line plot, otherwise, don't.
hoffset = 0.15; % +/- this much


tst_prepostepoch =[1 3];
trn_prepostepoch =[1 3];
stats_st =[];
modelcomprev =  figure('name','modelcomprev');
strk = 1:length(subj);
fig_h=[];
mrk ={'o','s'};
%for dt = 1 %:2 % differentiate btwn datatype regu v shuf 
for vv = 1:length(ROIs)
    % thisd_shuf = nan(length(subj),1000);
    for ee = 1:length(tst_prepostepoch)
        
        for aa =1:length(trn_prepostepoch)
            
            for ss = 1:length(subj)
                
                this_gat_idx = all_subj_gat==ss & all_ROIs_gat==vv & all_conds_gat(:,1)==2;
                thisd(aa,ee,ss) = mean(all_fidelity_gat{trn_prepostepoch(aa),tst_prepostepoch(ee),1}(this_gat_idx,:)); % aa = TRN idx, rows; ee = tst idx, blue = 1, red =2 , yellow =3
                stats_st = [stats_st; squeeze(thisd(aa,ee,ss)) aa ee vv  ss];
                
            end
            
            clear  this_gat_idx this_gat_shuf_idx
        end
        
        
    end
    
    
    %if dt ==1
    
    subplot(1,length(ROIs),vv)
    hold on;
    %first plot TRN PRE, TEST POST, aa = 1, ee = 2
    fig_h(1,1) = plot(1,  mean(squeeze(thisd(1,2,:))),'o','color',LORO_color,'markerfacecolor',LORO_color,'markersize',9,'linewidth',1); % CHANGING THIS FROM -- black to - black, may be confusing w previous figs!!
    my_trnpre_sem = std(squeeze(thisd(1,2,:))',[],2)/sqrt((length(subj)));
    plot([1 1],[mean(squeeze(thisd(1,2,:)))+1.*my_trnpre_sem mean(squeeze(thisd(1,2,:)))-1.*my_trnpre_sem],'-','color',LORO_color,'linewidth',1)
    hold on;
    clear my_trnpre_sem
    
    subplot(1,length(ROIs),vv)
    hold on;
    %now plot TRN POST, TEST PRE aa = 2, ee = 1
    fig_h(1,2) = plot(2, mean(thisd(2,1,:)),mrk{2},'color',LORO_color,'markerfacecolor',LORO_color,'markersize',9,'linewidth',1);
    
    my_trnpost_sem = std(squeeze(thisd(2,1,:)).',[],2)/sqrt((length(subj)));
    plot([2 2],[mean(squeeze(thisd(2,1,:)))+1.*my_trnpost_sem mean(squeeze(thisd(2,1,:)))-1.*my_trnpost_sem],'-','color',LORO_color,'linewidth',1)
    hold on;
    clear  my_trnpost_sem
    
    
    %[h,p(vv),~,stats] =
    %ttest(squeeze(thisd(1,2,:)),squeeze(thisd(2,1,:))); old
    [h,p_pre(vv),~,stats_pre] =ttest(squeeze(thisd(1,2,:)));
    t_pre(vv) = stats_pre.tstat;
    [h,p_post(vv),~,stats_post] =ttest(squeeze(thisd(2,1,:)));
    t_post(vv) = stats_post.tstat;
    %  [f_val(vv),p_val(vv)]= RMAOV1_gh([[squeeze(thisd(1,2,:));squeeze(thisd(2,1,:))], [ones(7,1);ones(7,1)+1] [strk';strk']],0.05);
    
    
    
    
    
    
    if plot_indiv == 1
                plot(1      ,squeeze(thisd(1,2,:)),'o',  'Color',LORO_color,'MarkerSize',3,'markerfacecolor',LORO_color);
        plot(2       ,squeeze(thisd(2,1,:)),'s',  'Color',LORO_color,'MarkerSize',3,'markerfacecolor',LORO_color); % 255,165,0
        
    end
    
    
    
    if vv ==1
        
        ylabel({'WM target Fidelity'})
        set(gca,'Xtick',[0 1 2 3],'Xticklabel',{'',' TRN PRE','TRN POST',''},'XTickLabelRotation',45,'TickDir','out');
        
    else
        set(gca,'Xtick',[0 1 2 3],'Xticklabel',{'','','','',''},'TickDir','out');
    end
    
    title(ROIs{vv});
    xlim([0.5 2.5])
    match_ylim(get(gcf,'Children'));
    
    
    
    
    
    
    
    
end


%% LORO non-matched TRN/TST against 1000 fidelity null 



% these need to be updated for loro 1000 stats, only 2 epochs were run 
tst_prepostepoch =[1 2]; % also changing order here 
trn_prepostepoch =[1 2];

shuf_T = nan(length(ROIs),length(trn_prepostepoch),length(trn_prepostepoch),1000);

for vv = 1:length(ROIs)
         
     
    for ee = 1:length(tst_prepostepoch)
         
        for aa =1:length(trn_prepostepoch)
            
            thisd_shuf = nan(length(subj),1000); 
            
            for ss = 1:length(subj)
                
            for ii = 1:1000
                
                this_gat_shuf_idx = all_subj_gat_shuf==ss & all_ROIs_gat_shuf==vv & all_conds_gat_shuf(:,1)==2;
                thisd_shuf(ss,ii) = mean(all_fidelity_gat_shuf{trn_prepostepoch(aa),tst_prepostepoch(ee),ii}(this_gat_shuf_idx,:));
  
            end 
                clear  this_gat_idx this_gat_shuf_idx
            end
            
            for ww =1:1000
                [~,~,~,stats] = ttest(thisd_shuf(:,ww));
                shuf_T(vv,aa,ee,ww) = stats.tstat; 
            end 
            
        end
        
    end 
    
         
           %2 * min ( mean( store_Tquart(cc,vv,dd,:) <= realT_quart(cc,vv,dd) ), mean ( store_Tquart(cc,vv,dd,:) >= realT_quart(cc,vv,dd)))
            p_pre_shuf(vv) = 2 * min(mean( squeeze(shuf_T(vv,1,2,:)) >=t_pre(vv)), mean( squeeze(shuf_T(vv,1,2,:)) <=t_pre(vv))); % TRN PRE / TST POST 
      
            p_post_shuf(vv) = 2 * min(mean( squeeze(shuf_T(vv,2,1,:)) >=t_post(vv)), mean( squeeze(shuf_T(vv,2,1,:)) <=t_post(vv))); % TRN POST / TST PRE
    
    

end 
% correct across ROIs 

fdr_thresh_pre = nan(1,size(p_pre_shuf,2));
fdr_thresh_post = nan(1,size(p_post_shuf,2));
for ii = 1:size(p_pre,1)
    [fdr_thresh_pre(ii,:), fdr_mask_pre(ii,:)] = fdr(p_pre_shuf(ii,:),0.05);
    [fdr_thresh_post(ii,:),fdr_mask_post(ii,:)] = fdr(p_post_shuf(ii,:),0.05);
end



end

