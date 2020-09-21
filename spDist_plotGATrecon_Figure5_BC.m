function spDist_plotGATrecon_fig5c(subj,sess,ROIs)
% GAT data loaded here, _GATdist_fig4TPTS, was created with the script spDist_channelRespAmp_GATdist_gh_082520
% and the tpts delay_tpt_range = [3.75 5.25; 8 9.5; 10.5 12];
root = '/share/data/spDist/';

task_dir = 'spDist';

if nargin < 1 || isempty(subj)
    subj = {'CC','KD','AY','MR','XL','SF','EK'};

    
end

if nargin < 2 || isempty(sess)

    sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
    
    
end

if nargin < 3 || isempty(ROIs)
    ROIs ={'V1','V2','V3','V3AB','hV4','LO1','IPS0','IPS1','IPS2','IPS3','sPCS'}; %ORIG
end

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


n_files= [1 2]; % how many data files do we care about?

%% load data


for yy = 1:length(n_files)
    startidx = 1;
    for ss = 1:length(subj)
        
        for vv = 1:length(ROIs)
            
            if cat_mode == 1
                % just one file to load
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
            
                else  
                    
                    fn = sprintf('%sspDist_reconstructions/%s_%s_%s_%s_%ichan%s_GATdist_fig4TPTS.mat',root,subj{ss},horzcat(sess{ss}{:}),ROIs{vv},func_suffix,nchan,vox_str);

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
                
                all_fn_gat(thisidx) = yy;
                
                startidx = thisidx(end)+1;
                
                clear data;
                
                
                
            end
            
        end
    end
    
end
end
%% plot only like trn/tst combindations 

trn_epoch =[1 2 3];
tst_epoch =[1 2 3];
gat_align ={'trn/tst:target/target'};

%cond_colors =cbrewer('qual','Set1',3);
cond_colors = [ 0 0 1; 0 0 1; 0 0 1]; 




for pg=1:length(gat_align) % can be length 1 - target aligned recon or length 2, targ aling and dist align. for now, need only target
    figure ('Name','trnlines')
    
    for vv = 1:length(ROIs)
        
        thisd = nan(length(subj),size(all_recons_gat{1},2)); %recon data
        h=[];
        for aa =1:length(trn_epoch)
            %for ee =1:length(tst_epoch)
                if aa==1
                    ee =1;
                elseif aa ==2
                    ee =2;
                elseif aa ==3
                    ee=3;
                end 
                
                thisd = nan(length(subj),size(all_recons_gat{1},2));
                for ss = 1:length(subj)
                    
                    subplot(size(all_recons_gat,1),length(ROIs),vv+(ee-1)*length(ROIs));hold on;
                    
                    thisidx = all_subj_gat==ss & all_ROIs_gat==vv & all_conds_gat(:,1)==2;
                    thisd(ss,:) = mean(all_recons_gat{aa,ee,pg}(thisidx,:)); %aa = TRN idx, rows; ee = TST idx, col; blue = 1, red =2 , yellow =3
                  
                end
                
                my_sem = std(thisd,1)/(length(subj));
                
                h(aa) = plot(linspace(-180,180,90),mean(thisd,1) - min(mean(thisd,1)),'-','LineWidth',1,'color',cond_colors(aa,:))% ,'LineWidth',2,'color',cond_colors(aa,:))
                hold on;
                %plot(linspace(-180,180,90),(mean(thisd,1) - min(mean(thisd,1)))+1.*my_sem,'-','LineWidth',.1,'color',cond_colors(aa,:),'HandleVisibility','off')
                
                %plot(linspace(-180,180,90),(mean(thisd,1) - min(mean(thisd,1)))-1.*my_sem,'-','LineWidth',.1,'color',cond_colors(aa,:),'HandleVisibility','off')
                
                btwn_fill = [(mean(thisd,1) - min(mean(thisd,1)))+1.*my_sem fliplr((mean(thisd,1) - min(mean(thisd,1)))-1.*my_sem)];     
                fill([linspace(-180,180,90) fliplr(linspace(-180,180,90))],btwn_fill,cond_colors(aa,:),'linestyle','none','facealpha',0.3);
                            
                hold on;
                line([min(xlim) max(xlim)], [0 0], 'color',[.2 .2 .2],'linewidth',0.1,'linestyle','-')
                ylim([-0.05 1.75])
                if ee== 1 && aa==1
                    title(ROIs{vv});
                else
                end
                
                if ee==1 && aa==1 && vv ==1
                    ylabel('Test Epoch 1')
                   set(gca,'XTick',-180:90:180,'Xticklabel', {'-180','-90','0','90','180'},'Xticklabelrotation',45,'TickDir','out')
                    set(gca,'YTick',0:0.5:1.5,'Yticklabel',{'0','0.5','1.0','1.5'},'TickDir','out')
                elseif ee==2 && aa==2 && vv ==1
                    ylabel('Test Epoch 2')
                     set(gca,'XTick',-180:90:180,'Xticklabel',{'-180','-90','0','90','180'},'Xticklabelrotation',45,'TickDir','out')
                    set(gca,'YTick',0:0.5:1.5,'Yticklabel',{'','','',''},'TickDir','out')
                elseif ee==3 && aa==3 && vv ==1
                    ylabel('Test Epoch 3')
                    set(gca,'XTick',-180:180:180,'Xticklabel', {'-180','0','180'},'Xticklabelrotation',45,'TickDir','out')
                    set(gca,'YTick',0:0.5:1.5,'Yticklabel',{'','','',''},'TickDir','out')
                    xlabel('Polar angle (\circ)');
                else
                 set(gca,'XTick',-180:180:180,'Xticklabel',{'','',''},'TickDir','out');
                end
                
            %end
            
        end
        
    end
    
    
    xlim([-180 180]);
    clear thisidx
    clear thisd
    
    set(gcf,'Position',[-132         503        2651         495])
    match_ylim(get(gcf,'Children'));
    match_xlim(get(gcf,'Children'));
    legend(h, {'Train 1','Train 2','Train 3'})
   % sgtitle(gat_align{pg})
    
end




%% plot matched trn/tst GAT fidelty data, on same plot, plot independently trained fidelity data 
% plot fidelity plotting TRAIN as dv 
% which tpts are we plotting throughout?

delay_tpt_range = [3.75 5.25; 8 9.5; 10.5 12];
myTR = 0.75;

delay_tpts = cell(size(delay_tpt_range,1),1);
for dd = 1:size(delay_tpt_range,1)
    delay_tpts{dd} = (tpts*myTR) >= delay_tpt_range(dd,1) & (tpts*myTR) < delay_tpt_range(dd,2);
end


trn_epoch =[1 2 3];
tst_epoch =[1 2 3];

gat_align ={'Target aligned'};
for n_files =1:2
    
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
                        end
                        
                    end
                end
                hold on;
             %  h(n_files) = plot([1 2 3], [mean(thisd(1,1,:),3) mean(thisd(2,2,:),3) mean(thisd(3,3,:),3)],'k-','linewidth',.5) % CHANGING THIS FROM -- black to - black, may be confuding w previous figs!!
                               h(n_files) = plot([1 2 3], [mean(thisd(1,1,:),3) mean(thisd(2,2,:),3) mean(thisd(3,3,:),3)],'b-','linewidth',.5) % CHANGING THIS FROM -- black to - black, may be confuding w previous figs!!

                for ii=1:3
                    my_sem = std(thisd(ii,ii,:),[],3)/(length(subj));
                    
                    hold on;
                    %plot([ii ii],[mean(thisd(ii,ii,:),3)+1.*my_sem mean(thisd(ii,ii,:),3)-1.*my_sem],'k-','linewidth',.5)
                     plot([ii ii],[mean(thisd(ii,ii,:),3)+1.*my_sem mean(thisd(ii,ii,:),3)-1.*my_sem],'b-','linewidth',.5)

                    clear my_sem
                end
                match_ylim(get(gcf,'Children'));
                
                if vv ==1
                    title(ROIs{vv})
                    xlim([0.5 3.4])
                    ylabel('WM target Fidelity')
                    set(gca,'Xtick',[0 1 2 3 4],'Xticklabel',{'','Epoch 1','Epoch 2','Epoch 3',''},'XTickLabelRotation',45,'TickDir','out');
                    
                else
                    set(gca,'Xtick',[0 1 2 3 4],'Xticklabel',{'','','','',''},'XTickLabelRotation',45,'TickDir','out');
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
                                                         % subplot(2,9,vv);hold on;

                                thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==2;
                                thisdata(dd,ss) = mean(mean(all_fidelity(thisidx,delay_tpts{dd},1),2),1);
                            end
                            
                        end
                hold on;
                 % h(n_files) = plot([1 2 3], [mean(thisdata(1,:),2) mean(thisdata(2,:),2) mean(thisdata(3,:),2)],'-','color', [0.5 0.5 0.5],'linewidth',.5) %collect n_files plot handle for legend use 
                 h(n_files) = plot([1 2 3], [mean(thisdata(1,:),2) mean(thisdata(2,:),2) mean(thisdata(3,:),2)],'-.','color', [0 0 1],'linewidth',.5) %collect n_files plot handle for legend use 

                for ii=1:3
                    my_sem = std(thisdata(ii,:),[],2)/(length(subj));
                    
                    hold on;
                   % plot([ii ii],[mean(thisdata(ii,:),2)+1.*my_sem mean(thisdata(ii,:),2)-1.*my_sem],'-','color', [0.5 0.5 0.5],'linewidth',.5)
                   plot([ii ii],[mean(thisdata(ii,:),2)+1.*my_sem mean(thisdata(ii,:),2)-1.*my_sem],'-','color', [0 0 1],'linewidth',.5)
                    clear my_sem
                end   
                        
                    clear thisidx     
                    xlim([0.5 3.5])
                    ylim([-0.05 .6])
                    if vv ==1
                        
                        ylabel('WM target Fidelity')
                        set(gca,'Xtick',[0 1 2 3 4],'Xticklabel',{'','Epoch 1','Epoch 2','Epoch 3',''},'XTickLabelRotation',45,'TickDir','out');
                        
                    else
                        
                        set(gca,'Xtick',[0 1 2 3 4],'Xticklabel',{'','','','',''},'XTickLabelRotation',45,'TickDir','out');
                    end
                    
                        
                    end
                    
         
                    
                    title(ROIs{vv});
                    
                    ylim([-0.05 .6])
                    xlim([0.5 3.5])
                    
                    
                    
                    
                end
                
        end
end     
        clear thisidx
        
        clear thisd
        set(gcf,'Renderer','painters')
        set(gcf,'Position',[-132         503        2651         495])
        legend(h, 'LORO  model', 'Ind Model')
       % sgtitle('Model Comparison')%% do this, over average delay epochs, for each condition
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

   clear thisroiidx 
   
   

%%  2-way permuation ANOVA w/_thruTime1 data 

thisy = y;
thisepoch = epoch_var;
thisroi= roi_var;
thissubj=subj_var;
iter = 20;
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

ta_anovan1 = table(ROIs',exact_store_1(:,1));
ta_anovan1.Properties.VariableNames={'ROI','Epoch'}

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
%% plot sigs on the subplots - based on anovan output
sig_colors = lines(3);
y_mod =[.2];

 
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
            text(max(xlim)-(.2*max(xlim)),max(ylim)-(y_mod(ee)*max(ylim)),sprintf('%s',sig_mrkr{ee}),'color','k','fontsize',15) %filled
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
which_effect = [3]; %3 epoch 4 cond 7epoch*cond 
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
ta_gat_anovan1.Properties.VariableNames={'ROI','Epoch'}
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

sig_colors = lines(3);
y_mod =[.3];

[p_fdr_perm_gat, p_masked_perm_gat] = fdr(exact_store_gat_1,0.05); 
     

sig_mrkr ={'+'};
y_mod = [.15];
for vv = 1:length(ROIs)
    for ee =1
        figure(modelcomp); hold on;
        subplot(1,length(ROIs),vv)
        
        
        if exact_store_gat_1(vv,ee) > p_fdr_perm_gat(ee) && exact_store_gat_1(vv,ee) <= 0.05 
            text(max(xlim)-(.2*max(xlim)),max(ylim)-(y_mod(ee)*max(ylim)),sprintf('%s',sig_mrkr{ee}),'color',[.5 .5 .5],'fontsize',15) %unfilled grey for IND
        elseif  exact_store_gat_1(vv,ee) <= p_fdr_perm_gat(ee) 
            text(max(xlim)-(.2*max(xlim)),max(ylim)-(y_mod(ee)*max(ylim)),sprintf('%s',sig_mrkr{ee}),'color','k','fontsize',15) %filled
        else
        end
        
       
    end
    
end


fprintf('the end')
end 

