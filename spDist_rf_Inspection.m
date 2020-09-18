function spDist_rf_Inspection(subj,sess,ROIs)


root = '/share/data/spDist/'

task_dir = 'spDist';

if nargin < 1 || isempty(subj)
    subj = {'CC','KD','AY','MR','XL','EK','SF'};
end

if nargin < 2 || isempty(sess)
    
    sess = {'spDist1','spDist2'}; %{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
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

t_range_to_plot = [-inf 12]; % plot b/w these (s)

trn_tpts = 7:15; % if blank, load files w/ no _trn%ito%i, otherwise,

plot_tpts = 7:15; % for a plot where we average reconstructions over a fixed time window

dist_time = 4.5; % onset at 4 s


% set up file loading strings for below
% if smooth_by == 1
%     smooth_str = '';
% else
%     smooth_str = sprintf('_smooth%i',smooth_by);
% end

% if isempty(trn_tpts)
%     trn_str = '';
% else
%     trn_str = sprintf('_trn%ito%i',trn_tpts(1),trn_tpts(end));
% end

% if which_vox < 1
%     vox_str = sprintf('_VE%03.f',100*which_vox);
% else
%     vox_str = sprintf('_%ivox',which_vox);
% end



%
% % for fidelity timecourses
% tmpcolors = lines(7);
%
% ROI_colors = [repmat(tmpcolors(1,:),3,1); % V1, V2, V3
%     tmpcolors(4,:);             % V3AB
%     tmpcolors(1,:);             % hV4
%     repmat(tmpcolors(3,:),1,1); % VO1
%     repmat(tmpcolors(1,:),2,1); % LO1/2
%
%     repmat(tmpcolors(2,:),2,1); % TO1-2
%
%     repmat(tmpcolors(5,:),2,1); % IPS0-1
%     repmat(tmpcolors(6,:),2,1); % IPS2-3
%     tmpcolors(7,:);             % sPCS
%     ];             % iPCS (color 1...)
%

%% load data
startidx = 1;
for ss = 1:length(subj)
    
    for vv = 1:length(ROIs)
        for sh =1:length(sess)
            if cat_mode == 1
                % just one file to load
                fn = sprintf('%s/%s_trialData/%s_%s_%s_%s_trialData.mat',root,task_dir,subj{ss},sess{sh},ROIs{vv},func_suffix);
                
                fprintf('loading %s...\n',fn);
                data = load(fn);
                
                
                if vv == 1 && ss == 1 && sh ==1
                    % initialize variables...
          
                    nblankt = size(data.rf.exp,2);
                    %  all_recons = cell(size(data.recons));
                    
                    all_exp = nan(nblankt,size(data.rf.exp,1));
                    all_ve = nan(nblankt,size(data.rf.exp,1));
                    all_x0 = nan(nblankt,size(data.rf.exp,1));
                    all_y0 = nan(nblankt,size(data.rf.exp,1));
                    
                    all_sigma = nan(nblankt,size(data.rf.exp,1));
                    
                    all_subj = nan(nblankt,1);
                    all_ROIs = nan(nblankt,1);
                    all_sess = nan(nblankt,1);
                    
                    
                    
                    
                    % ugh have to do this in a multi-D array...
                    %all_r2 = nan(length(ROIs),length(tpts),length(subj));
                    
                end
                
                
                
                thisidx = startidx:(startidx+size(data.rf.exp,2)-1);
                
                all_exp(thisidx,:) = data.rf.exp';
                all_ve(thisidx,:) = data.rf.ve';
                all_x0(thisidx,:) = data.rf.x0';
                all_y0(thisidx,:) = data.rf.y0';
                
                all_sigma(thisidx,:) = data.rf.sigma';
                
                all_subj(thisidx) = ss;
                
                
                all_ROIs(thisidx) = vv;
                
                all_sess(thisidx) = sh;
                
                
                startidx = thisidx(end)+1;
                
                clear data;
                
            else
                % NOT SUPPORTED YET!!!!
                
                %             for sess_idx = 1:length(sess{ss})
                %                 % build fn
                %                 fn = sprintf('%swmChoose_reconstructions/%s_%s_%s_%s_%ichan%s%s%s_recon_cv_thruTime1.mat',root,subj{ss},sess{ss}{sess_idx},ROIs{vv},func_suffix,nchan,vox_str,smooth_str,trn_str);
                %
                %                 fprintf('loading %s...\n',fn);
                %                 data = load(fn);
                %
                %
                %                 if vv == 1 && ss == 1
                %                     % initialize variables...
                %
                %
                %                     nblankt = length(ROIs)*numel(sess)*size(data.recons,1);
                %
                %                     all_recons = nan(nblankt,size(data.recons,2),size(data.recons,3));
                %                     all_conds = nan(nblankt,size(data.c_map,2));
                %
                %                     all_fidelity = nan(nblankt,size(data.recons,3)); % timecoruse of fidelity
                %
                %
                %                     all_subj = nan(nblankt,1);
                %                     all_ROIs = nan(nblankt,1);
                %                     all_sess = nan(nblankt,1);
                %
                %                     angs = data.angs;
                %                     tpts = data.delay_tpts;
                %
                %                     all_r2 = nan(length(ROIs),length(tpts),length(subj));
                %
                %                 end
                %
                %                 % set up our variable used to compute R2
                %                 if sess_idx == 1
                %                     tmp_r2 = nan(length(tpts),length(sess{ss})); % average acorss sessions...
                %                 end
                % %
                %                 thisidx = startidx:(startidx+size(data.rf.exp,1)-1);
                %
                %
                %                 all_recons(thisidx,:,:) = data.recons;
                %                 all_fidelity(thisidx,:) = squeeze(mean(cosd(angs) .* data.recons,2));
                %
                %                 all_conds(thisidx,:) = data.c_map;
                %
                %
                %
                %                 all_subj(thisidx) = ss;
                %
                %
                %                 all_ROIs(thisidx) = vv;
                %
                %                 all_sess(thisidx) = sess_idx;
                %
                %                 tmp_r2(:,sess_idx) = squeeze(mean(mean(data.r2_all,1),2));
                %
                %                 startidx = thisidx(end)+1;
                %
                %                 clear data;
                
            end
            
            %all_r2(vv,:,ss) = mean(tmp_r2,2);
            %clear tmp_r2;
            
            
        end
    end
end
%% plot and inspect bu subj

for ss = 1:length(subj)
    figure
    for vv = 1:length(ROIs)
        subplot(1,length(ROIs),vv); hold on
        %for sh =1:length(sess)
           thisidx = all_subj ==ss & all_ROIs ==vv ;
           histogram(all_exp(thisidx,:))
        %end
    end
    if vv==1    
    ylabel('Frequency');
    xlabel('Value')
    else
    end
    match_ylim(get(gcf,'Children'));
    match_xlim(get(gcf,'Children'));

    %sgtitle('%s',subj{ss})
    
end

%% for ss = 1:length(subj)
    figure
    for vv = 1:length(ROIs)
        subplot(1,length(ROIs),vv); hold on
      
           thisidx =  all_ROIs ==vv & all_ve > .1 ;
           histogram(all_exp(thisidx,:))
           n_vox = sum(all_exp(thisidx,:) > 1);
           n_pro = (sum(all_exp(thisidx,:) > 1) / length(all_exp(thisidx,:)))*100;
           text(max(xlim)-(.5*max(xlim)),max(ylim)-(.25*max(ylim)),sprintf('total %i',n_vox),'FontSize',12)
            text(max(xlim)-(.5*max(xlim)),max(ylim)-(.3*max(ylim)),sprintf('pro %.2f',n_pro),'FontSize',12)
 
    if vv==1    
    ylabel('Frequency');
    xlabel('Value')
    else
    end
    title(ROIs{vv})
   % match_ylim(get(gcf,'Children'));
    match_xlim(get(gcf,'Children'));

    clear thisidx
    
    end
    
  sgtitle('All subj RF-fit exp value per ROI')

end



