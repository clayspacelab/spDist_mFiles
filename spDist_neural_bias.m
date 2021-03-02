function spDist_neural_bias(subj,sess,ROIs)

root = spDist_loadRoot;


task_dir = 'spDist';

if nargin < 1 || isempty(subj)
   subj = {'AY','CC','EK','KD','MR','SF','XL'}; %alph


end

if nargin < 2 || isempty(sess)

    sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
  

end

if nargin < 3 || isempty(ROIs)
 ROIs = {'V1V2V3','V3AB','hV4','LO1','IPS0IPS1','IPS2IPS3','sPCS'};
end


func_suffix = 'surf';


nchan = 8;
which_vox = 0.1; % if > 1 , look for WHICH_VOXvox string; otherwise, look for VE<100*WHICH_VOX>

smooth_by = 1; % if this is 1, use regular files, otherwise, load smooth_by files

myTR = 0.75;

t_range_to_plot = [-inf 12]; % plot b/w these (s)

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

% seed random number generator 

rng(spDist_randSeed);

%% load neural data
startidx = 1;
bidx =1;
WHICH_EXCL = [13 20 21 22]; 
all_data_beh = [];
all_subj_beh= [];

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
                
                all_fidelity = nan(nblankt,size(data.recons{1},3),length(data.recons)); % timecourse of fidelity for each alignment condition
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
            all_r(thisidx,1) = (data.sess_all*100+data.r_all)';
            
            %%%inserting 
            for sessidx = 1:length(sess{ss})
              
                fn = sprintf('%s/spDist_behav_92220/%s_%s_scored.mat',root,subj{ss},sess{ss}{sessidx});
                %fprintf('Loading scored eye data from %s\n',fn);
                this_scored = load(fn);
                this_data.s_all = this_scored.ii_sess;
                
                thisbidx = bidx:(bidx+size(this_scored.ii_sess.trialinfo,1)-1);
                  
                this_data.sess_all = sessidx;
                all_data_beh= cat_struct(all_data_beh,this_data);
                all_subj_beh(thisbidx,1) = ss;
                
                bidx =thisbidx(end)+1;
            end
  
            startidx = thisidx(end)+1;
            
            clear data;
 
    end
    
end

all_data_beh.use_trial = ~cellfun( @any, cellfun( @(a) ismember(a, WHICH_EXCL), all_data_beh.s_all.excl_trial, 'UniformOutput',false));
fprintf('neural and behave loaded...\n')


%% load recons, no flipping

%all_recons_noflip = all_recons{1}; % simply reiterate, we are NOT flipping a thing 

tmprel =  all_angs(:,2) - all_angs(:,1);
this_rel = mod((tmprel+180), 360)-180;
% sign of this rel wil only match 0 bin trials, bc ang difference there is
% only the jitter, not bin dist 

flipidx = this_rel<0;
all_recons_flipped = all_recons{1}; 
all_recons_flipped(flipidx,:,:) = fliplr(all_recons_flipped(flipidx,:,:)); % updated this oct 20 


delay_tpt_range = [3.75 5.25; 8.25 9.75; 10.5 12]; %updated on nov 092020
myTR = 0.75;

delay_tpts = cell(size(delay_tpt_range,1),1);
for dd = 1:size(delay_tpt_range,1)
    delay_tpts{dd} = (tpts*myTR) >= delay_tpt_range(dd,1) & (tpts*myTR) < delay_tpt_range(dd,2);
end



cond = [1 2]; % do no distractor, then distractor 
store_b = nan(length(cond),length(ROIs),length(delay_tpts),length(subj),36);

test_store= [];

figure;

for cc= 1:length(cond)
    
    thisd = nan(length(cond),length(ROIs),length(delay_tpts),length(subj),36, length(angs));
    thisb = nan(1,length(subj));
    
   for vv = 1:length(ROIs) 
       
       for dd =1:length(delay_tpts)
           
        for ss = 1:length(subj)
  
            nruns = 0;
            
            % get one example set of trials
            ru = unique(all_r(all_subj==ss & all_ROIs==vv));
            
            % figure out how many runs we have for that subj, add to total
            nruns = nruns+length(ru);
            for nr =1:nruns
                if cc ==1 %%%%%% NOTE THIS COLLECTS THREE TRIALS IN MGS CONDITION %%%% NEED TO FIX!!!!
                    thisidx = all_data_beh.use_trial==1 &  all_r ==ru(nr) & all_subj==ss & all_ROIs==vv & all_conds(:,1)==1;
                else
                    thisidx = all_data_beh.use_trial==1 &  all_r ==ru(nr) & all_subj==ss & all_ROIs==vv & all_conds(:,1)==2 & all_conds(:,6)==0 ;
                end
                
                thisd(cc,vv,dd,ss,nr,:) = squeeze(mean(mean(all_recons_flipped(thisidx,:,delay_tpts{dd}),1),3));
                
                thisb(ss) = atan2d(sum(mean(mean(all_recons_flipped(thisidx,:, delay_tpts{dd}),1),3).*sind(angs)),...
                    sum(mean(mean(all_recons_flipped(thisidx,:,delay_tpts{dd}),1),3).*cosd(angs)));
                
                store_b(cc,vv,dd,ss,nr) = thisb(ss);
                test_store =[test_store; thisb(ss) dd vv cc ss]; 
                
           clear thisidx
           clear thisb
            end
        end
        end
         
        
    end
 
end

% do stats
clear thisidx the_y

for ss =1:length(subj)
thisidx = test_store(:,2)==3 & test_store(:,4)==2 & test_store(:,3)==1 & ~isnan(test_store(:,1)) & test_store(:,5)==ss ;
nb_mean(ss) = mean(test_store(thisidx,1));
clear thisidx 
end 

[h,p] = ttest(nb_mean); 
fprintf('neural bias computed ...\n')
end

