function spDist_run_all(subj,sess,ROIs)

% set-up default arguments 
if nargin < 1 || isempty(subj)
    subj = {'AY','CC','EK','KD','MR','SF','XL'};
end

if nargin < 2 || isempty(sess)
    
    sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
    
    
end

if nargin < 3 || isempty(ROIs)
    ROIs = {'V1V2V3','V3AB','hV4','LO1','IPS0IPS1','IPS2IPS3','sPCS'};
end
% takes ~45 minutes in total 
fprintf('\n Starting behavioral analysis .. @ %s \n',datestr(now,'HH:MM:SS'))
spDist_behavioralAnalysis(subj,sess) % takes ~ 1 minute
fprintf('\n Behavioral analysis complete @ %s \n',datestr(now,'HH:MM:SS'))

fprintf('\n Behavioral results done (Figure 1)\n')

fprintf('\n Beginning univariate analysis (Figure 2) @ %s \n',datestr(now,'HH:MM:SS'))
spDist_plotHRFs_ERA_pRFvoxSelection_btwnRFstats(subj,sess,ROIs) % extraneous figs are produced; takes ~ 5 minutes 
fprintf('\n Univariate analysis complete (Figure 2) @ %s \n',datestr(now,'HH:MM:SS'))

fprintf('\n Creating example trial (Figure 3) @ %s \n',datestr(now,'HH:MM:SS'))
spDist_plotBasis
spDist_exampleTrial %takes < 1 minute
fprintf('\n Example trial done (Figure 3) @ %s \n',datestr(now,'HH:MM:SS'))


fprintf('\n Starting reconstructions (Figure4) @ %s \n',datestr(now,'HH:MM:SS'))
spDist_plotReconstructions_condAlign(subj,sess,ROIs) % takes ~1 minutes
spDist_fidelity_stats_shuf(subj,sess,ROIs) % takes ~15 minutes 
fprintf('\n Reconstructions generated (Figure4) @ %s \n',datestr(now,'HH:MM:SS'))

fprintf('\n Starting Fidelity over epochs quantification (Figure 5) @ %s \n',datestr(now,'HH:MM:SS'))
spDist_epochFidelity_statsShuf(subj,sess,ROIs) % takes ~20 minutes 
fprintf('\n Fidelity over epochs quantified (Figure 5) @ %s \n',datestr(now,'HH:MM:SS'))

fprintf('Starting model comparison  (Figure 6) @ %s \n',datestr(now,'HH:MM:SS'))
spDist_plotGAT1
spDist_modelCompare(subj,sess,ROIs)
fprintf('Model comparison complete (Figure 6) @ %s \n',datestr(now,'HH:MM:SS'))

fprintf('Starting neural-behave correlation quantification (Figure 7) @ %s \n',datestr(now,'HH:MM:SS'))
spDist_neuralBehavCorr(subj,sess,ROIs) % takes ~3 minutes
fprintf('Neural-behave correlation complete (Figure 7) @ %s \n',datestr(now,'HH:MM:SS'))


end

