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


spDist_behavioralAnalysis(subj,sess) 

fprintf('\nBehavioral results done (Figure 1)\n')

spDist_exampleTrial

fprintf('\n example trial done (Figure 2)\n')

% takes > 15 minutes
%spDist_plotHRFs_ERA_pRFvoxSelection(subj,sess,ROIs)

fprintf('\n univariate analysis done (Figure 3)\n')

%spDist_plotReconstructions_condAlign(subj,sess,ROIs)

fprintf('\n reconstructions generated (Figure4)\n')

%spDist_epochFidelity_statsShuf(subj,sess,ROIs)

fprintf('\n fidelity over epochs quantified (Figure 5)\n')

%spDist_modelCompare(subj,sess,ROIs)

fprintf('\n model comparisons quantified (Figure 6)\n')

spDist_neuralBehavCorr(subj,sess,ROIs)

fprintf('\n neural-behave correlation quantified (Figure 7)\n')

end

