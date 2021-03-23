# spDist_mFiles
The MIT License (MIT)
Copyright (c) 2021 Grace E. Hallenbeck, Thomas C. Sprague, Masih Rahmati, Kartik Sreenivasan, & Clayton E. Curtis

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

dependencies
within spDist_mFiles:
RMAOV1_gh.m
RMAOV2_gh.m
RMAOV33_gh.m
RMAOV1.m
RMAOV2.m
RMAOV33.m
spDist_condColors.m
spDist_randSeed.m
spDist_concatBehav.m
spDist_loadRoot (should be modified to point to where neural data is stored)
to recreate the Figures shown in Hallenbeck et al, 2021, run the following code
dependencies included in spDist_mFiles

software:
Matlab2018b (or higher, though untested outside Matlab2018b)

hardware: n/a
install time: download

Instructions & Expected output
*PLEASE NOTE, at this time (03/02/2021), each script produces additional plots.
Upon acceptance, these will be relocated. The plots used to create figures are noted in the figure titles.
Expected run-time: 45 minutes

1. spDist_behavioralAnalysis (Figure 1)
Sections denote portions of Figure 1. Ensure behavioral data is accessible by 'root'. Eye data shared herein has been extracted, preprocessed, and scored by iEye_ts (github/tommysprague/iEye_ts). To plot aligned eye traces (1B) and scatter plot (1C; takes > 10 mins), set scatterplot_1BC = 1. Supplementary Figure 1 data is also given by this script.
2. spDist_plotHRFs_ERA_pRFvoxSelection_btwnRFstats (Figure 2)
Change spDist_loadRoot to location where extracted neural data has been stored locally. 4th input argument, 'alignment' defaults to target position. To recreate Figure 2c (distractor aligned Rf in vs out), set alignment argument to 'dist_ang_all'.
3. spDist_plotBasis,  spDist_exampleTrial.m (Figure 3)
3a is a screenshot from a retinotopy sesssion in SUMA. plotBasis recreates 3b (left), _exampleTrial creates  3b (right).
4. spDist_plotReconstructions_condAlign. Will plot 3 separate plots for Figure 4 A (distractor absent, target aligned) B (distractor present, target aligned), & C (distractor present, distractor aligned)(Figure 4). To obtain statistical results examining fidelity on a TR basis, compared to a null-distribution, run spDist_fidelity_stats_shuf.m. Requires _reconthruTime1_shuff1000.mat (in spDist_reconstructions dir) (Figure 4DE). takes ~15 minutes
5. spDist_epochFidelity_statsShuf.m (Figure 5)
6. spDist_plotGaAT1 & spDist_modelCompare (Figure 6)
spDist_plotGAT1 uses GAT/LORO model with name _GATdist_gh(Figure 6B)
relies on GAT trained data generated by spDist_channelRespAmp_GATdist, saved with _GATdist_gh (IMPORTANT spDist_channelRespAmp_GATdist.m can take up to 24 hours to run.)
7. spDist_neuralBehavCorr (Figure 7)
