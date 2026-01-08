The branch "Marginals" of this repository provides code scripts enabling stochastic scaling analysis of metabolism in biological populations based on the collapse of marginal probability distributions of the chosen metabolic proxy (_e.g._, basal metabolic rate).

To perform your analysis, download the full set of scripts and run the main function ("BMRCollapse..."). Use the function "BMRCollapseRes" to perform marginal collapses by means of the residual method (literature standard) or "BMRCollapseDist" for the distance method (newly developed approach). You can exploit the algorithm to either estimate scaling parameters from your dataset ("... estrazioneParam") or assess their sensitivity to the number of datapoints used ("... NumSamples").

Computer-generated datasets consisting of basal metabolic rates of cell-laden spheroids predicted through finite element modelling for different cell types are also provided for testing. 
