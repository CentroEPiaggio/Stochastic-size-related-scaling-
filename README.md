The branch "Marginals" of this repository provides code scripts enabling stochastic scaling analysis of metabolism in biological populations based on the collapse of marginal probability distributions of the chosen metabolic proxy (_e.g._, basal metabolic rate).

To perform your analysis, download the full set of scripts and run the main function ("BMRCollapse..."). Use the function "BMRCollapseRes" to perform marginal collapses by means of the residual method (literature standard) or "BMRCollapseDist" for the distance method (newly developed approach). You can exploit the algorithm to either estimate scaling parameters from your dataset ("... estrazioneParam") or assess their sensitivity to the number of datapoints used ("... NumSamples").
The main function also allows for the evaluation of the statistical significance of the obtained collapse (_i.e._, whether collapsed marginal distributions are statistically indistinguishable) based on the Anderson-Darling test.

Computer-generated datasets consisting of basal metabolic rates of cell-laden spheroids predicted through finite element modelling for different cell types are also provided for testing.

------------------------------------------------------

The branch "Joints" of this repository provides code scripts enabling stochastic scaling analysis of metabolism in biological populations based on the collapse of joint probability distributions of the chosen size and metabolic proxies (_e.g._, mass and basal metabolic rate).

To perform your analysis, download the full set of scripts and run the main function ("MultiParametricFittingDist..."). This function performs joint collapses by means of the the distance method (first approach for multi-dimensional collapse reported so far). You can exploit the algorithm to either estimate scaling parameters from your dataset ("... estrazioneParam") or assess their sensitivity to the number of datapoints used ("... NumSamples").
To evaluate the statistical significance of the obtained collapse (_i.e._, whether collapsed joint distributions are statistically indistinguishable) based on adapted 3-ways Anova, run the function "statMain".

Computer-generated datasets consisting of masses and basal metabolic rates of cell-laden spheroids predicted through finite element modelling for different cell types are also provided for testing.
