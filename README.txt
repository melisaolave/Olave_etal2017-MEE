DESCRIPTION:

This approach involves inferring a distribution of expected extra lineages given a model and the fit with real data (i.e. gene trees) to evaluate the likelihood of the different proposed models. 
The original article describing the method:
Olave, M., Avila, L.J., Sites, J.W. Jr and Morando, M (2017). Detecting Hybridization by Likelihood Calculation of Gene Tree Extra Lineages given Explicit Models. Methodos in Ecology and Evolution

The sequence of steps from (i) to (iii) is shown in (Fig. 2; Olave et al. in press) and corresponds to:
(i) Proposing a model including a species tree topology and branch length in coalescent units (CU = t/Ne; where t is generations and Ne is the effective population size). Tip: Species tree programs like ASTRAL (Mirarab et al. 2014) and MP-EST (Liu et al. 2010) will now calculate the species tree in CU.
(ii) Inference of expected extra lineages given the model. A set of H gene trees is simulated following the model described in (i). Extra lineages are counted between each simulated gene tree and the species tree S proposed in (i), providing a distribution of the expected number of extra lineages under the model in (i). The method currently uses the software ms (Hudson 2002). Then the package Phylonet (Than and Nakhleh 2009) is used to count extra lineages, and provide the distribution of expected discordance.
(iii) Likelihood of the empirical data. The empirical gene trees are compared to the distribution generated in (ii) and the likelihood calculations are obtained.
Steps i - iii should be repeated with as many models as is desired.
(iv) Model selection. A likelihood ratio test is then used to select the best model among the proposed ones.

REQUIREMENTS
functions.R --> github.com/melisaolave/Olave_etal2017-MEE
ape library installed in R
ms (Hudson 2002) located at the same folder than the data. https://uchicago.app.box.com/s/l3e5uf13tikfjm7e1il1eujitlsjdx13
Phylonet (Than and Nakhleh 2009). https://bioinfocs.rice.edu/PhyloNet

EXECUTION
See the tutorial R script (tutorial.R) for an example of how to use the method.

COMMENTS
This method uses most of its time estimating the probability vector, as described in (iii), thus time is very dependent on the number of H gene trees simulated. The method description calculates the distribution using H = 10,000, and it is recommended to use a similar number, specially when very similar models wants to be compared. However, in order to reduce the time consumption, it is possible to estimate the distribution using a small set of simulated gene trees for estimated the p vector (e.g. H = 5,000). Note that this might decrease the method accuracy, and will only be correct if models proposed are very different.

RECOMMENDATIONS BEFORE USING
1. Note that this method is designed to be applied only to a small number of species, and we strongly recommend restricting the number of outgroups.
If possible, we recommend that users incorporate the strategy described by Olave et al. (in press) for real data analyses (section 2.3), where individuals with a strong hybridization signature are removed prior to the species tree estimation, and later incorporated for preforming this test. This will prevent high underestimation of species divergence, and then results are more likely to be accurate.

2. If there is missing data (an individual is not present a gene tree), the distribution of expected extra lineages for those gene trees with missing data should be estimated separated for this case (without this individual) and compare to this particular gene tree. 
You should have some R program skills for doing this. Note that this will substantially increase the time consumption (since multiple distributions need to be estimated). Note that the current example does not accept missing data.


Any question, please refer to:
Email: melisa.olave@uni-konstanz.de





