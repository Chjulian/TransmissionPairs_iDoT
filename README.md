# HIV Direction of Transmission

This is the code and data repository for the HIV Direction of Transmission project.

This README includes:

- How to cite the paper
- An overview of repository structure


## Please remember to check the following:
- Villabona-Arenas CJ, Hué S, Baxter JAC, Hall M, Lythgoe KA, Bradley J, Atkins KE (2022) "**Using Phylogenetics to infer HIV-1 Transmission Direction between known Transmission Pairs**", https://doi.org/10.1101/2021.05.12.21256968 [*To be Updated*]

## Repository overview and organisation
Our repository is structured in three parts: data, src and models.

### 1. Data.
This is a csv (`pairdata.csv`) with the following fields:

1. **LANLdb_cluster_name**: Cluster name collected from the Los Alamos HIV Sequence Database.
2. **risk**: Risk group for the transmission pair.
3. **status.at.transmission.source**: HIV infection status of the transmitter at transmission to the recipient.
4. **sampling.calc.source.mid**: (days) Sampling of the transmitter partner relative to transmission time.
5. **sampling.calc.recipient.mid**: (days) Sampling of the recipient partner relative to transmission time.
6. **alignment.length**: (nucleotides) Sequence alignment length.
7. **source.nuc.div**: (substitutions per site) Nucleotide diversity of the transmitter virus population.
8. **rec.nuc.div**: (substitutions per site) Nucleotide diversity of the recipient virus population.
9. **hap.n.source**: Number of unique HIV-1 sequences collected from the transmitter partner.
10. **hap.n.rec**: Number of unique HIV-1 sequences collected from the transmitter partner.
11. **pFr1**: Probability of one founder unique sequence in the recipient. Estimated using the approach of Villabona-Arenas CJ, Hall M, Lythgoe KA, Gaffney SG, Regoes RR, Hué S, Atkins KE (2020) "**Number of HIV-1 founder variants is determined by the recency of the source partner infection**", Science, https**://doi.org/10.1126/science.aba5443.
12. **topology**: The topology class (either Paraphyletic-polyphyletic (PP), Paraphyletic-monophyletic (PM) or Monophyletic-monophyletic (MM))
13. **p.ancestral.pars**: The Ancestral State probability of the transmitter calculated using Parsimony methods.
14. **p.ancestral.ER**: The Ancestral State probability of the transmitter calculated using Maximum Likelihood methods.
15. **PD.UEH.source**: The phylogenetic diversity (PD, sensu Faith 1992) of the transmitter partner' subtree.
16. **PD.UEH.recipient**: The phylogenetic diversity (PD, sensu Faith 1992) of the recipient partner' subtree.
17. **tree.MinRootToTip.source**:  (substitutions per site) The minimum root-to-tip distance of the transmitter partner' taxa.
18. **tree.MinRootToTip.recipient**:  (substitutions per site) The minimum root-to-tip distance of the recipient partner' taxa.
19. **identities**: The identity of the tip(s) that minimizes the number of internal nodes along the paths between itself and the root (either Transmitter, Recipient or Both).
20. **min.inter.TBL**: (substitutions per site) The shortest patristic distance between tips from the partners.
21. **sacc**: GenBank accession of the selected hit from BLAST (related genetic sequence).
22. **pident**: The Blast' percentage of identical matches.
23. **site.model**: The best-fit model of nucleotide substitution using ModelFinder.
 
### 2. Source Code (src)
`BinomialLasso.R` is the script to fit *binomial* statistical models using Lasso Regression and infer the AUC (Area Under The Curve) ROC (Receiver Operating Characteristics) curve. 

`MultinomialLasso.R` is the script to fit *multi-categorical ordinal* statistical models using Lasso Regression and infer the AUC (Area Under The Curve) ROC (Receiver Operating Characteristics) curve. 

`Figures.R` generates the figures and supplementary figures in the manuscript.
 
### 3. Models
The model data generated with `BinomialLasso.R` or `MultinomialLasso.R` is stored here. The following files are saved for each approach—that is, binomial or multi-class ordinal; Probabilistic or via Maximum Parsimony; using the identity of the transmitter or not, when defining the covariates):

- [filename]**\_matriX.rds**: A list with the response matrices coded as either binomial or ordinal (multi-categorical) for each model.
- [filename]**\_matriY.rds**: A list with the predictor matrices for each model.
- [filename]**\_models.rds**: a list with each of the models with the class "lognet", "glmnet" or "ordinalNet".
- [filename]**\_rocs.rds**: a list with the ROC (Receiver Operating Characteristics) for each model.
- [filename]**\.csv**: a data frame with the following information for each model:
         
         - type: The model class (either E,S,G,P or any combination of them).
         - non0Terms: The model terms with non-zero coefficients.
         - thresholds: The threshold used for the definition of the response.
         - auc: The AUC (Area Under The Curve).
         - auc.min: The lower confidence interval (CI) of the AUC.
         - auc.max: The upper confidence interval (CI) of the AUC.
         - lamdba: The value of the tuning parameter lambda.
- [filename]**\.pdf**: Plots of the AUROC (Area Under the Receiver Operating Characteristics)
