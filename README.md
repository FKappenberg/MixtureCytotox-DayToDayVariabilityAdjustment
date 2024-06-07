# Paper Adjustment for day-to-day variability in mixture models

Humans are exposed to a multitude of substances at the same time, thus, toxicity assessment not only for the individual substances but for mixtures of several substances is critical.
An important challenge in testing the mixture effects of many substances in vitro is the variability between individual experiments, also named 'day-to-day variability’.
Two steps, each consisting of independent experiments, are needed when mixtures are tested. 
In the first step, for each substance, the concentration-response relationship is measured usually with 6-10 concentration for at least 3 experiments, and EC20 values, i.e. the concentration where a fitted curve intersects with a viability value of 80%, are determined.
These values are used as reference values and the mixture of several substances is based on the individual EC20 values.
Since the mixture experiments are new, independent experiments, typically conducted on different days, day-to-day variability can lead to deviations in the cytotoxicity at the reference EC20 concentration, i.e. to higher or lower observed viabilities.
In this work, a procedure is proposed, how a single concentration per substance tested in the same experiment as the mixture can be reliably used to adjust for day-to-day variability.
In the here-established procedure, a specific additive model denoted as 'budget approach' is introduced as a reference model to explore potential positive or negative interaction effects between compounds.

This repository contains all data and code to reproduce the results from the paper and to perform the reference-curve fitting, the fitting of the mixture curves and the evaluation of the mixture curves after adjusting for day-to-day variability

## Fitting of the reference curves

- Cytotox-HepG2.RData: Contains the appropriately pre-processed viability data (i.e. normalized with respect to the mean of the experiment-wise control values, transformed to percentages)
- FittingReferenceCurves.R: Contains the code to fit concentration-response curves, including a model selection and a re-normalization step. Reference EC values are calculated as well.
- Models-CTB-HepG2.RData: Result of the code above, contains the estimated EC values, and the overall curve fit

## Fitting of the mixture curves, Evaluation 

- MixtureData.RData: Contains the mixture data, i.e. the data for the entire mixture curve, the repeated measurement at the EC20 concentration, and the EC20 concentration per substance
- FittingAnalysingMixtureCurve: This contains the entire code needed to adjust for day-to-day variability in mixture toxicity, alongside the code for fitting a mixture curve and evaluating the curve with respect to additivity/positive interaction/negative interaction

