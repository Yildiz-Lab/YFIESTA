# Bayesian 2-D Step-Fitting Package

Run `runMoleculeBayes2DPipeline.m` after loading `Molecule` in MATLAB.
It uses the same Molecule input, QC/EFO truncation, plotting, FIESTA export,
and FIONA export workflow as the SIC package, but replaces the greedy SIC
splitter with `BayesStepFinder2D`.

The Bayesian finder computes a globally optimal MAP segmentation by dynamic
programming. A shared, path-oriented noise covariance makes transverse noise
less likely to be mistaken for a genuine step. `ExpectedDwellPoints` is the
main conservativeness control: increase it to favor fewer, longer dwells.

Outputs are prefixed `Molecule_Bayesian2D` and FIONA files are written to
`FIONA_Bayesian2D/`.
