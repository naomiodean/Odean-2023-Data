# Data Repository for "Transient oscillations of neural firing rate associated with routing of evidence in a perceptual decision"

This repository contains data from the paper by Odean, Sanayei, and Shadlen, published in the Journal of Neuroscience in 2023. It includes the following MATLAB files:

- `SampleCode.m`: A MATLAB script providing examples of how to load and analyze the data from `cuedAttentionData.mat` and `variableLocationData.mat`. This is the best place to start.
- `cuedAttentionData.mat`: Contains spike data from a cued attention task described in the study. Data are formatted where each row represents a trial, and columns contain task parameters and spike times.
- `variableLocationData.mat`: Contains spike data from a variable location task. Each row corresponds to a trial with columns detailing task parameters and spike times.
- `BigSBigT.m`: A MATLAB function used by `SampleCode.m` to generate plots of Peri-Stimulus Time Histograms (PSTHs).

