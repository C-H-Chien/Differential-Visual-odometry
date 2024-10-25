# Differential Visual Odometry (Pose Estimation)
Internal Research @ LEMS, Brown University

## Introduction
![intro-Image](https://i.ibb.co/Fbp0bZ1/Intro-Image.png) <br />
(To be updated ...)

## Dependencies
(1) VLFeat (version 0.9.21 or higher). Check its [official webpage](https://www.vlfeat.org/) <br />
(2) Spatial Math Toolbox; already included in this repository.

## Code Usage
- Feature Tracks: In the folder ``FeatureTracks``
	- ``feature_track_veridicality.m`` generates feature tracks in a fixed window size sliding over the entire KITTI dataset. 
	- ``draw_numOfTracks_in_Lengths.m`` plots inlier ratio of feature tracks versus feature track lengths and the frequency of feature track lengths in the dataset.
- Camera Dynamic Model Fitting
- Camera Geometry Model Fitting
- Depths from Feature Tracks

## Datasets Used for the Experiments
(1) KITTI dataset: can be downloaded from [KITTI Vision Benchmark Suite](https://www.cvlibs.net/datasets/kitti/eval_odometry.php), or can be accessed through LEMS portable hard drive. <br />
(2) EuRoC dataset: can be downloaded from [Official EuRoC MAV Dataset](https://projects.asl.ethz.ch/datasets/doku.php?id=kmavvisualinertialdatasets). <br />
(3) BPOD dataset: can be accessed through LEMS protable hard drive.

## Contributors
Chiang-Heng Chien (chiang-heng_chien@brown.edu) <br />
Qiwu Zhang (qiwu_zhang@brown.edu)

