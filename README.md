# hybrid-sensor2023-scripts

This repository contains Matlab scripts for generating the figures in the manuscript

Yu Mao, Harry Dankowicz (2023) *On a Principle for Mass Sensing Using Self-Excited Template Dynamics of Coupled Oscillators and Root-Finding Algorithms*, Journal of Sound and Vibration


## Content
This repository contains:

- `coco_2020Mar22`: March 22, 2020, release of COCO 
  [sourceforge.net/cocotools](https://sourceforge.net/projects/cocotools/)
  used for demos
- `add_recording.m`: COCO-compatible utility for recording Newton iterates
- `hybridsensor_delay.m`: COCO-compatible zero problem for actuator with constant delay
- `hybridsensor_PD.m`: COCO-compatible zero problem for actuator with feedback
- `fourier_coeff.m`: function for extracting Fourier coefficients from periodic signal
- `fourier_func.m`: function for reconstructing periodic signal from Fourier coefficients
- `Figure2.m`: script to generate the panels in Figure 2 describing variations
  in the amplitude ratio against mass ratio for the passive sensor design 
  in Section 2.1
- `Figure4.m`: script to generate Figure 4 describing variations
  in the feedback parameter against mass ratio for the active sensor design 
  in Section 2.2
- `Figure5.m`: script to generate the panels of Figure 5 showing projections of
  Hopf bifurcation curves for the active sensor design in Section 2.2
- `Figure6.m`: script to generate the panels of Figure 6 showing projections of
  Hopf bifurcation curves for the active sensor design in Section 2.2
- `Figures8and9.m`: script to generate the panels of Figures 8 and 9 describing
  variations in the amplitude ratio against mass ratio for the nonlinear 
  sensor design in Section 3
- `Figure10.m`: script to generate Figure 10 describing variations in the 
  amplitude ratio against mass ratio for the nonlinear sensor design
  in Section 3
- `Figure12.m`: script to generate the panels of Figure 12 describing
  iterates of a Newton algorithm applied to the hybrid realization of the
  nonlinear sensor design in the case of a constant-delay actuator
  considered in Section 4.2
- `Figure14.m`: script to generate the panels of Figure 14 describing
  iterates of a Newton algorithm applied to the hybrid realization of the
  nonlinear sensor design in the case of an actuator with feedback control
  considered in Section 4.3
- `Figure15a.m`: script to generate the left panel of Figure 15 showing
  convergence of a Newton algorithm applied to the hybrid realization of the
  nonlinear sensor design in the case of a constant-delay actuator
  considered in Section 4.2
- `Figure15b.m`: script to generate the right panel of Figure 15 showing
  convergence of a Newton algorithm applied to the hybrid realization of the
  nonlinear sensor design in the case of an actuator with feedback control
  considered in Section 4.3


## Usage

1. execute script `startup.m` in folder `coco_2020Mar22/coco` to initialize `coco`
2. execute each of the scripts `Figure*.m` to generate the corresponding figure(s)
