# Low-dose-CDI
# Supplementary Data Codes 

**Computational Microscopy beyond Perfect Lenses**

Xingyuan Lu<sup>1,2*</sup>, Minh Pham<sup>1,3*</sup>, Elisa Negrini<sup>3</sup>, Damek Davis<sup>4</sup>, Stanley J. Osher<sup>3</sup> & Jianwei Miao<sup>1†</sup>    

*<sup>1</sup>Department of Physics & Astronomy and California NanoSystems Institute, University of California, Los Angeles, CA 90095, USA.*    
*<sup>2</sup>School of Physical Science and Technology, Soochow University, Suzhou 215006, China.*    
*<sup>3</sup>Department of Mathematics, University of California, Los Angeles, CA 90095, USA.*      
*<sup>4</sup>School of Operations Research and Information Engineering, Cornell University. Ithaca, NY 14850, USA.*   
**These authors contributed equally to this work.*     
*†Correspondence and requests for materials should be addressed to J.M. (miao@physics.ucla.edu).*  

## Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Repositary Contents](#repositary-contents)

# Overview

Lens-based microscopy has played an important role in the evolution of modern science and technology. The development of coherent diffractive imaging (CDI) in 1999 that replaces the lens with coherent illumination and computational algorithms has transformed our conventional view of microscopy. Here we demonstrate through mathematical analysis and numerical experiments that in situ CDI, which harnesses the coherent interference between a strong beam illuminating a static structure and a weak beam illuminating a dynamic structure, could be the most dose-efficient imaging method. At low doses, in situ CDI can achieve higher resolution than perfect lenses with the point spread function as a delta function. We expect that computational microscopy based on in situ CDI can be implemented in different imaging modalities with photons and electrons for low-dose imaging of radiation-sensitive materials and biological samples.

# System Requirements

We recommend a computer with 16G DRAM, standard i7 4-core CPU, and a GPU to run most data analysis source codes.

## Software Requirements

### OS Requirements

This package has been tested on Windows 11.   

### Matlab Version Requirements

This package has been tested with `Matlab` R2022b. All the codes have to run in their own folders. We recommend the use of `Matlab` version R2021a or higher to test the data and source codes.

# Repositary Contents

### 1. Simulations of Perfect lens imaging and Insitu CDI for different dose.

Folder: [Perfect lens and insitu CDI](./Fig3_Perfectlens_and_iCDI)

This folder contains the main simulation code and subfuctions needed. 
Step-1: run the code 'step1_PLs_and_iCDI_dp_sim_data300nm_QE08.m' to creat data for perfect lens imaging and diffraction patterns of insitu CDI; One can change thickness in line 18;
Step-2: run the code 'step2_refinelacey_300nm.m' to do the insitu CDI reconstruction. Focus on the refinement of the reference structure in this step; 
Step-3: run the code 'step3_finalrecons_300nm.m' to do the final reconstruction.

### 2. Simulations of Perfect lens imaging and Low Dose CDI Scanning for different dose.

Folder: [Perfect lens and low dose CDI scanning](./Fig4_Perfectlens_and_LoCDI)

This folder contains the main simulation code and subfuctions needed. 
Step-1: run the code 'main1_ptycho_simdp.m' to creat data for perfect lens imaging and diffraction patterns of insitu CDI; 
Step-2: run the code 'main2_ptycho_recons.m' to do the low dose scanning CDI reconstruction; 
Step-3: run the code 'main3_ptycho_FRC.m' to plot the Fourier Ring correlation curves.
