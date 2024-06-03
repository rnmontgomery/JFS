# JFS
Contains code used in the paper "Jointly assessing multiple endpoints in pilot and feasibility studies". This repository contains R functions to propsectively define a cut point using the JFS method (JFScut.R) and code (JFS_N.R) to approximate the required sample size to achieve a specified $P(P \mid F)$ and $P(P \mid I)$

Here is a brief example to use the JFScut function to estimate the required sample size for the motivating example, which includes recruitment rate, to achieve  $P(P \mid F) = 0.80$ and $P(P \mid I) = 0.05$.


1) Use max sample size, e.g., N = 200

JFScut(F =  c(5.75, 0.8),
       I = c(4.44, 0.775),
       nsims = 1000,
       N = 200, 
       PpF = 0.80,
       nprelim = 165,
       recruitment = TRUE,
       alpha = c(1),
       beta = c(1),
       alpha_rec = 0.01,
       beta_rec = 0.01,
       recwindow = 36)
       
PpF will always be met, so we focus on PpI. With N = 200, it was 0.003.

2) OCs were better than reqquired, decrease sample to N = 100

JFScut(F =  c(5.75, 0.8),
       I = c(4.44, 0.775),
       nsims = 1000,
       N = 100, 
       PpF = 0.80,
       nprelim = 165,
       recruitment = TRUE,
       alpha = c(1),
       beta = c(1),
       alpha_rec = 0.01,
       beta_rec = 0.01,
       recwindow = 36)

 With N = 200, PpI = 0.054, slightly increase sample size.

 3) OCs were slighlty worse than required, increase sample, N = 105
 
 JFScut(F =  c(5.75, 0.8),
        I = c(4.44, 0.775),
        nsims = 1000,
        N = 105, 
        PpF = 0.80,
        nprelim = 165,
        recruitment = TRUE,
        alpha = c(1),
        beta = c(1),
        alpha_rec = 0.01,
        beta_rec = 0.01,
        recwindow = 36)
 
 With N = 105, PpI = 0.045, so approximately N = 105 participants are needed.
 This can be confirmed with a larger run (e.g., nsims = 10000).
