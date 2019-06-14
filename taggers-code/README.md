# Jets @ Les Houches 2019

## Basic code for q/g taggers

The file [qg-taggers.hh](qg-taggers.hh) contains code for 3 classifiers:

 - Angularity(alpha, R, optional_groomer):

    - this first grooms the jet if the groomer is set

    - the shape is computed as
    
      ![\Large \frac{1}{\sum_i p_{ti} R^\alpha}\sum_i p_{ti}\theta_i^\alpha](https://latex.codecogs.com/svg.latex?%5CLarge%20%5Cfrac%7B1%7D%7B%5Csum_i%20p_%7Bti%7D%20R%5E%5Calpha%7D%5Csum_i%20p_%7Bti%7D%5Ctheta_i%5E%5Calpha)

      where angles are measures wrt to the jet axis after re-clustering witn anti-kt and the winner-takes-all recombination scheme

 - EnergyCorrelationFunction(alpha, R, optional_groomer)

    - this first grooms the jet if the groomer is set

    - the shape is computed as
    
      ![\Large \frac{1}{(\sum_i p_{ti})^2 R^\alpha}\sum_{i<j} p_{ti}p_{tj}\theta_{ij}^\alpha](https://latex.codecogs.com/svg.latex?%5CLarge%20%5Cfrac%7B1%7D%7B%28%5Csum_i%20p_%7Bti%7D%29%5E2%20R%5E%5Calpha%7D%5Csum_%7Bi%3Cj%7D%20p_%7Bti%7Dp_%7Btj%7D%5Ctheta_%7Bij%7D%5E%5Calpha)

 - LesHouchesMultiplicity(ktcut, R, optional_zcut, optional_beta)

    - This is similar to an IteratedSoftDrop multiplicity with a kt
      cut and an optional SoftDrop condition

    - In practice, ths counts the primary Lund declusterings which
      have a kt above the ktcut and, when present, pass the SoftDrop
      condition for the specified zcut and beta

# Notes

The other files are just quick tests obtained with Pythia
