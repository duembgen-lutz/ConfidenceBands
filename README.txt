This project is about nonparametric confidence bands for a (continuous) distribution function.

The R file "GoF_CritValues.R" contains programs to compute the distribution function and exact critical values of four types of goodness-of-fit tests:
-  the Kolmogorov-Smirnov test,
-  a weighted Kolmogorov-Smirnov test proposed by
   Stepanova & Pavlenko (2018, Journal of Theoretical Probability),
-  the generalized Berk-Jones tests introduced by
   Jager & Wellner (2007, Annals of Statistics)
-  the goodness-of-fit tests introduced by
   Duembgen & Wellner (2022, arxiv:1402.2918)
The R file "GoF_CritValues_Demo.R" contains examples on how to use the functions in "GoF_CritValues.R" as well as some critical values.

The computation of the distribution function is based on a version of Noe's (1972) recursion, see NoesRecursion.pdf. It is implemented in the R file "NoeRecursion.R" and illustrated in "NoeRecursion_Demo.R". For large sample sizes, function is rather slow, but the numerical accuracy seems to be very good for sample sizes up to at least 8000. This is achieved by working with log-probabilities instead of probabilities in a suitable way.