# eddington

Uses the output of the [Wilman et al. (2008)](http://s-cubed.physics.ox.ac.uk/s3_sex) extragalactic radio simulation to estimate correction factors for Eddington bias in Euclidean-normalized differential source counts from radio surveys.

A FITS image of the survey sensitivity (e.g. as produced by [http://www.astron.nl/citt/pybdsm/](PyBDSF)) must be provided.

A suitable excerpt from the simulation database for this script can be downloaded from [here](http://www-astro.physics.ox.ac.uk/~ianh/e207dfe12f12c123f003651be8c224a626e9907d.result) (209 MB), which contains total intensity measurements at four frequencies for sources with 1400 MHz values greater than 0.1 uJy within a 100 square degree region.

This approach is not as good as the 'proper' way of doing this, namely injecting artificial sources into the survey map and recovering them using the same source finding software as was used to generate the catalogue. However it is cheaper, and probably good enough in observational regimes where resolution bias is not a major issue.

