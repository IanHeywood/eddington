# eddington_bias.py

Uses the output of the [Wilman et al. (2008)](http://s-cubed.physics.ox.ac.uk/s3_sex) extragalactic radio simulation to estimate correction factors for Eddington bias in Euclidean-normalized differential source counts from radio surveys.

A FITS image of the survey sensitivity (e.g. as produced by [PyBDSF](http://www.astron.nl/citt/pybdsm/)) must be provided.

A suitable excerpt from the simulation database for this script can be downloaded from [here](http://www-astro.physics.ox.ac.uk/~ianh/e207dfe12f12c123f003651be8c224a626e9907d.result) (209 MB), which contains total intensity measurements at four frequencies for sources with 1400 MHz values greater than 0.1 uJy within a 100 square degree region.

This approach is not as good as the 'proper' way of doing this, namely injecting artificial sources into the survey map and recovering them using the same source finding software as was used to generate the catalogue. However it is cheaper, and probably good enough in observational regimes where resolution bias is not a major issue.

```
Usage: eddington_bias.py [options] rms_map

Options:
  -h, --help         show this help message and exit
  --s_min=S_MIN      Lower flux density limit for count determination [Jy]
                     (default = 20e-6)
  --s_max=S_MAX      Upper flux density limit for count determination [Jy]
                     (default = 0.1)
  --n_bins=N_BINS    Number of logarithmically-spaced bins between s_min and
                     s_max (default = 25)
  --sigma=SIGMA      Sigma cut above which a source is considered to be
                     detected (default = 5)
  --n_iter=N_ITER    Number of iterations to simulate noisy detections
                     (default = 1)
  --s_cubed=S_CUBED  S-cubed query .result file (default =
                     e207dfe12f12c123f003651be8c224a626e9907d.result)
  --area=AREA        Area of the S-cubed query [square degrees] (default =
                     100)
  --col=COL          Flux density column to use from S-cubed query (default =
                     itot_1400)
  --doplot           Plot the true and simulated source counts (default =
                     True)
  --pngname=PNGNAME  Name of output PNG file
```
