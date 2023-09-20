# watch_it_public
Scripts to perform analysis of WATCH-IT study survey data and clinical factor modeling

# notes
- 'comorbidities.R' generates clinical factors/diseases with respect to WATCH-IT survey dates, in a manner analagous to the JEDI framework (see https://github.com/broadinstitute/jedi-public)
- 'combine_comorbidities.R' merges the comorbidities with the survey data
- 'main.R' performs the main analyses (multivariable logistic regression models, visualizations)
- '/functions/long_file_generator.R' is a helper function called for comorbidity generation

# dependencies
sktools (devtools::install_github:git@github.com:shaankhurshid/sktools.git)
