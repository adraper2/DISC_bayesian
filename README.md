## Notre Dame REU: Project 2

This is a repository for a Bayesian Poisson point process model we are recreating from a paper written by Malcolm Itter. Malcolm Itter is a PhD candidate at Michigan State University studying forest demographic processes. Read more about his work <a href= "https://www.mitter-forestecology.com/about.html">here</a>.

The Bayesian model attempts to predict local fires surrounding lakes by investigating charcoal counts from the past as an indicator of a regional or local fire. The Poisson process is a fitting addition to the model given the use of charcoal counts. With that being said, the original model was developed on the notion that charcoal data could be separated into a background and a foreground where the background represents regional fire noise and the foreground represents peak charcoal amounts that are caused by a local fire. This idea held true in Alaskan lakes, but now we are interested as to whether or not this is the case with lakes residing in the Midwest. We have 18 lakes of interest that we will model separately first and then, jointly.

Read through the walkthrough of our code and processes in the <a href="https://github.com/adraper2/DISC_bayesian/blob/master/walkthrough.Rmd">walkthrough.Rmd</a> file.
