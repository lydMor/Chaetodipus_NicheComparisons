# Chaetodipus_NicheComparisons
Scripts and data used to run pairwise niche comparisons for all focal species in a pocket mouse genus
Here I outline the steps I took to conduct pairwise niche comparisons for an upcoming manuscript. I've provided the data so it's easy to follow along. 

This project exists in two steps: 
1) clean your data and extract important environmental variables for comparison
   NOTE: Most of the record cleaning/verification for this project was done by field experts and many occurrence points were collected in the field over years of study. For projects that use mostly or entirely publicly available data (e.g. GBIF), I would suggest more thorough cleaning and processing steps than what was done here.
   
3) Use those variables to run pairwise niche similarity tests for all members of your group.
   NOTE: This process can optionally include a variable selection step. I didn't use it here, because I wanted all comparisons across all species pairs to be consistent and translatable. If you're just interested in one pairwise comparison, or if consistency isn't a concern, check out Jason Brown's Humboldt page for instructions on how to automate variable selection as part of this step:[ https://github.com/jasonleebrown/humboldt/tree/master]( https://github.com/jasonleebrown/humboldt/tree/master)
