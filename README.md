# MS Annika 3.0 Non-Cleavable Crosslink Search Results

This repository contains scripts used for the analysis and visualization of the
results of the [MS Annika 3.0](https://github.com/hgb-bin-proteomics/MSAnnika) non-cleavable crosslink search.

The structure of this repository is as follows:
- `CElegansFasta`: Creation of the filtered *C. elegans* protein database
- `CElegansFastaDecoy`: Creation of the filtered *C. elegans* protein database which also accounts for decoy PSMs
- `CElegansResults`: Analysis of the results for the *C. elegans* crosslink search using the protein database from `CElegansFasta`
- `CElegansResultsDecoy`: Analysis of the results for the *C. elegans* crosslink search using the protein database from `CElegansFastaDecoy`
- `CElegansPPIs`: Protein-protein-interaction analysis of the *C. elegans* results including visualization with xiView
- `Lenz`: Analysis of the results from the crosslink search using the dataset by Lenz and co-workers ([JPST000845](https://repository.jpostdb.org/entry/JPST000845))
- `Peplib_Beveridge`: Analysis of the results form the crosslink search using the dataset by Beveridge and co-workers ([PXD014337](https://www.ebi.ac.uk/pride/archive/projects/PXD014337))
- `Peplib_Matzinger`: Analysis of the results form the crosslink search using the dataset by Matzinger and co-workers ([PXD029252](https://www.ebi.ac.uk/pride/archive/projects/PXD029252))
- `Xi_Annika_scripts`: Miscellaneous scripts to aid in the processing of xi output files

For details please refer to the [MS Annika 3.0 publication](https://doi.org/10.1038/s42004-024-01386-x)!

## Contact

- [micha.birklbauer@fh-hagenberg.at](mailto:micha.birklbauer@fh-hagenberg.at)
