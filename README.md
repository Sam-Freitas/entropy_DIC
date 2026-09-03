# entropy_DIC
Using a 2D entropy filter to segment DIC (Differential Interference Contrast) imaging of cells for automated sctrach assay analysis

Our implementation is for scratch assays (ref: https://www.nature.com/articles/nprot.2007.30)

------------------------------
A gif example of the systems output (reload page if gif froze)

<img src="https://github.com/Sam-Freitas/entropy_DIC/blob/master/scripts/output/Example_output/RENCA_A28_10FBS_3HAA_1.gif" alt="GIF showing an example export" width="885" height="704" loop=infinite>


### setup: File organization

* this expects the experimental folders to be in a specific organizational structure
    * Experiment folder (containing_folder)
        * Individual conditions with repeats
            * images sorted by timepoint (TIFF)

<img src="https://github.com/Sam-Freitas/entropy_DIC/blob/master/file_setup.png" alt="PNG showing the file setup" width="2040" height="1122">


### usage: MATLAB

* Read and change the settings if necessary in **/scripts/Entropy_segment_cells_Batch.m**
* most important changes are
    * **containing_folder** (path to where the TIFF folders are)
    * **experiment_name** (what the export should be called)
* Run
```
Entropy_segment_cells_Batch.m
```
