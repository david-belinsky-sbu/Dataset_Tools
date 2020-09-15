# Fully_Annotated_Dataset_Creation
Generates a fully annotated dataset from selected region annotation, segmentation mask annotation, and class dot annotations.

Directory Structure

|-Fully_Annotated_Dataset_Creation/  
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-feature_extraction_from_manual_annotation.py  
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-misc/  
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-metrics/  
|-TCGA_Data/  
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-manifest.csv  
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-*.json  
|-Segmentation_masks/  
  

Before running modify the first few variables in feature_extraction_from_manual_annotation.py to point them at the manifest, json files, segmentation masks, and raw wsi data.
