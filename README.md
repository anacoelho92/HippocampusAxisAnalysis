# Hippocampus Axis Analysis
Extract centerline of hippocampus to perform analysis (volumetric, tractography, etc.) along the line 

### Step by step instructions

1. Run extract_hippo.sh to extract hippocampus from Freesurfer output
   - usage: bash extract_hippo.sh *dir* *subj_file*
   - dir: directory with output of Freesurfer of all subjects
   - subj_file: file with subjects IDs
   - output will be saved in mri folder of each subject
   
2. Run get_centerline_all.m in Matlab to extract first version of the centerline (script extracted from [here](https://github.com/garikoitz/hippovol):
   - need to change directory and hemi (lh - left; rh - right)
   - directory is where the output of Freesurfer of all subjects is stored (with hippocampus already extracted)
   - this will create 3 files: stl (3D rendering of hippocampus); bezier_pts.txt (coordinates of the bezier line); ctrl_pts.txt (coordinates of the control points used to define the bezier line)
 
3. Run centerline_analysis.py (centerline_analysis_rodent.py for rodent data) in python to improve the centerline, create the orthogonal planes and compute volumes of the sections
    - vtk library needs to be installed (pip install vtk)
    - usage: python centerline_analysis.py *working_dir* *subjects_file* *roi_name* *LoR* *nplanes*
    - working_dir: directory with hippocampus 3D stl file, centerline points and control points
    - subjects_file: file with subjects codes
    - roi_name: Hippo (same name in other files)
    - LoR: lh - left; rh - right
    - nplanes: number of planes to define (usually 20)
    - this will create 3 txt files: volumes (volume of each section for all subjects); areas (area of each section for all subjects); totalVolumes (total volume of the hippocampus for all subjects) 
  
  
