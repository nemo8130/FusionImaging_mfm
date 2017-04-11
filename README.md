# FusionImaging_mfm
Requires: `MATLAB (version R2016a or higher)`  

### (c) Indian Copyright (Diary Number. 9692/2016-CO/SW)

## Funding: The work was supported by the Department of Science and Technology – Science and Engineering Research Board (DST-SERB research grant PDF/2015/001079) and TEQIP, (R&D section) Centre of Excellence for Systems Biology and Nano Technology, University of Calcutta


Installation
```
git clone https://github.com/nemo8130/FusionImaging_mfm
cd FusionImaging_mfm
run the script FusionImaging_mfm.m in MATLAB command window as a function
```

```
=======
% FusionImaging_mfm is a MATLAB code for processing and analyses of Images 
% generated from Magnetic Force Microscopy (MFM) 
% The code was compiled in MATLAB Version: 9.0.0.341360 (R2016a) 
% Users are recommended to use this or higher version with the 
% 'Image Processing Toolbox' installed.
%
% The software converges a collection of raw MFM images for a given sample 
% collected at different Lift-heights into Fusion Images. 
% 
% The function needs two and exactly two input arguments
```

**Usage : FusionImaging_mfm (sampfold, datmm)**

```
% where 
% sampfold: Is the sample folder with full path containing raw MFM .tif images\n')
% and 
% datmm: Is the data file containing Minimum and Maximum voltages (in units of mV) 
% of each image in the 'sampfold' stored in the tabular format specified below. 
%
% Directory            Filename                     Min_Voltage(mV)         Max_Voltage(mV)
% APTES_1_P             ap_010.tif                   -1200.00                -1100.00
% APTES_1_P             ap_020.tif                   -1100.00                  979.80
% APTES_1_P             ap_030.tif                    -970.20                  875.00
% APTES_1_P             ap_040.tif                    -821.60                  763.50
% APTES_1_P             ap_050.tif                    -821.60                  763.50
% APTES_1_P             ap_060.tif                    -730.40                  681.00
% APTES_1_P             ap_070.tif                    -824.30                  750.50
% BSA_1_P               bs_060.tif                      -4.90                    4.90
% BSA_1_P               bs_080.tif                      -4.40                    4.50
% BSA_1_P               bs_090.tif                      -4.70                    4.60
% CITRATE_1_P           c1_010.tif                   -3200.00                 2500.00
% CITRATE_1_P           c1_020.tif                   -4600.00                 4000.00
% CITRATE_1_P           c1_030.tif                   -1900.00                 1500.00
% CITRATE_1_P           c1_040.tif                   -3400.00                 2600.00
% CITRATE_1_P           c1_050.tif                   -3300.00                 2400.00
% CITRATE_1_P           c1_060.tif                   -3700.00                 2600.00
% CITRATE_2_P           c2_010.tif                    -174.70                  135.30
% CITRATE_2_P           c2_020.tif                    -163.90                  129.70
% CITRATE_2_P           c2_030.tif                    -169.70                  126.00
%
% An example datmm file (Min_Max_Table_example_phase_data.format) is provided with this distribution.
% Note: The datmm file should contain the commented header starting with a '%'
% Images must have '.tif' extension
%
% Specification of Full path needs to be consistent with the OS (Unix / Windows: see example 'sampfold')
% 
% Full path may or may not end with a '/' (Linux) or '\' (Windows) : The program can handle both
% the datmm file can be directly called from the current directory by the file name alone
%
% logz=FusionImaging_mfm (sampfold, datmm) will return the (NxN) matrix containing 
% coefficients of the inverse-square term of the fitted polynomial; 
% (magnitudes scaled to natural logarithm, log_e)
% These coefficients are representative of the magnetic fields 
% corresponding to each small square grid in the 'two dimensional' image-space. 
%
% Example Input:  
%
% sampfold='/home/sankar/Dropbox/sankar-puja/Data/Data_Widout_VW/phase_data/' (Linux / Unix) 
% sampfold='C:\user\Sankar\Desktop\mfmdata\Data_Widout_VW\phase_data\'; (Windows / PC)
%
% datmm='/home/sankar/Dropbox/sankar-puja/Data/SOFTWARE/Min_Max_Table_example_phase_data.format'; (File with full path: Linux)
% datmm='Min_Max_Table_example_phase_data.format';  (File kept in the Current Directory)
%
% datmm='C:\user\Sankar\Data\SOFTWARE\Min_Max_Table_example_phase_data.format'; (File with full path: Windows)
%
%
```
