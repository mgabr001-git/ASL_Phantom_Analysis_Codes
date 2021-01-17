#!/bin/bash
############################################################################
##  This script creates a list of 3D volume files for each PLD, 
##  then uses created list file to merge them into 4D volume using fslmerge tool.
############################################################################


pld=(100 300 500 700 1000 1300 1700 2200 2700 3200 3500)
ld=1500
pump=(200 300 400)
#pump=(300)
#FILE=perf_list_LD$ld.txt

# rm perf_list_*.txt Variance_list_*.txt Weights_list_*.txt STD_list_*.txt meanCBF_list_*.txt  meanM0_list_*.txt DM_map_list_*.txt  Weighted_perf_list_*.txt  Mean_perf_smoothed_list_*.txt  



#if [ ! -f $FILE ]; then
  for ((iip=0; iip<${#pump[@]}; ++iip));
  do
   
     for ((idx=0; idx<${#pld[@]}; ++idx));
     do
        echo ../../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD${ld}_PLD${pld[idx]}_BS0_PUMP${pump[iip]}_ASL*/meanPerf_SUB01.nii  >>  perf_list_LD$ld.txt
	echo ../../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD${ld}_PLD${pld[idx]}_BS0_PUMP${pump[iip]}_ASL*/Variance_Perf_SUB01.nii  >>  Variance_list_LD$ld.txt
	echo ../../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD${ld}_PLD${pld[idx]}_BS0_PUMP${pump[iip]}_ASL*/Weights_SUB01.nii  >>  Weights_list_LD$ld.txt
	echo ../../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD${ld}_PLD${pld[idx]}_BS0_PUMP${pump[iip]}_ASL*/STD_Perf_SUB01.nii  >>  STD_list_LD$ld.txt
	echo ../../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD${ld}_PLD${pld[idx]}_BS0_PUMP${pump[iip]}_ASL*/meanCBF_SUB01.nii  >>  meanCBF_list_LD$ld.txt
	#	echo ../../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD${ld}_PLD${pld[idx]}_BS0_PUMP${pump[iip]}_ASL*/meanM0_SUB01.nii  >>  meanM0_list_LD$ld.txt
	#	echo ../../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD${ld}_PLD${pld[idx]}_BS0_PUMP${pump[iip]}_ASL*/Variance_M0_SUB01.nii  >>  Variance_M0_list_LD$ld.txt
	#	echo ../../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD${ld}_PLD${pld[idx]}_BS0_PUMP${pump[iip]}_ASL*/STD_M0_SUB01.nii  >>  STD_M0_list_LD$ld.txt
	
		
	#      echo $PWD/results/GK_PUMP300_BS0/Sigma3.0/DM_map_LD${ld}_PLD${pld[idx]}_sigma2.nii  >>  DM_map_list_LD$ld.txt
	echo $PWD/results/Smoothed_100M0_PUMP${pump[iip]}_BS0/Weighted_perf_LD${ld}_PLD${pld[idx]}_sigma3.nii  >>  Weighted_perf_list_LD$ld.txt
	echo $PWD/results/Smoothed_100M0_PUMP${pump[iip]}_BS0/Mean_perf_smoothed_LD${ld}_PLD${pld[idx]}_sigma3.nii  >>  Mean_perf_smoothed_list_LD$ld.txt
     done

     ## For pseudo M0, there is not PLD series,  100 M0 is created by merging unlabeled L/C images
     
     cp ../../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD${ld}_PLD1000_BS0_PUMP${pump[iip]}_ASL*/meanM0_SUB01.nii      $PWD/results/Smoothed_100M0_PUMP${pump[iip]}_BS0/
     cp ../../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD${ld}_PLD1000_BS0_PUMP${pump[iip]}_ASL*/Variance_M0_SUB01.nii $PWD/results/Smoothed_100M0_PUMP${pump[iip]}_BS0/
     cp ../../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD${ld}_PLD1000_BS0_PUMP${pump[iip]}_ASL*/STD_M0_SUB01.nii      $PWD/results/Smoothed_100M0_PUMP${pump[iip]}_BS0/

   
     echo "New file lists are created. Now merging files."
     fslmerge -tr   $PWD/results/Smoothed_100M0_PUMP${pump[iip]}_BS0/MeanPerf_Labeled_PLDseries_LD$ld.nii `cat  perf_list_LD$ld.txt`      6
     fslmerge -tr   $PWD/results/Smoothed_100M0_PUMP${pump[iip]}_BS0/Variance_Labeled_PLDseries_LD$ld.nii `cat  Variance_list_LD$ld.txt`  6
     fslmerge -tr   $PWD/results/Smoothed_100M0_PUMP${pump[iip]}_BS0/Weights_PLDseries_LD$ld.nii          `cat  Weights_list_LD$ld.txt`   6
     fslmerge -tr   $PWD/results/Smoothed_100M0_PUMP${pump[iip]}_BS0/STD_Labeled_PLDseries_LD$ld.nii      `cat  STD_list_LD$ld.txt`       6
     fslmerge -tr   $PWD/results/Smoothed_100M0_PUMP${pump[iip]}_BS0/MeanCBF_Labeled_PLDseries_LD$ld.nii  `cat  meanCBF_list_LD$ld.txt`   6
  #   fslmerge -tr   $PWD/results/Smoothed_100M0_PUMP${pump[iip]}_BS0/M0_Uabeled_LD$ld.nii `cat  meanM0_list_LD$ld.txt`    6

   #  fslmerge -tr   $PWD/results/Smoothed_100M0_PUMP${pump[iip]}_BS0/Variance_100M0_Unlabeled_LD$ld.nii   `cat  Variance_M0_list_LD$ld.txt`    6
  #   fslmerge -tr   $PWD/results/Smoothed_100M0_PUMP${pump[iip]}_BS0/STD_100M0_LD$ld.nii                  `cat  STD_M0_list_LD$ld.txt`         6
     
     # fslmerge -tr   $PWD/results/GK_PUMP300_BS0/Sigma3.0/DM_map_LD$ld.nii                     `cat  DM_map_list_LD$ld.txt`    6
     fslmerge -tr   $PWD/results/Smoothed_100M0_PUMP${pump[iip]}_BS0/Weighted_perf_LD$ld.nii              `cat  Weighted_perf_list_LD$ld.txt`       6
     fslmerge -tr   $PWD/results/Smoothed_100M0_PUMP${pump[iip]}_BS0/Mean_perf_smoothed_LD$ld.nii         `cat  Mean_perf_smoothed_list_LD$ld.txt`  6
     


# rm perf_list_*.txt Variance_list_*.txt Weights_list_*.txt STD_list_*.txt meanCBF_list_*.txt  meanM0_list_*.txt DM_map_list_*.txt  Weighted_perf_list_*.txt  Mean_perf_smoothed_list_*.txt  
rm perf_list_*.txt Variance_list_*.txt Weights_list_*.txt STD_list_*.txt meanCBF_list_*.txt  meanM0_list_*.txt  Weighted_perf_list_*.txt  Mean_perf_smoothed_list_*.txt  Variance_M0_list_*.txt  STD_M0_list_*.txt
     
 done
   
   

#else
    
#   echo "File $FILE already exists! (Re-)Creating MeanPerf_LD$ld.nii "
#   fslmerge -tr  MeanPerf_Labeled_Pump0_PLDseries_LD$ld.nii `cat  perf_list_LD$ld.txt`      6
#   fslmerge -tr  Variance_Labeled_Pump0_PLDseries_LD$ld.nii `cat  Variance_list_LD$ld.txt`  6
#   fslmerge -tr  STD_Labeled_Pump0_PLDseries_LD$ld.nii      `cat  STD_list_LD$ld.txt`       6
#   fslmerge -tr  MeanCBF_Labeled_Pump0_PLDseries_LD$ld.nii  `cat  meanCBF_list_LD$ld.txt`   6
#   fslmerge -tr  MeanM0_Labeled_Pump0_PLDseries_LD$ld.nii   `cat  meanM0_list_LD$ld.txt`    6
   
   
   
#fi







