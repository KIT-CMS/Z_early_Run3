#!/bin/sh

   #enter which channel you want (mm/ee/mmet/emet)
   channel='mmet'

   #if 0 then it sums the runs, if 1 you get a run by run plot
   runPlot=0

   #input the max luminosity you want runs to be grouped by (if 0 it will not group runs)
   groupRuns=0

   #if 1 it will seperate variables by things you can define in control_binning.py (useful for QCD template)
   seperateVariable=1

   #luminosity of the summed runs in fb-1 (if you use runPlot it will not matter what number it here because it will be weighted to 1 pb-1)
   totalLuminosity=59.832045316

   #Run numbers and correspoding luminosities you want to plot
   run_list="355443,355444,355445,355558,355559,355679,355680,355768,355769"
   lumi_list="1.971790835,1.243983634,1.537414347,11.7548517,6.878012741,1.078718626,46.78980232,1.222612315,16.16789967"

   #variables by channel (this can be changed to fit what variables you want plotted)
   if [ $channel == 'mm' ]
   then
      variable_list='met_uncorrected_barrel_iso5'
      #variable_list="pt_1,pt_2,eta_1,eta_2,phi_1,phi_2,m_vis,pt_vis,nTrackerLayers_1_barrel,nTrackerLayers_1_endcap,nTrackerLayers_2_barrel,nTrackerLayers_2_endcap,nStations_1_barrel,nStations_2_barrel,nStations_1_endcap,nStations_2_endcap,met_uncorrected,metphi_uncorrected,pfmet_uncorrected,pfmetphi_uncorrected"

   elif [ $channel == 'mmet' ]
   then
      #variable_list="ptOverMt,mt_uncorrected,met_uncorrected"
      #variable_list="mt_uncorrected,met_uncorrected,ptOverMt"
      variable_list="mt_uncorrected,met_uncorrected,ptOverMt"

   elif [ $channel == 'ee' ]
   then
      variable_list="mt_uncorrected,met_uncorrected,ptOverMt"
      
   elif [ $channel == 'emet' ]
   then
      variable_list="mt_uncorrected,met_uncorrected,ptOverMt"
   fi

   python shapes/produce_shapes.py --channels $channel --output-file output/earlyRun3_crown_2018_"$channel"\
   --directory /ceph/moh/CROWN_samples/EarlyRun3_V12/ntuples\
   --$channel-friend-directory /ceph/moh/CROWN_samples/EarlyRun3_V12/friends/crosssection\
   --era 2018 --num-processes 64 --num-threads 64 --optimization-level 1 --control-plots --run-plot $runPlot --group-runs $groupRuns --seperate-variables $seperateVariable\
   --control-plot-set $variable_list --total-lumi $totalLuminosity --run-list $run_list --lumi-list $lumi_list --ntuple_type crown --skip-systematic-variations

   