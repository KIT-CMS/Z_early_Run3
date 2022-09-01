#!/bin/sh

   # enter which channel you want (mm/ee/mmet/emet)
   channel=$1

   nthreads=16
   if [ ! -z "$2" ]
   then
      nthreads=$2
   fi

   # output postfix
   postfix="DYWNLO"
   # postfix="DYWLO"

   # if 0 then it sums the runs, if 1 you get a run by run plot
   runPlot=1

   # input the max luminosity you want runs to be grouped by (if 0 it will not group runs)
   groupRuns=50

   # if 1 it will seperate variables by things you can define in control_binning.py (useful for QCD template)
   seperateVariables=0

   # luminosity of the summed runs in fb-1 (if you use runPlot it will not matter what number it here because it will be weighted to 1 pb-1)
   totalLuminosity=0.656885772534

   # Run numbers and correspoding luminosities you want to plot
   run_list="355872,355892,355912,355913,355921,355933,355942,355988,355989,356043,356071,356074,356075,356076,356077,356135,356309,356316,356321,356322,356323,356371,356375,356378,356381,356383,356385,356386,356426,356428,356433,356434,356435,356446"
   lumi_list="21.029407135,1.976640736,9.007156928,5.725368440,22.873241336,0.724271840,1.371972162,2.780268121,1.265671750,5.567402786,15.248489527,2.999361674,13.902884694,17.133845292,49.787188261,2.612524012,10.957679862,9.292921771,7.684875420,2.060615508,67.033187162,1.569606629,6.515209807,41.464490983,131.203129037,2.867059944,2.467573112,8.706233768,2.563759568,45.822680251,47.969426483,1.724343482,0.441028997,92.536256058"

   # variables by channel (this can be changed to fit what variables you want plotted)
   if [ $channel == 'mm' ]
   then
      variable_list='pt_1,eta_1,phi_1,m_vis,met_uncorrected,metphi_uncorrected'
      # variable_list='pt_1,pt_2,eta_1,eta_2,phi_1,phi_2,m_vis,pt_vis,mt_uncorrected,met_uncorrected,metphi_uncorrected,pfmet_uncorrected,pfmetphi_uncorrected,nTrackerLayers_1,nTrackerLayers_2,nStations_1,nStations_2'

   elif [ $channel == 'mmet' ]
   then
      variable_list='pt_1,eta_1,phi_1,met_uncorrected,metphi_uncorrected'
      # variable_list='pt_1,eta_1,phi_1,mt_uncorrected,met_uncorrected,metphi_uncorrected,pfmet_uncorrected,pfmetphi_uncorrected,nTrackerLayers_1,nStations_1'

   elif [ $channel == 'ee' ]
   then
      variable_list='pt_1,eta_1,phi_1,m_vis,met_uncorrected,metphi_uncorrected'
      # variable_list='pt_1,pt_2,eta_1,eta_2,etaSC_1,etaSC_2,phi_1,phi_2,m_vis,pt_vis,mt_uncorrected,met_uncorrected,metphi_uncorrected,pfmet_uncorrected,pfmetphi_uncorrected,deltaetaSC_1,deltaetaSC_2,eInvMinusPInv_1,eInvMinusPInv_2,hoe_1,hoe_2,scEtOverPt_1,scEtOverPt_2,sieie_1,sieie_2,lostHits_1,lostHits_2'
      
   elif [ $channel == 'emet' ]
   then
      variable_list='pt_1,eta_1,phi_1,met_uncorrected,metphi_uncorrected'
      # variable_list='pt_1,eta_1,etaSC_1,phi_1,mt_uncorrected,met_uncorrected,metphi_uncorrected,pfmet_uncorrected,pfmetphi_uncorrected,deltaetaSC_1,eInvMinusPInv_1,hoe_1,scEtOverPt_1,sieie_1,lostHits_1'
   fi

   python shapes/produce_shapes.py --channels $channel --era 2022 \
   --output-file output/earlyRun3_2022_"$channel"_runPlot"$runPlot"_"$postfix" \
   --directory /ceph/moh/CROWN_samples/Run3V01/ntuples \
   --$channel-friend-directory /ceph/moh/CROWN_samples/Run3V01/friends/crosssection \
   --num-processes 8 --num-threads $nthreads \
   --optimization-level 1 --control-plots \
   --run-plot $runPlot --group-runs $groupRuns \
   --run-list $run_list --lumi-list $lumi_list \
   --seperate-variables $seperateVariables \
   --control-plot-set $variable_list --total-lumi $totalLuminosity \
   --ntuple_type crown --skip-systematic-variations
