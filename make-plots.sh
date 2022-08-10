#!/bin/sh

#channel options include (mm/ee/mmet/emet)
channel='mm'

#if this value is not 0 it will normalize the data to the MC
matchData=0

#if you summed the runs change this to match the total lumi for the label
lumiLabel="0.088645086184"

#if value is not 0 this will write all the plots you make to a latex file and create a slideshow with them
#Always adds on so you must go to (presentation/control-plots-slides.tex) if you want to delete/change the presentation
latex=0

#sets up a latex environment
export LATEX=/ceph/mlink/texlive/2022/
export PATH=$LATEX/bin/x86_64-linux:$PATH
export MANPATH=$LATEX/texmf-dist/doc/man:$MANPATH
export INFOPATH=$LATEX/texmf-dist/doc/info:$INFOPATH

#variables by channel
   if [ $channel == 'mm' ]
   then
      variable_list="pt_1,pt_2,eta_1,eta_2,phi_1,phi_2,m_vis,pt_vis,nTrackerLayers_1_endcap,nTrackerLayers_2_endcap,nStations_1_endcap,nStations_2_endcap,nTrackerLayers_1_barrel,nTrackerLayers_2_barrel,nStations_1_barrel,nStations_2_barrel,met_uncorrected,metphi_uncorrected,pfmet_uncorrected,pfmetphi_uncorrected"

   elif [ $channel == 'mmet' ]
   then
      variable_list="pt_1_pos,pt_1_neg,eta_1_pos,eta_1_neg,phi_1_pos,phi_1_neg,mt_uncorrected_pos,mt_uncorrected_neg,nTrackerLayers_1_barrel_pos,nTrackerLayers_1_barrel_neg,nTrackerLayers_1_endcap_pos,nTrackerLayers_1_endcap_neg,nStations_1_barrel_pos,nStations_1_barrel_neg,nStations_1_endcap_pos,nStations_1_endcap_neg,met_uncorrected_pos,met_uncorrected_neg,metphi_uncorrected_pos,metphi_uncorrected_neg,pfmet_uncorrected_pos,pfmet_uncorrected_neg,pfmetphi_uncorrected_pos,pfmetphi_uncorrected_neg"

   elif [ $channel == 'ee' ]
   then
      variable_list="pt_1,pt_2,eta_1,eta_2,phi_1,phi_2,m_vis,pt_vis,deltaetaSC_1_barrel_nv015,deltaetaSC_2_barrel_nv015,deltaetaSC_1_endcap_nv015,deltaetaSC_2_endcap_nv015,eInvMinusPInv_1_barrel_nv015,eInvMinusPInv_2_barrel_nv015,eInvMinusPInv_1_endcap_nv015,eInvMinusPInv_2_endcap_nv015,hoe_1_barrel_nv015,hoe_2_barrel_nv015,hoe_1_endcap_nv015,hoe_2_endcap_nv015,scEtOverPt_1_barrel_nv015,scEtOverPt_2_barrel_nv015,scEtOverPt_1_endcap_nv015,scEtOverPt_2_endcap_nv015,sieie_1_barrel_nv015,sieie_2_barrel_nv015,sieie_1_endcap_nv015,sieie_2_endcap_nv015,lostHits_1_barrel_nv015,lostHits_2_barrel_nv015,lostHits_1_endcap_nv015,lostHits_2_endcap_nv015,deltaetaSC_1_barrel_nv1530,deltaetaSC_2_barrel_nv1530,deltaetaSC_1_endcap_nv1530,deltaetaSC_2_endcap_nv1530,eInvMinusPInv_1_barrel_nv1530,eInvMinusPInv_2_barrel_nv1530,eInvMinusPInv_1_endcap_nv1530,eInvMinusPInv_2_endcap_nv1530,hoe_1_barrel_nv1530,hoe_2_barrel_nv1530,hoe_1_endcap_nv1530,hoe_2_endcap_nv1530,scEtOverPt_1_barrel_nv1530,scEtOverPt_2_barrel_nv1530,scEtOverPt_1_endcap_nv1530,scEtOverPt_2_endcap_nv1530,sieie_1_barrel_nv1530,sieie_2_barrel_nv1530,sieie_1_endcap_nv1530,sieie_2_endcap_nv1530,lostHits_1_barrel_nv1530,lostHits_2_barrel_nv1530,lostHits_1_endcap_nv1530,lostHits_2_endcap_nv1530,deltaetaSC_1_barrel_nv3045,deltaetaSC_2_barrel_nv3045,deltaetaSC_1_endcap_nv3045,deltaetaSC_2_endcap_nv3045,eInvMinusPInv_1_barrel_nv3045,eInvMinusPInv_2_barrel_nv3045,eInvMinusPInv_1_endcap_nv3045,eInvMinusPInv_2_endcap_nv3045,hoe_1_barrel_nv3045,hoe_2_barrel_nv3045,hoe_1_endcap_nv3045,hoe_2_endcap_nv3045,scEtOverPt_1_barrel_nv3045,scEtOverPt_2_barrel_nv3045,scEtOverPt_1_endcap_nv3045,scEtOverPt_2_endcap_nv3045,sieie_1_barrel_nv3045,sieie_2_barrel_nv3045,sieie_1_endcap_nv3045,sieie_2_endcap_nv3045,lostHits_1_barrel_nv3045,lostHits_2_barrel_nv3045,lostHits_1_endcap_nv3045,lostHits_2_endcap_nv3045,met_uncorrected,metphi_uncorrected,pfmet_uncorrected,pfmetphi_uncorrected"
   elif [ $channel == 'emet' ]
   then
      variable_list="pt_1_pos,pt_1_neg,eta_1_pos,eta_1_neg,phi_1_pos,phi_1_neg,mt_uncorrected_pos,mt_uncorrected_neg,deltaetaSC_1_barrel_pos_nv015,deltaetaSC_1_barrel_neg_nv015,eInvMinusPInv_1_barrel_pos_nv015,eInvMinusPInv_1_barrel_neg_nv015,hoe_1_barrel_pos_nv015,hoe_1_barrel_neg_nv015,scEtOverPt_1_barrel_pos_nv015,scEtOverPt_1_barrel_neg_nv015,sieie_1_barrel_pos_nv015,sieie_1_barrel_neg_nv015,lostHits_1_barrel_pos_nv015,lostHits_1_barrel_neg_nv015,deltaetaSC_1_barrel_pos_nv1530,deltaetaSC_1_barrel_neg_nv1530,eInvMinusPInv_1_barrel_pos_nv1530,eInvMinusPInv_1_barrel_neg_nv1530,hoe_1_barrel_pos_nv1530,hoe_1_barrel_neg_nv1530,scEtOverPt_1_barrel_pos_nv1530,scEtOverPt_1_barrel_neg_nv1530,sieie_1_barrel_pos_nv1530,sieie_1_barrel_neg_nv1530,lostHits_1_barrel_pos_nv1530,lostHits_1_barrel_neg_nv1530,deltaetaSC_1_barrel_pos_nv3045,deltaetaSC_1_barrel_neg_nv3045,eInvMinusPInv_1_barrel_pos_nv3045,eInvMinusPInv_1_barrel_neg_nv3045,hoe_1_barrel_pos_nv3045,hoe_1_barrel_neg_nv3045,scEtOverPt_1_barrel_pos_nv3045,scEtOverPt_1_barrel_neg_nv3045,sieie_1_barrel_pos_nv3045,sieie_1_barrel_neg_nv3045,lostHits_1_barrel_pos_nv3045,lostHits_1_barrel_neg_nv3045,deltaetaSC_1_endcap_pos_nv015,deltaetaSC_1_endcap_neg_nv015,eInvMinusPInv_1_endcap_pos_nv015,eInvMinusPInv_1_endcap_neg_nv015,hoe_1_endcap_pos_nv015,hoe_1_endcap_neg_nv015,scEtOverPt_1_endcap_pos_nv015,scEtOverPt_1_endcap_neg_nv015,sieie_1_endcap_pos_nv015,sieie_1_endcap_neg_nv015,lostHits_1_endcap_pos_nv015,lostHits_1_endcap_neg_nv015,deltaetaSC_1_endcap_pos_nv1530,deltaetaSC_1_endcap_neg_nv1530,eInvMinusPInv_1_endcap_pos_nv1530,eInvMinusPInv_1_endcap_neg_nv1530,hoe_1_endcap_pos_nv1530,hoe_1_endcap_neg_nv1530,scEtOverPt_1_endcap_pos_nv1530,scEtOverPt_1_endcap_neg_nv1530,sieie_1_endcap_pos_nv1530,sieie_1_endcap_neg_nv1530,lostHits_1_endcap_pos_nv1530,lostHits_1_endcap_neg_nv1530,deltaetaSC_1_endcap_pos_nv3045,deltaetaSC_1_endcap_neg_nv3045,eInvMinusPInv_1_endcap_pos_nv3045,eInvMinusPInv_1_endcap_neg_nv3045,hoe_1_endcap_pos_nv3045,hoe_1_endcap_neg_nv3045,scEtOverPt_1_endcap_pos_nv3045,scEtOverPt_1_endcap_neg_nv3045,sieie_1_endcap_pos_nv3045,sieie_1_endcap_neg_nv3045,lostHits_1_endcap_pos_nv3045,lostHits_1_endcap_neg_nv3045,met_uncorrected_pos,met_uncorrected_neg,metphi_uncorrected_pos,metphi_uncorrected_neg,pfmet_uncorrected_pos,pfmet_uncorrected_neg,pfmetphi_uncorrected_pos,pfmetphi_uncorrected_neg"
   fi

#creates the plots
source plotting/plot_shapes_control.sh 2018 output/earlyRun3_crown_2022_"$channel".root $variable_list None $matchData $lumiLabel $latex $channel earlyRun3$channel

#compiles control-plots-slides.tex into a pdf
if [ $latex != 0 ]
then
   pdflatex presentation/control-plots-slides.tex
fi