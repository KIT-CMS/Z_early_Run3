#! /bin/bash

OUT=Run3V01_outputs
IN=/ceph/moh/CROWN_samples/Run3V01
METVAR=met_uncorrected
METPHIVAR=metphi_uncorrected
LUMI=.657

####################################################################################################
#--------------------------------------------------------------------------------------------------#
#-            Z --> mm                                                                           -#
# -------------------------------------------------------------------------------------------------#
####################################################################################################

#Data
root -l -b -q fitRecoil_master.C\(3,3,1,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZmm_data_triple&
root -l -b -q fitRecoil_master.C\(0,0,1,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZmm_data_CB&

#Data and background
root -l -b -q fitRecoil_master.C\(3,3,0,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZmm_dataBckg_triple&
root -l -b -q fitRecoil_master.C\(0,0,0,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZmm_dataBckg_CB&

#Simulation
#root -l -b -q fitRecoil_master.C\(3,3,1,0,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZmm_sim_triple&
root -l -b -q fitRecoil_master.C\(0,0,1,0,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZmm_sim_CB&

####################################################################################################
#--------------------------------------------------------------------------------------------------#
#-            Z --> ee                                                                           -#
# -------------------------------------------------------------------------------------------------#
####################################################################################################

#Data
root -l -b -q fitRecoil_master.C\(3,3,1,1,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZee_data_triple&
root -l -b -q fitRecoil_master.C\(0,0,1,1,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZee_data_CB&

#Data and background
#root -l -b -q fitRecoil_master.C\(3,3,0,1,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZee_dataBckg_triple&
#root -l -b -q fitRecoil_master.C\(0,0,0,1,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZee_dataBckg_CB&

#Simulation
#root -l -b -q fitRecoil_master.C\(3,3,1,1,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZee_sim_triple&
root -l -b -q fitRecoil_master.C\(0,0,1,1,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZee_sim_CB&

####################################################################################################
#--------------------------------------------------------------------------------------------------#
#-            W+ --> mmet                                                                           -#
# -------------------------------------------------------------------------------------------------#
####################################################################################################

#Simulation
root -l -b -q fitRecoil_master.C\(3,3,1,0,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWpmmet_sim_triple&
root -l -b -q fitRecoil_master.C\(0,0,1,0,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWpmmet_sim_CB& 

####################################################################################################
#--------------------------------------------------------------------------------------------------#
#-            W- --> mmet                                                                           -#
# -------------------------------------------------------------------------------------------------#
####################################################################################################

#Simulation
root -l -b -q fitRecoil_master.C\(3,3,1,0,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWnmmet_sim_triple&
root -l -b -q fitRecoil_master.C\(0,0,1,0,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWnmmet_sim_CB&

####################################################################################################
#--------------------------------------------------------------------------------------------------#
#-            W- --> emet                                                                           -#
# -------------------------------------------------------------------------------------------------#
####################################################################################################

#Simulation
root -l -b -q fitRecoil_master.C\(3,3,1,1,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWnemet_sim_triple&
root -l -b -q fitRecoil_master.C\(0,0,1,1,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWnemet_sim_CB&

####################################################################################################
#--------------------------------------------------------------------------------------------------#
#-            W+ --> emet                                                                           -#
# -------------------------------------------------------------------------------------------------#
####################################################################################################

#Simulation
root -l -b -q fitRecoil_master.C\(3,3,1,1,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWpemet_sim_triple&
root -l -b -q fitRecoil_master.C\(0,0,1,1,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWpemet_sim_CB&