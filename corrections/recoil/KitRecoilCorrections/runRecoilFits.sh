#! /bin/bash

OUT=Run3V06_ptuncorr_outputs
IN=/storage/9/jdriesch/earlyrun3/samples/Run3V06
METVAR=met_uncorrected
METPHIVAR=metphi_uncorrected
LUMI=5.035650234254

# ####################################################################################################
# #--------------------------------------------------------------------------------------------------#
# #-            Z --> mm                                                                           -#
# # -------------------------------------------------------------------------------------------------#
# ####################################################################################################

#Data
root -l -b -q fitRecoil.C\(3,3,1,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZmm_data_triple&
root -l -b -q fitRecoil.C\(2,2,1,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZmm_data_double&
# root -l -b -q fitRecoil.C\(0,0,1,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZmm_data_CB&

#Data and background
root -l -b -q fitRecoil.C\(3,3,0,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZmm_dataBckg_triple&
root -l -b -q fitRecoil.C\(2,2,0,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZmm_dataBckg_double&
# root -l -b -q fitRecoil.C\(0,0,0,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZmm_dataBckg_CB&

root -l -b -q fitRecoil.C\(3,3,0,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},0\) &>logZmm_dataBckg_zrap0&
root -l -b -q fitRecoil.C\(3,3,0,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},1\) &>logZmm_dataBckg_zrap1&
root -l -b -q fitRecoil.C\(3,3,0,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},2\) &>logZmm_dataBckg_zrap2&

#Simulation
root -l -b -q fitRecoil.C\(3,3,1,0,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZmm_sim_triple&
root -l -b -q fitRecoil.C\(2,2,1,0,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZmm_sim_double&
# root -l -b -q fitRecoil.C\(0,0,1,0,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZmm_sim_CB&

root -l -b -q fitRecoil.C\(3,3,1,0,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},0\) &>logZmm_sim_zrap0&
root -l -b -q fitRecoil.C\(3,3,1,0,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},1\) &>logZmm_sim_zrap1&
# root -l -b -q fitRecoil.C\(3,3,1,0,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},2\) &>logZmm_sim_zrap2&

# ####################################################################################################
# #--------------------------------------------------------------------------------------------------#
# #-            Z --> ee                                                                           -#
# # -------------------------------------------------------------------------------------------------#
# ####################################################################################################

# #Data
# root -l -b -q fitRecoil.C\(3,3,1,1,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZee_data_triple&
# root -l -b -q fitRecoil.C\(0,0,1,1,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZee_data_CB&

# #Data and background
# root -l -b -q fitRecoil.C\(3,3,0,1,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZee_dataBckg_triple&
# root -l -b -q fitRecoil.C\(0,0,0,1,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZee_dataBckg_CB&

# #Simulation
# root -l -b -q fitRecoil.C\(3,3,1,1,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZee_sim_triple&
# root -l -b -q fitRecoil.C\(0,0,1,1,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logZee_sim_CB&

# ####################################################################################################
# #--------------------------------------------------------------------------------------------------#
# #-            W+ --> mmet                                                                           -#
# # -------------------------------------------------------------------------------------------------#
# ####################################################################################################

# Simulation
root -l -b -q fitRecoil.C\(3,3,1,0,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWpmmet_sim_triple&
root -l -b -q fitRecoil.C\(2,2,1,0,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWpmmet_sim_double&
# root -l -b -q fitRecoil.C\(0,0,1,0,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWpmmet_sim_CB&

root -l -b -q fitRecoil.C\(3,3,1,0,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},0\) &>logWpmmet_sim_zrap0&
root -l -b -q fitRecoil.C\(3,3,1,0,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},1\) &>logWpmmet_sim_zrap1&
root -l -b -q fitRecoil.C\(3,3,1,0,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},2\) &>logWpmmet_sim_zrap2&

# ####################################################################################################
# #--------------------------------------------------------------------------------------------------#
# #-            W- --> mmet                                                                           -#
# # -------------------------------------------------------------------------------------------------#
# ####################################################################################################

# Simulation
root -l -b -q fitRecoil.C\(3,3,1,0,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWnmmet_sim_triple&
root -l -b -q fitRecoil.C\(2,2,1,0,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWnmmet_sim_double&
# root -l -b -q fitRecoil.C\(0,0,1,0,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWnmmet_sim_CB&

root -l -b -q fitRecoil.C\(3,3,1,0,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},0\) &>logWnmmet_sim_zrap0&
root -l -b -q fitRecoil.C\(3,3,1,0,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},1\) &>logWnmmet_sim_zrap1&
root -l -b -q fitRecoil.C\(3,3,1,0,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},2\) &>logWnmmet_sim_zrap2&

# ####################################################################################################
# #--------------------------------------------------------------------------------------------------#
# #-            W- --> emet                                                                           -#
# # -------------------------------------------------------------------------------------------------#
# ####################################################################################################

# #Simulation
# root -l -b -q fitRecoil.C\(3,3,1,1,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWnemet_sim_triple&
# root -l -b -q fitRecoil.C\(0,0,1,1,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWnemet_sim_CB&

# ####################################################################################################
# #--------------------------------------------------------------------------------------------------#
# #-            W+ --> emet                                                                           -#
# # -------------------------------------------------------------------------------------------------#
# ####################################################################################################

# #Simulation
# root -l -b -q fitRecoil.C\(3,3,1,1,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWpemet_sim_triple&
# root -l -b -q fitRecoil.C\(0,0,1,1,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logWpemet_sim_CB&
