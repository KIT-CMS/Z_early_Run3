#! /bin/bash

OUT=Run3V02_outputs
IN=/ceph/moh/CROWN_samples/Run3V02
METVAR=pfmet_uncorrected
METPHIVAR=pfmetphi_uncorrected
LUMI=4.844307925632

# ####################################################################################################
# #--------------------------------------------------------------------------------------------------#
# #-            Z --> mm                                                                           -#
# # -------------------------------------------------------------------------------------------------#
# ####################################################################################################

#Data
# root -l -b -q fitRecoil.C\(3,3,1,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfZmm_data_triple&
# root -l -b -q fitRecoil.C\(2,2,1,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfZmm_data_double&
# root -l -b -q fitRecoil.C\(0,0,1,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfZmm_data_CB&

#Data and background
# root -l -b -q fitRecoil.C\(3,3,0,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfZmm_dataBckg_triple&
# root -l -b -q fitRecoil.C\(2,2,0,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfZmm_dataBckg_double&
# root -l -b -q fitRecoil.C\(0,0,0,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfZmm_dataBckg_CB&

root -l -b -q fitRecoil.C\(3,3,0,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},0\) &>logpfZmm_dataBckg_zrap0&
root -l -b -q fitRecoil.C\(3,3,0,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},1\) &>logpfZmm_dataBckg_zrap1&
root -l -b -q fitRecoil.C\(3,3,0,0,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},2\) &>logpfZmm_dataBckg_zrap2&

#Simulation
# root -l -b -q fitRecoil.C\(3,3,1,0,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfZmm_sim_triple&
# root -l -b -q fitRecoil.C\(2,2,1,0,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfZmm_sim_double&
# root -l -b -q fitRecoil.C\(0,0,1,0,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfZmm_sim_CB&

root -l -b -q fitRecoil.C\(3,3,1,0,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},0\) &>logpfZmm_sim_zrap0&
root -l -b -q fitRecoil.C\(3,3,1,0,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},1\) &>logpfZmm_sim_zrap1&
root -l -b -q fitRecoil.C\(3,3,1,0,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},2\) &>logpfZmm_sim_zrap2&

# ####################################################################################################
# #--------------------------------------------------------------------------------------------------#
# #-            Z --> ee                                                                           -#
# # -------------------------------------------------------------------------------------------------#
# ####################################################################################################

# #Data
# root -l -b -q fitRecoil.C\(3,3,1,1,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfZee_data_triple&
# root -l -b -q fitRecoil.C\(0,0,1,1,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfZee_data_CB&

# #Data and background
# root -l -b -q fitRecoil.C\(3,3,0,1,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfZee_dataBckg_triple&
# root -l -b -q fitRecoil.C\(0,0,0,1,0,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfZee_dataBckg_CB&

# #Simulation
# root -l -b -q fitRecoil.C\(3,3,1,1,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfZee_sim_triple&
# root -l -b -q fitRecoil.C\(0,0,1,1,1,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfZee_sim_CB&

# ####################################################################################################
# #--------------------------------------------------------------------------------------------------#
# #-            W+ --> mmet                                                                           -#
# # -------------------------------------------------------------------------------------------------#
# ####################################################################################################

# Simulation
# root -l -b -q fitRecoil.C\(3,3,1,0,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfWpmmet_sim_triple&
# root -l -b -q fitRecoil.C\(2,2,1,0,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfWpmmet_sim_double&
# root -l -b -q fitRecoil.C\(0,0,1,0,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfWpmmet_sim_CB&

root -l -b -q fitRecoil.C\(3,3,1,0,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},0\) &>logpfWpmmet_sim_zrap0&
root -l -b -q fitRecoil.C\(3,3,1,0,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},1\) &>logpfWpmmet_sim_zrap1&
root -l -b -q fitRecoil.C\(3,3,1,0,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},2\) &>logpfWpmmet_sim_zrap2&

# ####################################################################################################
# #--------------------------------------------------------------------------------------------------#
# #-            W- --> mmet                                                                           -#
# # -------------------------------------------------------------------------------------------------#
# ####################################################################################################

# Simulation
# root -l -b -q fitRecoil.C\(3,3,1,0,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfWnmmet_sim_triple&
# root -l -b -q fitRecoil.C\(2,2,1,0,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfWnmmet_sim_double&
# root -l -b -q fitRecoil.C\(0,0,1,0,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfWnmmet_sim_CB&

root -l -b -q fitRecoil.C\(3,3,1,0,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},0\) &>logpfWnmmet_sim_zrap0&
root -l -b -q fitRecoil.C\(3,3,1,0,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},1\) &>logpfWnmmet_sim_zrap1&
root -l -b -q fitRecoil.C\(3,3,1,0,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI},2\) &>logpfWnmmet_sim_zrap2&

# ####################################################################################################
# #--------------------------------------------------------------------------------------------------#
# #-            W- --> emet                                                                           -#
# # -------------------------------------------------------------------------------------------------#
# ####################################################################################################

# #Simulation
# root -l -b -q fitRecoil.C\(3,3,1,1,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfWnemet_sim_triple&
# root -l -b -q fitRecoil.C\(0,0,1,1,3,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfWnemet_sim_CB&

# ####################################################################################################
# #--------------------------------------------------------------------------------------------------#
# #-            W+ --> emet                                                                           -#
# # -------------------------------------------------------------------------------------------------#
# ####################################################################################################

# #Simulation
# root -l -b -q fitRecoil.C\(3,3,1,1,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfWpemet_sim_triple&
# root -l -b -q fitRecoil.C\(0,0,1,1,2,\"${IN}\",\"${METVAR}\",\"${METPHIVAR}\",\"${OUT}\",${LUMI}\) &>logpfWpemet_sim_CB&