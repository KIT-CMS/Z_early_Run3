from ntuple_processor import Histogram
import os 
import numpy as np

#corresponding lists are in the order [lower_bound, upper_bound, number_of_bins]
uniformBinningDict = {
  'pt': [25, 200, 40],
  'ptOverMt': [0, 120, 40],
  'eta': [(-2.4), 2.4, 40],
  'phi': [(-np.pi), (np.pi), 40],
  'm_vis': [60, 120, 40],
  'pt_vis': [0, 200, 40],
  'mt_uncorrected': [40, 200, 40],
  'nStations': [0, 6, 6],
  'nTrackerLayers': [0, 20, 20],
  'deltaetaSC': [(-0.05), 0.05, 40],
  'eInvMinusPInv': [(-0.005), (0.005), 40],
  'hoe': [0, 0.05, 40],
  'scEtOverPt': [(-0.2), 0.05, 40],
  'sieie': [],
  'lostHits': [0, 2, 2],
  'met_uncorrected': [0, 200, 40],
  'metphi_uncorrected': [(-np.pi), (np.pi), 40],
  'pfmet_uncorrected': [0, 200, 40],
  'pfmetphi_uncorrected':  [(-np.pi), (np.pi), 40]
  }


#function which allows you to add new variables not in the ntuple
def renamingVariable(var, channel):
  selection = "1"
  base = var.replace("_1", "").replace("_2", "").replace("_barrel", "").replace("_endcap", "").replace("_pos", "").replace("_neg", "").replace("_nv015", "").replace("_nv1530", "").replace("_nv3045", "").replace("_isoSR", "").replace("_iso5", "").replace("_iso6", "").replace("_iso7", "").replace("_iso8", "").replace("_iso9", "").replace("_iso10", "").replace("_iso11", "").replace("_iso12", "").replace("_iso13", "")
  idx = "2" if ("_2" in var) else "1"

  #number of vertices
  if "_nv015" in var: selection = selection + "*(0 < npvGood && npvGood <= 15)"
  elif "_nv1530" in var: selection = selection + "*(15 < npvGood && npvGood <= 30)"
  elif "_nv3045" in var: selection = selection + "*(30 < npvGood && npvGood < 45)"

  #mT
  if "mt_uncorrected" in var and "mm" in channel: base = "(sqrt(2.*pt_rc_1*met_uncorrected*(1.-cos(phi_1 - metphi_uncorrected))))"
  elif "mt_uncorrected" in var: base = "(sqrt(2.*pt_1*met_uncorrected*(1.-cos(phi_1 - metphi_uncorrected))))"

  #ptOvermT
  if "ptOverMt" in var: base = "(pt_rc_1/(sqrt(2.*pt_rc_1*met_uncorrected*(1.-cos(phi_1 - metphi_uncorrected)))))"

  #charge
  if "_pos" in var: selection = selection + "*(q_1 > 0)"
  elif "_neg" in var: selection = selection + "*(q_1 < 0)"

  #signal isolation
  if "_isoSR" in var and channel == "mmet": selection = selection + "*(iso_1 < 0.15)"
  elif "_isoSR" in var and "_barrel" in var: selection = selection + "*(iso_1 < 0.0478+0.506/pt_1)"
  elif "_isoSR" in var and "_endcap" in var: selection = selection + "*(iso_1 < 0.0658+0.963/pt_1)"

  #background isolation
  if "_iso5" in var: selection = selection + "*(iso_1 > 0.20 && iso_1 < 0.25)"
  elif "_iso6" in var: selection = selection + "*(iso_1 > 0.25 && iso_1 < 0.30)"
  elif "_iso7" in var: selection = selection + "*(iso_1 > 0.30 && iso_1 < 0.35)"
  elif "_iso8" in var: selection = selection + "*(iso_1 > 0.35 && iso_1 < 0.40)"
  elif "_iso9" in var: selection = selection + "*(iso_1 > 0.40 && iso_1 < 0.45)"
  elif "_iso10" in var: selection = selection + "*(iso_1 > 0.45 && iso_1 < 0.50)"
  elif "_iso11" in var: selection = selection + "*(iso_1 > 0.50 && iso_1 < 0.55)"
  elif "_iso12" in var: selection = selection + "*(iso_1 > 0.55 && iso_1 < 0.60)"
  elif "_iso13" in var: selection = selection + "*(iso_1 > 0.60 && iso_1 < 0.65)"

  #barrel and endcap for electrons
  if channel == "ee" or channel == "emet":
    #barrel or endcap
    if "_barrel" in var: selection = selection + "*(abs(eta_{idx} + deltaetaSC_{idx}) <= 1.479)".format(idx=idx)
    elif "_endcap" in var: selection = selection + "*(abs(eta_{idx} + deltaetaSC_{idx}) > 1.479)".format(idx=idx)

  #barrel and endcap for muons
  elif channel == "mm" or channel == "mmet":
    if "_barrel" in var: selection = selection + "*(abs(eta_{idx}) < 1.2)".format(idx=idx)
    elif "_endcap" in var: selection = selection + "*(abs(eta_{idx}) > 1.2)".format(idx=idx)

  full_var = "{base}*({sel}) -1.e9*(!({sel}))".format(base=base, sel=selection)
  #assert full_var != None
  return full_var

#function that seperates your variables by quantities like number of primary vertices, charge, etc.
def seperateVariables(input_list):
  variable_list = list()

  charge_list=["_pos", "_neg"]
  eta_list=["_barrel", "_endcap"]
  iso_list=["_isoSR", "_iso5", "_iso6", "_iso7", "_iso8", "_iso9", "_iso10", "_iso11", "_iso12", "_iso13"]
  npv_list=["_npv015", "_npv1530", "_npv3045"]

  for var in input_list:
    for charge in charge_list:
      for eta in eta_list:
        for iso in iso_list:
          variable_list.append(var+charge+eta+iso)

  return variable_list


def get_control_binning(channel, variable_list):
  control_binning = {channel:dict()}
  for i in range(len(variable_list)):
    var_name = renamingVariable(variable_list[i], channel)
    histo_bin_list=list()

    #makes the bins and bounds
    binning_base = variable_list[i].replace("_1", "").replace("_2", "").replace("_barrel", "").replace("_endcap", "").replace("_pos", "").replace("_neg", "").replace("_nv015", "").replace("_nv1530", "").replace("_nv3045", "").replace("_isoSR", "").replace("_iso5", "").replace("_iso6", "").replace("_iso7", "").replace("_iso8", "").replace("_iso9", "").replace("_iso10", "").replace("_iso11", "").replace("_iso12", "").replace("_iso13", "")
    
    for k in range(0, uniformBinningDict[binning_base][2]):
      var_bin_list=uniformBinningDict[binning_base][0]+k*(uniformBinningDict[binning_base][1] - uniformBinningDict[binning_base][0])/(uniformBinningDict[binning_base][2])
      histo_bin_list.append(var_bin_list)
    
    #control binning
    control_binning[channel][variable_list[i]]=Histogram(variable_list[i], var_name, histo_bin_list)
  return control_binning
