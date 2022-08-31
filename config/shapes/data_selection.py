from ntuple_processor.utils import Selection
import os

"""Base processes

List of base processes, mostly containing only weights:
    - triggerweight
    - triggerweight_emb
    - tau_by_iso_id_weight
    - ele_hlt_Z_vtx_weight
    - ele_reco_weight
    - aiso_muon_correction
    - lumi_weight
    - MC_base_process_selection
    - DY_base_process_selection
    - TT_process_selection
    - VV_process_selection
    - W_process_selection
    - HTT_base_process_selection
    - HTT_process_selection
    - HWW_process_selection
"""

def lumi_weight(era, run_numb, run_lumi):
    return ("{}".format(1./(run_lumi[run_numb])), "lumi")

def data_base_process_selection(channel, era, run_numb, run_lumi):
    if channel in ["mmet"]:
        data_base_process_weights = [
            lumi_weight(era, run_numb, run_lumi),
        ]
        return Selection(name="data base", weights=data_base_process_weights)
    elif channel in ["emet"]:
        data_base_process_weights = [
            lumi_weight(era, run_numb, run_lumi),
        ]
        return Selection(name="data base", weights=data_base_process_weights)
    elif channel in ["mm"]:
        data_base_process_weights = [
            lumi_weight(era, run_numb, run_lumi),
        ]
        return Selection(name="data base", weights=data_base_process_weights)
    elif channel in ["ee"]:
        data_base_process_weights = [
            lumi_weight(era, run_numb, run_lumi),
        ]
        return Selection(name="data base", weights=data_base_process_weights)

def group_lumi_weight(era, run_list, run_lumi):
    sum_lumi = 0
    for run_numb in run_list:
        sum_lumi+=run_lumi[run_numb]
    print("Summed Luminosity of Run Group: ", sum_lumi)
    return ("{}".format(1./sum_lumi), "lumi")

def data_base_group_process_selection(channel, era, run_list, run_lumi):
    if channel in ["mmet"]:
        data_base_process_weights = [
            group_lumi_weight(era, run_list, run_lumi),
        ]
        return Selection(name="data base", weights=data_base_process_weights)
    elif channel in ["emet"]:
        data_base_process_weights = [
            group_lumi_weight(era, run_list, run_lumi),
        ]
        return Selection(name="data base", weights=data_base_process_weights)
    elif channel in ["mm"]:
        data_base_process_weights = [
            group_lumi_weight(era, run_list, run_lumi),
        ]
        return Selection(name="data base", weights=data_base_process_weights)
    elif channel in ["ee"]:
        data_base_process_weights = [
            group_lumi_weight(era, run_list, run_lumi),
        ]
        return Selection(name="data base", weights=data_base_process_weights)

def data_process_selection(channel, era, run_numb, run_lumi):
    data_process_weights = data_base_process_selection(channel, era, run_numb, run_lumi).weights
    return Selection(name="data", weights=data_process_weights)

def data_group_process_selection(channel, era, run_list, run_lumi):
    data_process_weights = data_base_group_process_selection(channel, era, run_list, run_lumi).weights
    return Selection(name="data", weights=data_process_weights)

