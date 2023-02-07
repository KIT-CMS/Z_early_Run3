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

# weight the luminosity of the MC in fb-1
# def lumi_weight(era, runPlot, totalLumi):
#     if era == "2016":
#         lumi = "35.87"
#     elif era == "2017":
#         lumi = "41.529"
#     elif era == "2018":
#         # lumi = "59.8"
#         lumi = "31.75"  # Run2018D
#         # lumi = "0.001" if runPlot else str(totalLumi)
#     elif era == "2022":
#         lumi = "0.001" if runPlot else str(totalLumi)
#     else:
#         raise ValueError("Given era {} not defined.".format(era))
#     return ("{} * 1000.0".format(lumi), "lumi")

def lumi_weight(era, runPlot, totalLumi, run_list, run_lumi, norm1invpb):
    lumi = None
    if era == "2022":            
        if runPlot:
            if norm1invpb:
                lumi = "0.001"
            else:
                sum_lumi = 0
                for run_numb in run_list:
                    sum_lumi+=run_lumi[run_numb]
                lumi = sum_lumi * 1.e-3
        else:
            lumi = str(totalLumi)
    
    return ("{} * 1000.0".format(lumi), "lumi")

# weights for the monte carlo seperated by channel
def MC_base_process_selection(channel, era, runPlot, totalLumi, run_list, run_lumi, norm1invpb, applySF):
    if channel in ["mmet"]:
        MC_base_process_weights = [
            ("genweight*sumwWeight*crossSectionPerEventWeight", "normWeight"),
            # ("puweight", "puweight"),
            ("1.", "LHEScaleWeight"),
            ("1.", "LHEPdfWeight"),
            lumi_weight(era, runPlot, totalLumi, run_list, run_lumi, norm1invpb),
        ]
        if applySF:
            MC_base_process_weights += [("sf_trk*sf_sta*sf_id*sf_iso*sf_trg", "sfWeight")]
        return Selection(name="MC base", weights=MC_base_process_weights)
    elif channel in ["emet"]:
        MC_base_process_weights = [
            ("genweight*sumwWeight*crossSectionPerEventWeight", "normWeight"),
            # ("puweight", "puweight"),
            # ("id_wgt_ele_wpmedium_1", "idweight"),
            ("1.", "LHEScaleWeight"),
            ("1.", "LHEPdfWeight"),
            lumi_weight(era, runPlot, totalLumi, run_list, run_lumi, norm1invpb),
        ]
        return Selection(name="MC base", weights=MC_base_process_weights)
    elif channel in ["mm"]:
        MC_base_process_weights = [
            ("genweight*sumwWeight*crossSectionPerEventWeight", "normWeight"),
            # ("puweight", "puweight"),
            ("1.", "LHEScaleWeight"),
            ("1.", "LHEPdfWeight"),
            lumi_weight(era, runPlot, totalLumi, run_list, run_lumi, norm1invpb),
        ]
        if applySF:
            MC_base_process_weights += [("sf_trk*sf_sta*sf_id*sf_iso*sf_trg", "sfWeight")]
            # MC_base_process_weights += [("sf_trk*sf_sta*sf_id*sf_iso", "sfWeight")]  # Z counting
        return Selection(name="MC base", weights=MC_base_process_weights)
    elif channel in ["ee"]:
        MC_base_process_weights = [
            ("genweight*sumwWeight*crossSectionPerEventWeight", "normWeight"),
            # ("puweight", "puweight"),
            # ("id_wgt_ele_wpmedium_1*id_wgt_ele_wpmedium_2", "idweight"),
            ("1.", "LHEScaleWeight"),
            ("1.", "LHEPdfWeight"),
            lumi_weight(era, runPlot, totalLumi, run_list, run_lumi, norm1invpb),
        ]
        return Selection(name="MC base", weights=MC_base_process_weights)


def DY_process_selection(channel, era, runPlot, totalLumi, run_list, run_lumi, norm1invpb, applySF):
    DY_process_weights = MC_base_process_selection(channel, era, runPlot, totalLumi, run_list, run_lumi, norm1invpb, applySF).weights
    cuts = [
        ("(is_dy_tt != 1)", "TauFilter")
    ]
    return Selection(name="DY", cuts=cuts, weights=DY_process_weights)

def W_process_selection(channel, era, runPlot, totalLumi, run_list, run_lumi, norm1invpb, applySF):
    W_process_weights = MC_base_process_selection(channel, era, runPlot, totalLumi, run_list, run_lumi, norm1invpb, applySF).weights
    cuts = [
        ("(gen_match_1 != 15)", "TauFilter")
    ]
    return Selection(name="W", cuts=cuts, weights=W_process_weights)

def TT_process_selection(channel, era, runPlot, totalLumi, run_list, run_lumi, norm1invpb, applySF):
    TT_process_weights = MC_base_process_selection(channel, era, runPlot, totalLumi, run_list, run_lumi, norm1invpb, applySF).weights
    return Selection(name="TT", weights=TT_process_weights)

def ST_process_selection(channel, era, runPlot, totalLumi, run_list, run_lumi, norm1invpb, applySF):
    ST_process_weights = MC_base_process_selection(channel, era, runPlot, totalLumi, run_list, run_lumi, norm1invpb, applySF).weights
    return Selection(name="ST", weights=ST_process_weights)

def VV_process_selection(channel, era, runPlot, totalLumi, run_list, run_lumi, norm1invpb, applySF):
    VV_process_weights = MC_base_process_selection(channel, era, runPlot, totalLumi, run_list, run_lumi, norm1invpb, applySF).weights
    return Selection(name="VV", weights=VV_process_weights)

def DYtau_process_selection(channel, era, runPlot, totalLumi, run_list, run_lumi, norm1invpb, applySF):
    DYtau_process_weights = MC_base_process_selection(channel, era, runPlot, totalLumi, run_list, run_lumi, norm1invpb, applySF).weights
    cuts = [
        ("(is_dy_tt == 1)", "TauFilter")
    ]
    return Selection(name="DYtau", cuts=cuts, weights=DYtau_process_weights)

def Wtau_process_selection(channel, era, runPlot, totalLumi, run_list, run_lumi, norm1invpb, applySF):
    Wtau_process_weights = MC_base_process_selection(channel, era, runPlot, totalLumi, run_list, run_lumi, norm1invpb, applySF).weights
    cuts = [
        ("(gen_match_1 == 15)", "TauFilter")
    ]
    return Selection(name="Wtau", cuts=cuts, weights=Wtau_process_weights)