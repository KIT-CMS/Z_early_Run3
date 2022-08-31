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

#weight the luminosity of the MC in fb-1
def lumi_weight(era, runPlot, totalLumi):
    if era == "2016":
        lumi = "35.87"
    elif era == "2017":
        lumi = "41.529"
    elif era == "2018":
        # FIXME: testing with Run2018D only
        # lumi = "31.75"
        # lumi = "13.98"  # Run2018A
        lumi = "0.001" if runPlot else str(totalLumi)
        #lumi = "59.8"
    else:
        raise ValueError("Given era {} not defined.".format(era))
    return ("{} * 1000.0".format(lumi), "lumi")

#weights for the monte carlo seperated by channel
def MC_base_process_selection(channel, era, runPlot, totalLumi):
    if channel in ["mmet"]:
        MC_base_process_weights = [
            ("genweight*sumwWeight*crossSectionPerEventWeight", "normWeight"),
            ("puweight", "puweight"),
            ("id_wgt_mu_1", "idweight"),
            ("iso_wgt_mu_1", "isoweight"),
            lumi_weight(era, runPlot, totalLumi),
        ]
        return Selection(name="MC base", weights=MC_base_process_weights)
    elif channel in ["emet"]:
        MC_base_process_weights = [
            ("genweight*sumwWeight*crossSectionPerEventWeight", "normWeight"),
            ("puweight", "puweight"),
            ("id_wgt_ele_wpmedium_1", "idweight"),
            lumi_weight(era, runPlot, totalLumi),
        ]
        return Selection(name="MC base", weights=MC_base_process_weights)
    elif channel in ["mm"]:
        MC_base_process_weights = [
            ("genweight*sumwWeight*crossSectionPerEventWeight", "normWeight"),
            ("puweight", "puweight"),
            ("id_wgt_mu_1*id_wgt_mu_2", "idweight"),
            ("iso_wgt_mu_1*iso_wgt_mu_2", "isoweight"),
            lumi_weight(era, runPlot, totalLumi),
        ]
        return Selection(name="MC base", weights=MC_base_process_weights)
    elif channel in ["ee"]:
        MC_base_process_weights = [
            ("genweight*sumwWeight*crossSectionPerEventWeight", "normWeight"),
            ("puweight", "puweight"),
            ("id_wgt_ele_wpmedium_1*id_wgt_ele_wpmedium_2", "idweight"),
            lumi_weight(era, runPlot, totalLumi),
        ]
        return Selection(name="MC base", weights=MC_base_process_weights)


def DY_process_selection(channel, era, runPlot, totalLumi):
    DY_process_weights = MC_base_process_selection(channel, era, runPlot, totalLumi).weights
    DY_process_weights.append(
        ("6372.24/6077.22", "XS13p6TeV")
    )
    cuts = None
    if channel in ["mm", "ee"]:
        cuts = [
            ("(gen_match_1 != 15 && gen_match_2 != 15)", "tautauFilter")
        ]
    return Selection(name="DY", cuts=cuts, weights=DY_process_weights)


def TT_process_selection(channel, era, runPlot, totalLumi):
    TT_process_weights = MC_base_process_selection(channel, era, runPlot, totalLumi).weights
    TT_process_weights.append(
        ("920.92/831.76", "XS13p6TeV")
    )
    return Selection(name="TT", weights=TT_process_weights)


def VV_process_selection(channel, era, runPlot, totalLumi):
    VV_process_weights = MC_base_process_selection(channel, era, runPlot, totalLumi).weights
    return Selection(name="VV", weights=VV_process_weights)


def W_process_selection(channel, era, runPlot, totalLumi):
    W_process_weights = MC_base_process_selection(channel, era, runPlot, totalLumi).weights
    return Selection(name="W", weights=W_process_weights)


def ZL_process_selection(channel):
    return Selection(name="ZL", weights=[("1.0", "zl_weight")])


def TTL_process_selection(channel):
    return Selection(name="TTL", weights=[("1.0", "ttl_weight")])


def VVL_process_selection(channel):
    return Selection(name="VVL", weights=[("1.0", "vvl_weight")])
