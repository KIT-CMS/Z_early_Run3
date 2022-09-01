from pickletools import float8
from ntuple_processor.utils import Selection

# fill in with appropriate cuts!!
def channel_selection(channel, era):
    if "mmet" in channel:
        cuts = [
            ("(pt_1>25. && abs(eta_1) < 2.4)", "acceptance"),
            ("(trg_single_mu24_1)", "trg_matching"),
            ("(iso_1 < 0.15)", "lepton_iso"),
            ("(extramuon_veto == 1)", "lepton_veto"),

            # ("(sqrt(2.*pt_1*met_uncorrected*(1.-cos(phi_1 - metphi_uncorrected)))>40.)", "mt_cut"),
        ]
        return Selection(name="mmet", cuts=cuts)
    elif "emet" in channel:
        cuts = [
            ("(pt_1>29. && abs(eta_1) < 2.5)", "acceptance"),
            ("(trg_single_ele27_1)", "trg_matching"),
            ("(((abs(eta_1 + deltaetaSC_1) <= 1.479) && (iso_1 < 0.0478+0.506/pt_1)) || ((abs(eta_1 + deltaetaSC_1) > 1.479) && (iso_1 < 0.0658+0.963/pt_1)))", "lepton_iso"),
            ("(extraelec_veto == 1)", "lepton_veto"),

            # ("(sqrt(2.*pt_1*met_uncorrected*(1.-cos(phi_1 - metphi_uncorrected)))>40.)", "mt_cut"),
        ]
        return Selection(name="emet", cuts=cuts)
    elif "mm" in channel:
        cuts = [
            ("(pt_1>25. && pt_2>25. && abs(eta_1) < 2.4 && abs(eta_2) < 2.4)", "acceptance"),
            ("(q_1*q_2 < 0)", "opposite_charge"),
            ("(m_vis > 60. && m_vis < 120.)", "Z_mass_window"),
            ("(trg_single_mu24_1 || trg_single_mu24_2)", "trg_matching"),
        ]
        return Selection(name="mm", cuts=cuts)
    elif "ee" in channel:
        cuts = [
            ("(pt_1>29. && pt_2>29. && abs(eta_1) < 2.5 && abs(eta_2) < 2.5)", "acceptance"),
            ("(q_1*q_2 < 0)", "opposite_charge"),
            ("(m_vis > 60. &&  m_vis < 120.)", "Z_mass_window"),
            ("(trg_single_ele27_1 || trg_single_ele27_1)", "trg_matching"),
        ]
        return Selection(name="ee", cuts=cuts)