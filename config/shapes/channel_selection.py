from pickletools import float8
from dataclasses import asdict
from ntuple_processor.utils import Selection

def channel_selection(channel, era, doQCD = False, n_trg_match_mm = None, corr_postfix = ""):
    name_base = channel + corr_postfix
    if "mmet" in channel:
        cuts = [
            (f"(pt_1{corr_postfix}>25. && abs(eta_1) < 2.4)", "acceptance"),
            ("(trg_single_mu24_1)", "trg_matching"),
            ("(extramuon_veto == 1)", "lepton_veto"),
            ("(sqrt(2.*pt_1_corr*pfmet_corr*(1.-cos(phi_1 - pfmetphi_corr)))) > -1.", "pfmtcut"),
        ]
        if not doQCD:
            cuts.append(
                ("(iso_1 < 0.15)", "lepton_iso")
            )
        return Selection(name=name_base, cuts=cuts)
    elif "emet" in channel:
        cuts = [
            (f"(pt_1{corr_postfix}>37. && abs(eta_1) < 2.5)", "acceptance"),
            ("(trg_single_ele27_1 || trg_single_ele28_1 || trg_single_ele32_1 || trg_single_ele35_1)", "trg_matching"),
            ("(extraelec_veto == 1)", "lepton_veto"),
            ("(sqrt(2.*pt_1_corr*pfmet_corr*(1.-cos(phi_1 - pfmetphi_corr)))) > -1.", "pfmtcut"),
        ]
        if not doQCD:
            cuts.append(
                ("(((abs(eta_1 + deltaetaSC_1) <= 1.479) && (iso_1 < 0.0478+0.506/pt_1)) || ((abs(eta_1 + deltaetaSC_1) > 1.479) && (iso_1 < 0.0658+0.963/pt_1)))", "lepton_iso")
            )
        return Selection(name=name_base, cuts=cuts)
    elif "mm" in channel:
        cuts = [
            (f"(pt_1{corr_postfix}>25. && pt_2{corr_postfix}>25. && abs(eta_1) < 2.4 && abs(eta_2) < 2.4)", "acceptance"),
            ("(q_1*q_2 < 0)", "opposite_charge"),
            (f"(m_vis{corr_postfix} > 60. && m_vis{corr_postfix} < 120.)", "Z_mass_window"),
        ]
        name_mm = name_base
        if n_trg_match_mm != None:
            assert (n_trg_match_mm in [2, 1, 0, -1])
            if n_trg_match_mm == 2:
                cuts.append(("(trg_single_mu24_1 && trg_single_mu24_2)", "trg_matching"))
                name_mm = name_base+"-HLT2"
            elif n_trg_match_mm == 1:
                cuts.append(("((trg_single_mu24_1 && !trg_single_mu24_2) || (!trg_single_mu24_1 && trg_single_mu24_2))", "trg_matching"))
                name_mm = name_base+"-HLT1"
            elif n_trg_match_mm == 0:
                cuts.append(("(!trg_single_mu24_1 && !trg_single_mu24_2)", "trg_matching"))
                name_mm = name_base+"-HLT0"
            elif n_trg_match_mm == -1:
                cuts.append(("(trg_single_mu24_1 || trg_single_mu24_2)", "trg_matching"))
        else:
            cuts.append(
                ("(trg_single_mu24_1 || trg_single_mu24_2)", "trg_matching")
            )
        return Selection(name=name_mm, cuts=cuts)
    elif "ee" in channel:
        cuts = [
            (f"(pt_1{corr_postfix}>37. && pt_2{corr_postfix}>25. && abs(eta_1) < 2.5 && abs(eta_2) < 2.5)", "acceptance"),
            ("(q_1*q_2 < 0)", "opposite_charge"),
            (f"(m_vis{corr_postfix} > 60. &&  m_vis{corr_postfix} < 120.)", "Z_mass_window"),
            ("(trg_single_ele27_1 || trg_single_ele28_1 || trg_single_ele32_1 || trg_single_ele35_1)", "trg_matching"),
        ]
        return Selection(name=name_base, cuts=cuts)
