from ntuple_processor.variations import ReplaceWeight, ReplaceCut

# mu_sf_weight = [
#     ReplaceWeight("SFTrkUp", "sfWeight", ("sf_trk_up*sf_sta*sf_id*sf_iso*sf_trg", "sfWeightTrkUp")),
#     ReplaceWeight("SFTrkDn", "sfWeight", ("sf_trk_dn*sf_sta*sf_id*sf_iso*sf_trg", "sfWeightTrkDn")),

#     ReplaceWeight("SFStaUp", "sfWeight", ("sf_trk*sf_sta_up*sf_id*sf_iso*sf_trg", "sfWeightStaUp")),
#     ReplaceWeight("SFStaDn", "sfWeight", ("sf_trk*sf_sta_dn*sf_id*sf_iso*sf_trg", "sfWeightStaDn")),

#     ReplaceWeight("SFIDUp", "sfWeight", ("sf_trk*sf_sta*sf_id_up*sf_iso*sf_trg", "sfWeightIDUp")),
#     ReplaceWeight("SFIDDn", "sfWeight", ("sf_trk*sf_sta*sf_id_dn*sf_iso*sf_trg", "sfWeightIDDn")),

#     ReplaceWeight("SFIsoUp", "sfWeight", ("sf_trk*sf_sta*sf_id*sf_iso_up*sf_trg", "sfWeightIsoUp")),
#     ReplaceWeight("SFIsoDn", "sfWeight", ("sf_trk*sf_sta*sf_id*sf_iso_dn*sf_trg", "sfWeightIsoDn")),

#     ReplaceWeight("SFTrgUp", "sfWeight", ("sf_trk*sf_sta*sf_id*sf_iso*sf_trg_up", "sfWeightTrgUp")),
#     ReplaceWeight("SFTrgDn", "sfWeight", ("sf_trk*sf_sta*sf_id*sf_iso*sf_trg_dn", "sfWeightTrgDn")),
# ]
mu_sf_weight = [
    ReplaceWeight("SFTrkUp", "sfWeight", ("sf_sta*sf_id*sf_iso*sf_trg", "sfWeightTrkUp")),
    ReplaceWeight("SFTrkDn", "sfWeight", ("sf_sta*sf_id*sf_iso*sf_trg", "sfWeightTrkDn")),

    ReplaceWeight("SFStaUp", "sfWeight", ("sf_sta_up*sf_id*sf_iso*sf_trg", "sfWeightStaUp")),
    ReplaceWeight("SFStaDn", "sfWeight", ("sf_sta_dn*sf_id*sf_iso*sf_trg", "sfWeightStaDn")),

    ReplaceWeight("SFIDUp", "sfWeight", ("sf_sta*sf_id_up*sf_iso*sf_trg", "sfWeightIDUp")),
    ReplaceWeight("SFIDDn", "sfWeight", ("sf_sta*sf_id_dn*sf_iso*sf_trg", "sfWeightIDDn")),

    ReplaceWeight("SFIsoUp", "sfWeight", ("sf_sta*sf_id*sf_iso_up*sf_trg", "sfWeightIsoUp")),
    ReplaceWeight("SFIsoDn", "sfWeight", ("sf_sta*sf_id*sf_iso_dn*sf_trg", "sfWeightIsoDn")),

    ReplaceWeight("SFTrgUp", "sfWeight", ("sf_sta*sf_id*sf_iso*sf_trg_up", "sfWeightTrgUp")),
    ReplaceWeight("SFTrgDn", "sfWeight", ("sf_sta*sf_id*sf_iso*sf_trg_dn", "sfWeightTrgDn")),
]

pdf_weight = [
    ReplaceWeight(f"LHEPdfWeight{i}", "LHEPdfWeight", (f"LHEPdfWeight{i}", f"LHEPdfWeight{i}")) for i in range(1, 101)
] + [
    ReplaceWeight("LHEPdfWeightAlphaSUp",   "LHEPdfWeight", ("LHEPdfWeight102", "LHEPdfWeightAlphaSUp")),
    ReplaceWeight("LHEPdfWeightAlphaSDown", "LHEPdfWeight", ("LHEPdfWeight101", "LHEPdfWeightAlphaSDown")),
] + [
    ReplaceWeight("LHEScaleWeightMUFUp",      "LHEScaleWeight", ("LHEScaleWeight4", "LHEScaleWeightMUFUp")),
    ReplaceWeight("LHEScaleWeightMUFDown",    "LHEScaleWeight", ("LHEScaleWeight3", "LHEScaleWeightMUFDown")),
    ReplaceWeight("LHEScaleWeightMURUp",      "LHEScaleWeight", ("LHEScaleWeight6", "LHEScaleWeightMURUp")),
    ReplaceWeight("LHEScaleWeightMURDown",    "LHEScaleWeight", ("LHEScaleWeight1", "LHEScaleWeightMURDown")),
    ReplaceWeight("LHEScaleWeightMUFMURUp",   "LHEScaleWeight", ("LHEScaleWeight7", "LHEScaleWeightMUFMURUp")),
    ReplaceWeight("LHEScaleWeightMUFMURDown", "LHEScaleWeight", ("LHEScaleWeight0", "LHEScaleWeightMUFMURDown")),
]

pu_weight = [
    ReplaceWeight("puWeightUp", "puweight", ("puweightUp", "puWeightUp")),
    ReplaceWeight("puWeightDn", "puweight", ("puweightDn", "puWeightDn")),
]

ps_weight = [
    ReplaceWeight("PSWeightISRUp", "PSWeight", ("PSWeight0", "PSWeightISRUp")),
    ReplaceWeight("PSWeightISRDown", "PSWeight", ("PSWeight2", "PSWeightISRDown")),
    ReplaceWeight("PSWeightFSRUp", "PSWeight", ("PSWeight1", "PSWeightFSRUp")),
    ReplaceWeight("PSWeightFSRDown", "PSWeight", ("PSWeight3", "PSWeightFSRDown")),
]
# mt_cuts = [
#     ReplaceCut("pfmtcut10", "pfmtcut", ("(sqrt(2.*pt_1_corr*pfmet_corr*(1.-cos(phi_1 - pfmetphi_corr)))) > 10.", "pfmtcut10")),
#     ReplaceCut("pfmtcut20", "pfmtcut", ("(sqrt(2.*pt_1_corr*pfmet_corr*(1.-cos(phi_1 - pfmetphi_corr)))) > 20.", "pfmtcut20")),
#     ReplaceCut("pfmtcut30", "pfmtcut", ("(sqrt(2.*pt_1_corr*pfmet_corr*(1.-cos(phi_1 - pfmetphi_corr)))) > 30.", "pfmtcut30")),
#     ReplaceCut("pfmtcut40", "pfmtcut", ("(sqrt(2.*pt_1_corr*pfmet_corr*(1.-cos(phi_1 - pfmetphi_corr)))) > 40.", "pfmtcut40")),
# ]
