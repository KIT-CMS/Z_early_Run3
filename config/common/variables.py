def seperate_var(input_list, doZpt=False, doQCD=False, doLepCorrBins=False):
    variable_list = list()

    for var in input_list:
        variable_list.append(var)
        if "run" in var:
            continue

        # lepcorr bins
        if doLepCorrBins:
            # lepcorrbins = {'abseta': [0., 1.2, 2.4], 'pt': [25, 30, 35, 39, 42, 45, 49, 60, 120, 9999]}  # v06
            lepcorrbins = {'eta': [-2.4, -1.2, 0., 1.2, 2.4], 'pt': [25, 30, 35, 39, 42, 45, 49, 60, 120, 9999]}  # v05

            bin1, bin2 = list(lepcorrbins.keys())[0], list(lepcorrbins.keys())[1]
            nbins1, nbins2 = len(lepcorrbins[bin1])-1, len(lepcorrbins[bin2])-1
            for i in range(nbins1):
                for j in range(nbins2):
                    lepcorr_str = "_lepcorr_{}{}_{}{}".format(bin1, i, bin2, j)
                    variable_list.append(var+lepcorr_str)
        # zpt bins
        elif doZpt:
            if 'uP' in var or "met" in var or "mt" in var or "bosonpt" in var or "bosonrap" in var:
                zptbins = [0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,65,70,75,80,90,100,125,150,1000]
                zrapbins = [0.0, 0.5, 1.0, 1.e9]
                for izptbin, zptbin in enumerate(zptbins):
                    if izptbin == len(zptbins)-1:
                        continue
                    zptbin_str = "_zpt{}".format(izptbin)
                    variable_list.append(var+zptbin_str)
                    for izrapbin, zrapbin in enumerate(zrapbins):
                        if izrapbin == len(zrapbins)-1:
                            continue
                        zrapbin_str = "_zrap{}".format(izrapbin)
                        variable_list.append(var+zptbin_str+zrapbin_str)
        elif doQCD:
            charge_list = ["_pos", "_neg"]
            iso_list = ["_isoSR", "_iso5", "_iso6", "_iso7", "_iso8", "_iso9", "_iso10", "_iso11", "_iso12", "_iso13", "_iso14", "_iso15", "_iso16", "_iso17", "_iso18", "_iso19", "_iso20"]
            for charge in charge_list:
                for iso in iso_list:
                    if var+charge+iso not in variable_list:
                        variable_list.append(var+charge+iso)
                    # if var+charge+eta+iso not in variable_list:
                    #     variable_list.append(var+charge+eta+iso)
        else:
            charge_list = ["", "_pos", "_neg"]
            eta_list = ["", "_barrel", "_endcap"]
            for charge in charge_list:
                if var+charge not in variable_list:
                    variable_list.append(var+charge)
                for eta in eta_list:
                    if 'eta_' in var:
                        continue
                    if 'm_vis' in var or 'pt_vis' in var:
                        continue
                    if '_BB' in var or '_BE' in var or '_EE' in var:
                        continue
                    if "met" in var or "mt" in var or "ptOverMt" in var or "uP" in var:
                        continue
                    if var+eta not in variable_list:
                        variable_list.append(var+eta)
                    if var+charge+eta not in variable_list:
                        variable_list.append(var+charge+eta)

    if doQCD:
        variable_list_QCD = []
        for _var in variable_list:
            if "_iso" in _var:
                variable_list_QCD.append(_var)
        return variable_list_QCD

    return variable_list

def get_base_name(var):
    var_base = var.replace("_barrel", "").replace("_endcap", "")\
                  .replace("_BB", "").replace("_BE", "").replace("_EE", "")\
                  .replace("_pos", "").replace("_neg", "")\
                  .replace("_sfbin", "")\
                  .replace("_double", "")\
                  .replace("_sigOnly", "")

    var_base = var_base.split("_stat")[0]
    var_base = var_base.split("_iso")[0]
    var_base = var_base.split("_zpt")[0]
    var_base = var_base.split("_zrap")[0]
    var_base = var_base.split("_lepcorr")[0]

    return var_base

def get_all_variables(doZpt, doQCD, doLepCorrBins):
    if doLepCorrBins:
        variable_dict = {
            "mm": ["pt_1", "pt_1_corr", "eta_1", "m_vis", "m_vis_corr"],
            "mmet": [],
            "ee": [],
            "emet": [],
        }
        return variable_dict

    if doZpt:
        variable_dict = {
            "mm": [
                ## "uP1_uncorrected", "uP2_uncorrected", "bosonpt", "bosonrap"
                # "bosonpt", "bosonrap",
                # "pfuP1_uncorrected", "pfuP2_uncorrected",

                ## "pfmet_uncorrected", "pfmt_uncorrected", "pfuP1_uncorrected", "pfuP2_uncorrected",
                ## "pfmet_corr", "pfmt_corr", "pfuP1_corr", "pfuP2_corr",
            ],
            "mmet": [
                ## "uP1_uncorrected_pos", "uP2_uncorrected_pos", "bosonpt_pos", "bosonrap_pos",
                ## "uP1_uncorrected_neg", "uP2_uncorrected_neg", "bosonpt_neg", "bosonrap_neg",
                
                #"bosonpt_pos", "bosonrap_pos",
                "bosonpt_neg", "bosonrap_neg",
                
                #"pfuP1_uncorrected_pos", "pfuP2_uncorrected_pos",
                #"pfuP1_uncorrected_neg", "pfuP2_uncorrected_neg",

                ## "pfmet_uncorrected_pos", "pfuP1_uncorrected_pos", "pfuP2_uncorrected_pos",
                ## "pfmet_corr_pos", "pfmt_corr_pos", "pfuP1_corr_pos", "pfuP2_corr_pos",

                ## "pfmet_uncorrected_neg", "pfuP1_uncorrected_neg", "pfuP2_uncorrected_neg",
                ## "pfmet_corr_neg", "pfmt_corr_neg", "pfuP1_corr_neg", "pfuP2_corr_neg",
            ],
            "ee": [],
            "emet": [],
        }
        return variable_dict

    if doQCD:
        variable_dict = {
            "mm": [],
            "mmet": [
                # "eta_1", "iso_1",
                # "met_corr", "mt_corr",

                "pfmet_corr",
                "pfmt_corr",
                # "ptOverPfmt_corr"
            ],
            "ee": [],
            "emet": [],
        }
        return variable_dict

    common_vars_mm = [
        # "run",

        "pt_1", "eta_1", "phi_1",
        "pt_1_corr",

        "pt_2", "eta_2", "phi_2",
        "pt_2_corr",

        # "iso_1",

        "m_vis", "m_vis_BB", "m_vis_BE", "m_vis_EE",
        "m_vis_corr", "m_vis_corr_BB", "m_vis_corr_BE", "m_vis_corr_EE",
        "m_vis_corr_up", "m_vis_corr_dn",

        "pt_vis",
        "pt_vis_corr",

        "eta_1_sfbin",
        "eta_2_sfbin",

        "m_toFit"
    ]

    common_vars_mmet = [
        "pt_1", "eta_1", "phi_1",
        "pt_1_corr",
    ]

    met_vars = [
        # "met_uncorrected", "metphi_uncorrected", "mt_uncorrected",
        # "met_corr", "metphi_corr", "mt_corr",

        # "uP1_uncorrected", "uP2_uncorrected",
        # "uP1_corr", "uP2_corr",

        "pfmet_uncorrected", "pfmetphi_uncorrected", "pfmt_uncorrected",
        "pfmet_corr", "pfmetphi_corr",  "pfmt_corr",

        "pfuP1_uncorrected", "pfuP2_uncorrected",
        "pfuP1_corr", "pfuP2_corr",
    ]

    final_vars_mm = [
        "m_toFit",
        # "m_vis",
        "m_toFit_up",
        "m_toFit_dn"
        ]
    final_vars_mmet = [
        #"pfmet_corr_pos", "pfmet_corr_neg",
        "pfmt_corr_pos", "pfmt_corr_neg",
    ]

    final_vars_mmet_syst = [
        #"pfmet_corr_double_neg", "pfmet_corr_double_pos",
        "pfmt_corr_double_neg",  "pfmt_corr_double_pos",

        #"pfmet_corr_sigOnly_neg", "pfmet_corr_sigOnly_pos",
        "pfmt_corr_sigOnly_neg",  "pfmt_corr_sigOnly_pos",

        #"pfmet_corr_zrap0_neg", "pfmet_corr_zrap0_pos",
        "pfmt_corr_zrap0_neg",  "pfmt_corr_zrap0_pos",

        #"pfmet_corr_zrap1_neg", "pfmet_corr_zrap1_pos",
        "pfmt_corr_zrap1_neg",  "pfmt_corr_zrap1_pos",

        #"pfmet_corr_zrap2_neg", "pfmet_corr_zrap2_pos",
        "pfmt_corr_zrap2_neg",  "pfmt_corr_zrap2_pos",
    ] + [
        "pfmt_corr_stat{istat}{updn}_{posneg}".format(istat=istat, updn=updn, posneg=posneg) for istat in range(15) for updn in ["Up", "Down"] for posneg in ["pos", "neg"]
    ] + [
      #  "pfmet_corr_stat{istat}{updn}_{posneg}".format(istat=istat, updn=updn, posneg=posneg) for istat in range(15) for updn in ["Up", "Down"] for posneg in ["pos", "neg"]
    ]

    variable_dict = {
        "mm": final_vars_mm,
        # "mmet": final_vars_mmet,
        "mmet": final_vars_mmet_syst,

        # "mm": common_vars_mm + met_vars,
        # "mmet": common_vars_mmet + met_vars,

        # "mm": common_vars_mm + met_vars,  # + ["nTrackerLayers_1", "nStations_1"],
        # "mmet": common_vars_mmet + met_vars,  # + ["nTrackerLayers_1", "nStations_1"],
        # "ee": common_vars + met_vars,  # + ["eInvMinusPInv_1", "scEtOverPt_1", "sieie_1", "hoe_1"],
        # "emet": common_vars + met_vars,  # + ["eInvMinusPInv_1", "scEtOverPt_1", "sieie_1", "hoe_1"],
    }

    return variable_dict
