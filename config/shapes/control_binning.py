from ntuple_processor import Histogram
import os 
import numpy as np
from config.common.variables import seperate_var, get_base_name

#corresponding lists are in the order [lower_bound, upper_bound, number_of_bins]
uniformBinningDict = {
    'pt': [25, 200, 50],
    'pt_corr': [25, 200, 50],
    'eta': [(-2.5), 2.5, 50],
    'etaSC': [(-2.5), 2.5, 50],
    'phi': [(-np.pi), (np.pi), 40],
    'sf_sel': [0, 1.2, 120],
    'sf_trg': [0, 1.2, 120],
    'm_vis': [60, 120, 60],
    'm_vis_corr': [60, 120, 60],
    'pt_vis': [0, 200, 40],
    'pt_vis_corr': [0, 200, 40],
    'bosonpt': [0, 1000, 2000],
    'bosonrap': [-5, 5, 50],
    'nStations': [0, 6, 6],
    'nTrackerLayers': [0, 20, 20],
    'deltaetaSC': [(-0.05), 0.05, 40],
    'eInvMinusPInv': [(-0.005), (0.005), 40],
    'hoe': [0, 0.05, 40],
    'scEtOverPt': [(-0.2), 0.05, 40],
    'sieie': [0, 0.04, 40],
    'lostHits': [0, 2, 2],
    'uP1_uncorrected': [-300, 150, 50],
    'uP1_corr': [-300, 150, 50],
    'uP2_uncorrected': [-150, 150, 40],
    'uP2_corr': [-150, 150, 40],
    'pfuP1_uncorrected': [-300, 150, 50],
    'pfuP1_corr': [-300, 150, 50],
    'pfuP2_uncorrected': [-150, 150, 40],
    'pfuP2_corr': [-150, 150, 40],
    'iso': [0, 2, 200],
    'mt': [0, 120, 20],
    #'mt': [0, 20, 10],
    'mt_uncorrected': [0, 120, 20],
    #'mt_uncorrected': [0, 20, 10],
    'mt_corr': [0, 120, 20],
    #'mt_corr': [0, 20, 10],
    'mt_lepuncorr': [0, 120, 20],
    #'mt_lepuncorr': [0, 20, 10],
    'met_uncorrected': [0, 120, 20],
    'met_corr': [0, 120, 20],
    'metphi_uncorrected': [(-np.pi), (np.pi), 40],
    'metphi_corr': [(-np.pi), (np.pi), 40],
    'pfmt': [0, 120, 20],
    #'pfmt': [0, 20, 10],
    'pfmt_uncorrected': [0, 120, 20],
    'pfmt_corr': [0, 120, 20],
    'pfmt_lepuncorr': [0, 120, 20],
    #'pfmt_uncorrected': [0, 20, 10],
    #'pfmt_corr': [0, 20, 10],
    #'pfmt_lepuncorr': [0, 20, 10],
    'pfmet_uncorrected': [0, 120, 20],
    'pfmet_corr': [0, 120, 20],
    'pfmetphi_uncorrected': [(-np.pi), (np.pi), 40],
    'pfmetphi_corr': [(-np.pi), (np.pi), 40],
    'm_toFit': [60, 120, 30],
    'm_toFit_lepuncorr': [60, 120, 30],

    'ptOverPfmt_corr': [0, 3, 30],
    'm_vis_corr_up': [60, 120, 60],
    'm_vis_corr_dn': [60, 120, 60],

}

zptbins = [0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,65,70,75,80,90,100,125,150,1000]
zrapbins = [0.0, 0.5, 1.0, 1.e9]

#function which allows you to add new variables not in the ntuple
def renamingVariable(var, base, channel):
    selection = "1"
    idx = "2" if ("_2" in var) else "1"

    # fine binned mass
    if "m_toFit" in var:
        base = "m_vis_corr"
    if "m_toFit_lepuncorr" in var:
        base = "m_vis"

    # mT
    if "mt_uncorrected" in var:
        base = "(sqrt(2.*pt_1*met_uncorrected*(1.-cos(phi_1 - metphi_uncorrected))))"

    if "mt_corr" in var:
        base = "(sqrt(2.*pt_1_corr*met_corr*(1.-cos(phi_1 - metphi_corr))))"
    
    if "mt_lepuncorr" in var:
        base = "(sqrt(2.*pt_1*met_corr*(1.-cos(phi_1 - metphi_corr))))"

    if "pfmt_uncorrected" in var:
        base = "(sqrt(2.*pt_1*pfmet_uncorrected*(1.-cos(phi_1 - pfmetphi_uncorrected))))"

    if "pfmt_corr" in var:
        base = "(sqrt(2.*pt_1_corr*pfmet_corr*(1.-cos(phi_1 - pfmetphi_corr))))"

    if "pfmt_lepuncorr" in var:
        base = "(sqrt(2.*pt_1*pfmet_corr*(1.-cos(phi_1 - pfmetphi_corr))))"

    if "pfmt_corr_double" in var:
        base = "(sqrt(2.*pt_1_corr*pfmet_corr_double*(1.-cos(phi_1 - pfmetphi_corr_double))))"
    
    if "pfmet_corr_double" in var:
        base = "pfmet_corr_double"

    if "pfmt_corr_sigOnly" in var:
        base = "(sqrt(2.*pt_1_corr*pfmet_corr_sigOnly*(1.-cos(phi_1 - pfmetphi_corr_sigOnly))))"
    
    if "pfmet_corr_sigOnly" in var:
        base = "pfmet_corr_sigOnly"

    if "pfmt_corr_zrap0" in var:
        base = "(sqrt(2.*pt_1_corr*pfmet_corr_zrap0*(1.-cos(phi_1 - pfmetphi_corr_zrap0))))"
    
    if "pfmet_corr_zrap0" in var:
        base = "pfmet_corr_zrap0"

    if "pfmt_corr_zrap1" in var:
        base = "(sqrt(2.*pt_1_corr*pfmet_corr_zrap1*(1.-cos(phi_1 - pfmetphi_corr_zrap1))))"
    
    if "pfmet_corr_zrap1" in var:
        base = "pfmet_corr_zrap1"

    if "pfmt_corr_zrap2" in var:
        base = "(sqrt(2.*pt_1_corr*pfmet_corr_zrap2*(1.-cos(phi_1 - pfmetphi_corr_zrap2))))"
    
    if "pfmet_corr_zrap2" in var:
        base = "pfmet_corr_zrap2"

    # ptOvermT
    if "ptOverMt_uncorrected" in var:
        base = "(pt_1/(sqrt(2.*pt_1*met_uncorrected*(1.-cos(phi_1 - metphi_uncorrected)))))"

    if "ptOverMt_corr" in var:
        base = "(pt_1_corr/(sqrt(2.*pt_1_corr*met_corr*(1.-cos(phi_1 - metphi_corr)))))"

    if "ptOverPfmt_corr" in var:
        base = "(pt_1_corr/(sqrt(2.*pt_1_corr*pfmet_corr*(1.-cos(phi_1 - pfmetphi_corr)))))"

    if "ptOverPfmt_corr_double" in var:
        base = "(pt_1_corr/(sqrt(2.*pt_1_corr*pfmet_corr_double*(1.-cos(phi_1 - pfmetphi_corr_double)))))"
    
    if "ptOverPfmt_corr_sigOnly" in var:
        base = "(pt_1_corr/(sqrt(2.*pt_1_corr*pfmet_corr_sigOnly*(1.-cos(phi_1 - pfmetphi_corr_sigOnly)))))"
    
    if "ptOverPfmt_corr_zrap0" in var:
        base = "(pt_1_corr/(sqrt(2.*pt_1_corr*pfmet_corr_zrap0*(1.-cos(phi_1 - pfmetphi_corr_zrap0)))))"

    if "ptOverPfmt_corr_zrap1" in var:
        base = "(pt_1_corr/(sqrt(2.*pt_1_corr*pfmet_corr_zrap1*(1.-cos(phi_1 - pfmetphi_corr_zrap1)))))"
    
    if "ptOverPfmt_corr_zrap2" in var:
        base = "(pt_1_corr/(sqrt(2.*pt_1_corr*pfmet_corr_zrap2*(1.-cos(phi_1 - pfmetphi_corr_zrap2)))))"

    for istat in range(15):
        found = False
        for updn in ["Up", "Down"]:
            if f"pfmt_corr_stat{istat}{updn}" in var:
                base = f"(sqrt(2.*pt_1_corr*pfmet_corr_stat{istat}{updn}*(1.-cos(phi_1 - pfmetphi_corr_stat{istat}{updn}))))"
                found = True
                break
            if f"pfmet_corr_stat{istat}{updn}" in var:
                base = f"pfmet_corr_stat{istat}{updn}"
                found = True
                break
            if f"ptOverPfmt_corr_stat{istat}{updn}" in var:
                base = f"(pt_1_corr/(sqrt(2.*pt_1_corr*pfmet_corr_stat{istat}{updn}*(1.-cos(phi_1 - pfmetphi_corr_stat{istat}{updn})))))"
                found = True
                break
        if found:
            break

    # charge
    if "_pos" in var:
        selection = selection + "*(q_1 > 0)"
    elif "_neg" in var:
        selection = selection + "*(q_1 < 0)"

    # isolation
    if "_isoSR" in var:
        if channel == "mmet":
            selection = selection + "*(iso_1 < 0.15)"
        elif channel == "emet":
            if "_barrel" in var:
                selection = selection + "*(iso_1 < 0.0478+0.506/pt_1)"
            elif "_endcap" in var:
                selection = selection + "*(iso_1 < 0.0658+0.963/pt_1)"
    elif "_iso5" in var:
        selection = selection + "*(iso_1 > 0.20 && iso_1 < 0.25)"
    elif "_iso6" in var:
        selection = selection + "*(iso_1 > 0.25 && iso_1 < 0.30)"
    elif "_iso7" in var:
        selection = selection + "*(iso_1 > 0.30 && iso_1 < 0.35)"
    elif "_iso8" in var:
        selection = selection + "*(iso_1 > 0.35 && iso_1 < 0.40)"
    elif "_iso9" in var:
        selection = selection + "*(iso_1 > 0.40 && iso_1 < 0.45)"
    elif "_iso10" in var:
        selection = selection + "*(iso_1 > 0.45 && iso_1 < 0.50)"
    elif "_iso11" in var:
        selection = selection + "*(iso_1 > 0.50 && iso_1 < 0.55)"
    elif "_iso12" in var:
        selection = selection + "*(iso_1 > 0.55 && iso_1 < 0.60)"
    elif "_iso13" in var:
        selection = selection + "*(iso_1 > 0.60 && iso_1 < 0.65)"
    elif "_iso14" in var:
        selection = selection + "*(iso_1 > 0.65 && iso_1 < 0.70)"
    elif "_iso15" in var:
        selection = selection + "*(iso_1 > 0.70 && iso_1 < 0.75)"
    elif "_iso16" in var:
        selection = selection + "*(iso_1 > 0.75 && iso_1 < 0.80)"
    elif "_iso17" in var:
        selection = selection + "*(iso_1 > 0.80 && iso_1 < 0.85)"
    elif "_iso18" in var:
        selection = selection + "*(iso_1 > 0.85 && iso_1 < 0.90)"
    elif "_iso19" in var:
        selection = selection + "*(iso_1 > 0.90 && iso_1 < 0.95)"
    elif "_iso20" in var:
        selection = selection + "*(iso_1 > 0.95 && iso_1 < 1.00)"

    # lepcorr bins
    if "_lepcorr" in var:
        # lepcorrbins = {'abseta': [0., 1.2, 2.4], 'pt': [25, 30, 35, 39, 42, 45, 49, 60, 120, 9999]}  # v06
        lepcorrbins = {'eta': [-2.4, -1.2, 0., 1.2, 2.4], 'pt': [25, 30, 35, 39, 42, 45, 49, 60, 120, 9999]}  # v05

        bin1, bin2 = list(lepcorrbins.keys())[0], list(lepcorrbins.keys())[1]
        nbins1, nbins2 = len(lepcorrbins[bin1])-1, len(lepcorrbins[bin2])-1

        the_idx_1 = None
        the_idx_2 = None
        for i in range(nbins1):
            for j in range(nbins2):
                lepcorr_str = "_lepcorr_{}{}_{}{}".format(bin1, i, bin2, j)
                if var.endswith(lepcorr_str):
                    the_idx_1 = i
                    the_idx_2 = j
                    break
            if the_idx_2 is not None:
                break
        if the_idx_1 is not None and the_idx_2 is not None:
            bin1_l, bin1_r = lepcorrbins[bin1][the_idx_1], lepcorrbins[bin1][the_idx_1+1]
            bin2_l, bin2_r = lepcorrbins[bin2][the_idx_2], lepcorrbins[bin2][the_idx_2+1]

            filter_template = "({bin}_{n}{corr} > {bin_l} && {bin}_{n}{corr} < {binr})"

            corr_str = ""
            if "_corr" in var:
                corr_str = "_corr"

            filter1a = filter_template.format(bin=bin1, n=1, bin_l=bin1_l, binr = bin1_r, corr="")
            filter1b = filter_template.format(bin=bin2, n=1, bin_l=bin2_l, binr = bin2_r, corr=corr_str)
            
            filter2a = filter_template.format(bin=bin1, n=2, bin_l=bin1_l, binr = bin1_r, corr="")
            filter2b = filter_template.format(bin=bin2, n=2, bin_l=bin2_l, binr = bin2_r, corr=corr_str)

            # filters = "({} && {}) && ({} && {})".format(filter1a, filter1b, filter2a, filter2b)  # v06
            filters = "({} && {}) || ({} && {})".format(filter1a, filter1b, filter2a, filter2b)  # v05
            filters = filters.replace("abseta_1", "abs(eta_1)").replace("abseta_2", "abs(eta_2)")
            selection = selection + "*({})".format(filters)

    # zpt bins
    if "_zpt" in var:
        if 'uP' in var or "met" in var or "mt" in var or "bosonpt" in var or "bosonrap" in var:
            # Z pt bins
            the_zpt_str = None
            the_zpt_idx = None
            for izptbin, zptbin in enumerate(zptbins):
                if izptbin == len(zptbins)-1:
                    continue
                zptbin_str = "_zpt{}".format(izptbin)
                if var.split("_zrap")[0].endswith(zptbin_str):
                    the_zpt_str = zptbin_str
                    the_zpt_idx = izptbin
                    break
            if the_zpt_str is not None:
                selection = selection + "*(bosonpt > {} && bosonpt <= {})".format(
                    zptbins[the_zpt_idx],
                    zptbins[the_zpt_idx+1]
                )

            # Z rapidity bins
            the_zrap_str = None
            the_zrap_idx = None
            for izrapbin, zrapbin in enumerate(zrapbins):
                if izrapbin == len(zrapbins)-1:
                    continue
                zrapbin_str = "_zrap{}".format(izrapbin)
                if var.endswith(zrapbin_str):
                    the_zrap_str = zrapbin_str
                    the_zrap_idx = izrapbin
                    break
            if the_zrap_str is not None:
                selection = selection + "*(abs(bosonrap) > {} && abs(bosonrap) <= {})".format(
                    zrapbins[the_zrap_idx],
                    zrapbins[the_zrap_idx+1]
                )

    # barrel and endcap
    if channel == "ee" or channel == "emet":
        #barrel or endcap
        if "_barrel" in var:
            selection = selection + "*(abs(eta_{idx} + deltaetaSC_{idx}) <= 1.479)".format(idx=idx)
        elif "_endcap" in var:
            selection = selection + "*(abs(eta_{idx} + deltaetaSC_{idx}) > 1.479)".format(idx=idx)
    elif channel == "mm" or channel == "mmet":
        if "_barrel" in var:
            selection = selection + "*(abs(eta_{idx}) <= 1.2)".format(idx=idx)
        elif "_endcap" in var:
            selection = selection + "*(abs(eta_{idx}) > 1.2)".format(idx=idx)

    if channel == "ee":
        if "_BB" in var:
            selection = selection + "*(abs(eta_1 + deltaetaSC_1) <= 1.479 && abs(eta_2 + deltaetaSC_2) <= 1.479)"
        elif "_BE" in var:
            selection = selection + "*((abs(eta_1 + deltaetaSC_1) <= 1.479 && abs(eta_2 + deltaetaSC_2) > 1.479) || (abs(eta_1 + deltaetaSC_1) > 1.479 && abs(eta_2 + deltaetaSC_2) <= 1.479))"
        elif "_EE" in var:
            selection = selection + "*(abs(eta_1 + deltaetaSC_1) > 1.479 && abs(eta_2 + deltaetaSC_2) > 1.479)"
    elif channel == "mm":
        if "_BB" in var:
            selection = selection + "*(abs(eta_1) <= 1.2 && abs(eta_2) <= 1.2)"
        elif "_BE" in var:
            selection = selection + "*((abs(eta_1) <= 1.2 && abs(eta_2) > 1.2) || (abs(eta_1) > 1.2 && abs(eta_2) < 1.2))"
        elif "_EE" in var:
            selection = selection + "*(abs(eta_1) > 1.2 && abs(eta_2) > 1.2)"

    full_var = base
    if selection != "1":
        full_var = "{base}*(int(bool({sel}))) -1e9*(int(!(bool({sel}))))".format(base=base, sel=selection)

    return full_var

def get_control_binning(channel, variable_list):
    print("")
    print("="*50)
    print(f"Channel: {channel}")
    

    control_binning = {channel:dict()}
    for i in range(len(variable_list)):
        var_base = get_base_name(variable_list[i])
        var_name = renamingVariable(variable_list[i], var_base, channel)
        print(f"\t{variable_list[i]}: {var_name}")

        # makes the bins and bounds
        histo_bin_list=list()
        if 'run' in variable_list[i]:
            histo_bin_list = [355862,355863,355870,355871,355872,355892,355912,355913,355921,355933,355942,355988,355989,356043,356071,356074,356075,356076,356077,356135,356309,356316,356321,356322,356323,356371,356375,356378,356381,356383,356385,356386,356426,356428,356433,356434,356435,356446,356523,356531,356563,356568,356569,356570,356576,356578,356580,356582,356614,356615,356619,356811,356812,356813,356814,356815,356824,356908,356919,356937,356946,356947,356948,356949,356951,356954,356955,356956,356968,356969,356970,356998,356999,357000,357001,357079,357080,357081,357101,357102,357104,357106,357112,357268,357271,357328,357329,357330,357331,357332,357333,357401,357406,357438,357440,357441,357442,357447,357472,357478,357479,357482,357538,357542,357550,357610,357611,357612,357613,357688,357696,357697,357698,357699,357700,357720,357732,357734,357735,357754,357756,357757,357758,357759,357766,357777,357778,357779,357781,357802,357803,357804,357805,357806,357807,357808,357809,357812,357813,357814,357815,357898,357899,357900,357901]
        elif "_sfbin" in variable_list[i]:
            histo_bin_list = [-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4]
        elif "uP" in variable_list[i] and "_zpt" in variable_list[i]:
            zptbinidx = variable_list[i].split("_zpt")[-1]
            if "_zrap" in zptbinidx:
                zptbinidx = zptbinidx.split("_zrap")[0]
            zptbinidx = int(zptbinidx)
            window = 250.
            if zptbins[zptbinidx] > 140:
                window = 500.
            if 'uP1_' in variable_list[i]:
                histo_bin_list = [-window - zptbins[zptbinidx] + float(x)*(window/100.) for x in range(0, 201)]
            elif 'uP2_' in variable_list[i]:
                histo_bin_list = [-window + float(x)*(window/100.) for x in range(0, 201)]
        else:
            binning_base = var_base.replace("_1", "").replace("_2", "")
            for k in range(0, uniformBinningDict[binning_base][2]+1):
                var_bin_list=uniformBinningDict[binning_base][0]+k*(uniformBinningDict[binning_base][1] - uniformBinningDict[binning_base][0])/(uniformBinningDict[binning_base][2])
                histo_bin_list.append(var_bin_list)

        # control binning
        control_binning[channel][variable_list[i]]=Histogram(variable_list[i], var_name, histo_bin_list)
    
    print("="*50)
    print("")

    return control_binning
