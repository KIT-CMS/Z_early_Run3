#!/usr/bin/env python
import argparse
import logging

from ntuple_processor import Histogram
from ntuple_processor import dataset_from_artusoutput, Unit, UnitManager, GraphManager, RunManager

from config.shapes.channel_selection import channel_selection
from config.shapes.file_names import files
from config.shapes.process_selection import DY_process_selection, TT_process_selection, VV_process_selection, W_process_selection, ZTT_process_selection, ZL_process_selection, ZJ_process_selection, TTT_process_selection, TTL_process_selection, TTJ_process_selection, VVT_process_selection, VVJ_process_selection, VVL_process_selection, ggH125_process_selection, qqH125_process_selection, ZTT_embedded_process_selection, ZH_process_selection, WH_process_selection, ggHWW_process_selection, qqHWW_process_selection, ZHWW_process_selection, WHWW_process_selection, ttH_process_selection
from config.shapes.process_selection import SUSYbbH_process_selection, SUSYggH_process_selection, SUSYggH_Ai_contribution_selection, SUSYggH_At_contribution_selection, SUSYggH_Ab_contribution_selection, SUSYggH_Hi_contribution_selection, SUSYggH_Ht_contribution_selection, SUSYggH_Hb_contribution_selection, SUSYggH_hi_contribution_selection, SUSYggH_ht_contribution_selection, SUSYggH_hb_contribution_selection
from config.shapes.category_selection import categorization
# Variations for estimation of fake processes
from config.shapes.variations import same_sign, same_sign_em, anti_iso_lt, anti_iso_tt
# Energy scale uncertainties
from config.shapes.variations import tau_es_3prong, tau_es_3prong1pizero, tau_es_1prong, tau_es_1prong1pizero, emb_tau_es_3prong, emb_tau_es_3prong1pizero, emb_tau_es_1prong, emb_tau_es_1prong1pizero, jet_es, mu_fake_es_1prong, mu_fake_es_1prong1pizero, ele_es, ele_res, emb_e_es, ele_fake_es_1prong, ele_fake_es_1prong1pizero
# MET related uncertainties.
from config.shapes.variations import met_unclustered, recoil_resolution, recoil_response
# efficiency uncertainties
from config.shapes.variations import tau_id_eff_lt, tau_id_eff_tt, emb_tau_id_eff_lt, emb_tau_id_eff_tt
# fake rate uncertainties
from config.shapes.variations import jet_to_tau_fake, zll_et_fake_rate_2016, zll_et_fake_rate_2017, zll_et_fake_rate_2018, zll_mt_fake_rate_2016, zll_mt_fake_rate_2017, zll_mt_fake_rate_2018
# trigger efficiencies
from config.shapes.variations import tau_trigger_eff_tt, tau_trigger_eff_emb_tt, lep_trigger_eff_mt_2016, lep_trigger_eff_et_2016, lep_trigger_eff_et_emb_2016, lep_trigger_eff_mt_emb_2016, tau_trigger_eff_et_2016, tau_trigger_eff_mt_2016, tau_trigger_eff_et_emb_2016, tau_trigger_eff_mt_emb_2016, lep_trigger_eff_et_2017, lep_trigger_eff_mt_2017, lep_trigger_eff_et_emb_2017, lep_trigger_eff_mt_emb_2017, tau_trigger_eff_et_2017, tau_trigger_eff_mt_2017, tau_trigger_eff_et_emb_2017, tau_trigger_eff_mt_emb_2017, lep_trigger_eff_mt_2018, lep_trigger_eff_et_2018, lep_trigger_eff_et_emb_2018, lep_trigger_eff_mt_emb_2018, tau_trigger_eff_et_2018, tau_trigger_eff_mt_2018, tau_trigger_eff_et_emb_2018, tau_trigger_eff_mt_emb_2018
from config.shapes.variations import prefiring, btag_eff, mistag_eff, ggh_acceptance, qqh_acceptance, zpt, top_pt, emb_decay_mode_eff
from config.shapes.variations import ff_variations_lt, ff_variations_tt, qcd_variations_em
from config.shapes.control_binning import control_binning, minimal_control_plot_set

logger = logging.getLogger("")


def setup_logging(output_file, level=logging.DEBUG):
    logger.setLevel(level)
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    file_handler = logging.FileHandler(output_file, "w")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)


def parse_arguments():
    parser = argparse.ArgumentParser(
            description="Produce shapes for the legacy MSSM analysis.")
    parser.add_argument(
        "--era",
        required=True,
        type=str,
        help="Experiment era."
    )
    parser.add_argument(
        "--channels",
        default=[],
        type=lambda channellist: [channel for channel in channellist.split(',')],
        help="Channels to be considered, seperated by a comma without space"
    )
    parser.add_argument(
        "--directory",
        required=True,
        type=str,
        help="Directory with Artus outputs."
    )
    parser.add_argument(
        "--et-friend-directory",
        type=str,
        default=[],
        nargs='+',
        help=
        "Directories arranged as Artus output and containing a friend tree for et."
    )
    parser.add_argument(
        "--mt-friend-directory",
        type=str,
        default=[],
        nargs='+',
        help=
        "Directories arranged as Artus output and containing a friend tree for mt."
    )
    parser.add_argument(
        "--tt-friend-directory",
        type=str,
        default=[],
        nargs='+',
        help=
        "Directories arranged as Artus output and containing a friend tree for tt."
    )
    parser.add_argument(
        "--em-friend-directory",
        type=str,
        default=[],
        nargs='+',
        help=
        "Directories arranged as Artus output and containing a friend tree for em."
    )
    parser.add_argument(
        "--optimization-level",
        default=2,
        type=int,
        help="Level of optimization for graph merging."
    )
    parser.add_argument(
        "--num-processes",
        default=1,
        type=int,
        help="Number of processes to be used."
    )
    parser.add_argument(
        "--num-threads",
        default=1,
        type=int,
        help="Number of threads to be used."
    )
    parser.add_argument(
        "--skip-systematic-variations",
        action="store_true",
        help="Do not produce the systematic variations."
    )
    parser.add_argument(
        "--output-file",
        required=True,
        type=str,
        help="ROOT file where shapes will be stored."
    )
    parser.add_argument(
        "--control-plots",
        action="store_true",
        help="Produce shapes for control plots. Default is production of analysis shapes."
    )
    parser.add_argument(
        "--control-plot-set",
        default=minimal_control_plot_set,
        type=lambda varlist: [variable for variable in varlist.split(',')],
        help="Variables the shapes should be produced for."
    )
    return parser.parse_args()


def main(args):
    # Parse given arguments.
    friend_directories = {
        "et": args.et_friend_directory,
        "mt": args.mt_friend_directory,
        "tt": args.tt_friend_directory,
        "em": args.em_friend_directory,
    }
    if ".root" in args.output_file:
        output_file = args.output_file
        log_file = args.output_file.replace(".root", ".log")
    else:
        output_file = "{}.root".format(args.output_file)
        log_file = "{}.log".format(args.output_file)

    m_sv_hist = Histogram("m_sv_puppi", "m_sv_puppi", [i for i in range(0, 255, 5)])
    mt_tot_hist = Histogram("mt_tot_puppi", "mt_tot_puppi", [i for i in range(0, 3200, 10)])
    nominals = {}
    nominals[args.era] = {}
    nominals[args.era]['datasets'] = {}
    nominals[args.era]['units'] = {}

    susy_masses = {
        "2016": {
            "bbH": [ 80, 90, 100, 110, 120, 130, 140, 160, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1500, 1600, 1800, 2000, 2300, 2600, 2900, 3200],
            "bbH_nlo": [ 80, 90, 110, 120, 130, 140, 160, 180, 200, 250, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2300, 2600, 2900, 3200],
            "ggH": [ 80, 90, 100, 110, 120, 130, 140, 160, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1500, 1600, 1800, 2000, 2300, 2600, 2900, 3200],
        },
        "2017": {
            "bbH": [  80,   90,  100,  110,  120,  130,  140,  160,  180,  200, 250,  300,  350,  400,  600,  700,  800,  900, 1200, 1400, 1500, 1600, 1800, 2000, 2300, 2600, 2900, 3200],
            "bbH_nlo": [  80,   90,  110,  120,  125,  130,  140,  160,  180,  200, 250,  300,  350,  400,  500,  600,  700,  800,  900, 1000, 1200, 1400, 1600, 1800, 2000, 2300, 2600, 2900, 3200],
            "ggH": [  80,   90,  100,  110,  120,  130,  140,  180,  200, 250,  300,  350,  400,  450,  600,  700,  800,  900, 1200, 1400, 1500, 1600, 1800, 2000, 2300, 2600, 2900, 3200],
        },
        "2018": {
            "bbH": [  80,   90,  100,  110,  120,  130,  140,  160,  180,  200, 250,  300,  350,  400,  450,  600,  700,  800,  900, 1200, 1400, 1500, 1600, 1800, 2000, 2300, 2600, 2900, 3200],
            "bbH_nlo": [  80,   90,  100,  110,  120,  125,  130,  140,  160,  180,  200, 250,  300,  350,  400,  450,  500,  600,  700,  800,  900, 1000, 1200, 1400, 1600, 1800, 2000, 2300, 2600, 2900, 3200, 3500],
            "ggH": [  80,   90,  100,  110,  120,  130,  140,  160,  180,  200, 250,  300,  350,  400,  450,  600,  700,  800,  900, 1200, 1400, 1500, 1600, 1800, 2000, 2300, 2600, 2900, 3200],
        },
    }

    def get_nominal_datasets(era, channel):
        datasets = dict()
        for key, names in files[era][channel].items():
            datasets[key] = dataset_from_artusoutput(
                    key, names, channel + '_nominal', args.directory, friend_directories[channel])
        return datasets

    def get_analysis_units(channel, era, datasets):
        return {
                "data" : [Unit(
                            datasets["data"], [
                                channel_selection(channel, era),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "emb": [Unit(
                            datasets["EMB"], [
                                channel_selection(channel, era),
                                ZTT_embedded_process_selection(channel, era),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "ztt" : [Unit(
                            datasets["DY"], [
                                channel_selection(channel, era),
                                DY_process_selection(channel, era),
                                ZTT_process_selection(channel),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "zl" :  [Unit(
                            datasets["DY"], [
                                channel_selection(channel, era),
                                DY_process_selection(channel, era),
                                ZL_process_selection(channel),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "zj" :  [Unit(
                            datasets["DY"], [
                                channel_selection(channel, era),
                                DY_process_selection(channel, era),
                                ZJ_process_selection(channel),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "ttt" : [Unit(
                            datasets["TT"], [
                                channel_selection(channel, era),
                                TT_process_selection(channel, era),
                                TTT_process_selection(channel),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "ttl" : [Unit(
                            datasets["TT"], [
                                channel_selection(channel, era),
                                TT_process_selection(channel, era),
                                TTL_process_selection(channel),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "ttj" : [Unit(
                            datasets["TT"], [
                                channel_selection(channel, era),
                                TT_process_selection(channel, era),
                                TTJ_process_selection(channel),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "vvt" : [Unit(
                            datasets["VV"], [
                                channel_selection(channel, era),
                                VV_process_selection(channel, era),
                                VVT_process_selection(channel),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "vvl" : [Unit(
                            datasets["VV"], [
                                channel_selection(channel, era),
                                VV_process_selection(channel, era),
                                VVL_process_selection(channel),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "vvj" : [Unit(
                            datasets["VV"], [
                                channel_selection(channel, era),
                                VV_process_selection(channel, era),
                                VVJ_process_selection(channel),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "w"   : [Unit(
                            datasets["W"], [
                                channel_selection(channel, era),
                                W_process_selection(channel, era),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "ggh" : [Unit(
                            datasets["ggH"], [
                                channel_selection(channel, era),
                                ggH125_process_selection(channel, era),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "qqh" : [Unit(
                            datasets["qqH"], [
                                channel_selection(channel, era),
                                qqH125_process_selection(channel, era),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "wh"  : [Unit(
                            datasets["WH"], [
                                channel_selection(channel, era),
                                WH_process_selection(channel, era),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "zh"  : [Unit(
                            datasets["ZH"], [
                                channel_selection(channel, era),
                                ZH_process_selection(channel, era),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "tth"  : [Unit(
                            datasets["ttH"], [
                                channel_selection(channel, era),
                                ttH_process_selection(channel, era),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "gghww"  : [Unit(
                            datasets["ggHWW"], [
                                channel_selection(channel, era),
                                ggHWW_process_selection(channel, era),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "qqhww"  : [Unit(
                            datasets["qqHWW"], [
                                channel_selection(channel, era),
                                qqHWW_process_selection(channel, era),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "zhww"  : [Unit(
                            datasets["ZHWW"], [
                                channel_selection(channel, era),
                                ZHWW_process_selection(channel, era),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                "whww"  : [Unit(
                            datasets["WHWW"], [
                                channel_selection(channel, era),
                                WHWW_process_selection(channel, era),
                                category_selection], [m_sv_hist,mt_tot_hist]) for category_selection in categorization[channel]],
                **{"ggh{}".format(mass): [Unit(
                                            datasets["susyggH_{}".format(mass)], [
                                                channel_selection(channel, era),
                                                SUSYggH_process_selection(channel, era),
                                                contribution_selection(channel),
                                                category_selection], [m_sv_hist, mt_tot_hist]) for category_selection in categorization[channel]
                                                                                               for contribution_selection in [
                                                                                                                              SUSYggH_Ai_contribution_selection,
                                                                                                                              SUSYggH_At_contribution_selection,
                                                                                                                              SUSYggH_Ab_contribution_selection,
                                                                                                                              SUSYggH_Hi_contribution_selection,
                                                                                                                              SUSYggH_Ht_contribution_selection,
                                                                                                                              SUSYggH_Hb_contribution_selection,
                                                                                                                              SUSYggH_hi_contribution_selection,
                                                                                                                              SUSYggH_ht_contribution_selection,
                                                                                                                              SUSYggH_hb_contribution_selection]]
                                            for mass in susy_masses[era]["ggH"]},
                **{"bbh{}".format(mass): [Unit(
                                            datasets["susybbH_{}".format(mass)], [
                                                channel_selection(channel, era),
                                                SUSYggH_process_selection(channel, era),
                                                contribution_selection(channel),
                                                category_selection], [m_sv_hist, mt_tot_hist]) for category_selection in categorization[channel]]
                                            for mass in susy_masses[era]["bbH"]},
                **{"bbh{}_nlo".format(mass): [Unit(
                                                datasets["susybbH_nlo_{}".format(mass)], [
                                                    channel_selection(channel, era),
                                                    SUSYggH_process_selection(channel, era),
                                                    contribution_selection(channel),
                                                    category_selection], [m_sv_hist, mt_tot_hist]) for category_selection in categorization[channel]]
                                            for mass in susy_masses[era]["bbH_nlo"]},
        }

    def get_control_units(channel, era, datasets):
        return {
               'data' : [Unit(
                   datasets['data'],[
                       channel_selection(channel, era)],
                       [control_binning[channel][v] for v in set(control_binning[channel].keys()) & set(args.control_plot_set)])],
               'emb' : [Unit(
                   datasets['EMB'],[
                       channel_selection(channel, era),
                       ZTT_embedded_process_selection(channel, era)],
                       [control_binning[channel][v] for v in set(control_binning[channel].keys()) & set(args.control_plot_set)])],
               'ztt' : [Unit(
                   datasets['DY'], [
                       channel_selection(channel, era),
                       DY_process_selection(channel, era),
                       ZTT_process_selection(channel)],
                       [control_binning[channel][v] for v in set(control_binning[channel].keys()) & set(args.control_plot_set)])],
               'zl' : [Unit(
                   datasets['DY'], [
                      channel_selection(channel, era),
                      DY_process_selection(channel, era),
                      ZL_process_selection(channel)],
                      [control_binning[channel][v] for v in set(control_binning[channel].keys()) & set(args.control_plot_set)])],
               'zj' : [Unit(
                   datasets['DY'], [
                       channel_selection(channel, era),
                       DY_process_selection(channel, era),
                       ZJ_process_selection(channel)],
                       [control_binning[channel][v] for v in set(control_binning[channel].keys()) & set(args.control_plot_set)])],
               'ttl' : [Unit(
                   datasets['TT'], [
                       channel_selection(channel, era),
                       TT_process_selection(channel, era),
                       TTL_process_selection(channel)],
                       [control_binning[channel][v] for v in set(control_binning[channel].keys()) & set(args.control_plot_set)])],
               'ttt' : [Unit(
                   datasets['TT'], [
                       channel_selection(channel, era),
                       TT_process_selection(channel, era),
                       TTT_process_selection(channel)],
                       [control_binning[channel][v] for v in set(control_binning[channel].keys()) & set(args.control_plot_set)])],
               'ttj' : [Unit(
                   datasets['TT'], [
                       channel_selection(channel, era),
                       TT_process_selection(channel, era),
                       TTJ_process_selection(channel)],
                       [control_binning[channel][v] for v in set(control_binning[channel].keys()) & set(args.control_plot_set)])],
               'vvl' : [Unit(
                   datasets['VV'], [
                       channel_selection(channel, era),
                       VV_process_selection(channel, era),
                       VVL_process_selection(channel)],
                       [control_binning[channel][v] for v in set(control_binning[channel].keys()) & set(args.control_plot_set)])],
               'vvt' : [Unit(
                   datasets['VV'], [
                       channel_selection(channel, era),
                       VV_process_selection(channel, era),
                       VVT_process_selection(channel)],
                       [control_binning[channel][v] for v in set(control_binning[channel].keys()) & set(args.control_plot_set)])],
               'vvj' : [Unit(
                   datasets['VV'], [
                       channel_selection(channel, era),
                       VV_process_selection(channel, era),
                       VVJ_process_selection(channel)],
                       [control_binning[channel][v] for v in set(control_binning[channel].keys()) & set(args.control_plot_set)])],
               'w' :   [Unit(
                   datasets['W'], [
                       channel_selection(channel, era),
                       W_process_selection(channel, era)],
                       [control_binning[channel][v] for v in set(control_binning[channel].keys()) & set(args.control_plot_set)])],
               'ggh' : [Unit(
                   datasets['ggH'], [
                       channel_selection(channel, era),
                       ggH125_process_selection(channel, era)],
                       [control_binning[channel][v] for v in set(control_binning[channel].keys()) & set(args.control_plot_set)])],
               'qqh' : [Unit(
                   datasets['qqH'], [
                       channel_selection(channel, era),
                       qqH125_process_selection(channel, era)],
                       [control_binning[channel][v] for v in set(control_binning[channel].keys()) & set(args.control_plot_set)])],
                }
    # Step 1: create units and book actions
    for channel in args.channels:
        nominals[args.era]['datasets'][channel] = get_nominal_datasets(args.era, channel)
        if args.control_plots:
            nominals[args.era]['units'][channel] = get_control_units(channel, args.era, nominals[args.era]['datasets'][channel])
        else:
            nominals[args.era]['units'][channel] = get_analysis_units(channel, args.era, nominals[args.era]['datasets'][channel])

    um = UnitManager()

    jetFakesDS = {
        "et": {"zj", "ttj", "vvj", "w"},
        "mt": {"zj", "ttj", "vvj", "w"},
        "tt": {"zj", "ttj", "vvj", "w"},
        "em": {"w"}
    }
    leptonFakesS = {"zl", "ttl", "vvl"}
    trueTauBkgS = {"ztt", "ttt", "vvt"}
    sm_signalsS = {"ggh", "qqh", "tth", "zh", "wh", "gghww", "qqhww", "zhww", "whww"}
    mssm_signalsS = set("ggh{}".format(mass) for mass in susy_masses[args.era]["ggH"]) | \
                    set("bbh{}_nlo".format(mass) for mass in susy_masses[args.era]["bbH_nlo"])
    signalsS = sm_signalsS | mssm_signalsS
    simulatedProcsDS = {
        chname_: jetFakesDS[chname_] | leptonFakesS | trueTauBkgS | signalsS for chname_ in ["et", "mt", "tt", "em"]
    }
    if args.control_plots:
        signalsS = signalsS & {"ggh", "qqh"}
        simulatedProcsDS = {
            chname_: jetFakesDS[chname_] | leptonFakesS | trueTauBkgS | signalsS for chname_ in ["et", "mt", "tt", "em"]
        }


    for ch_ in args.channels:
        if ch_ == 'mt':
            um.book([unit for d in {'data', 'emb'} | trueTauBkgS | leptonFakesS for unit in nominals[args.era]['units'][ch_][d]], [same_sign, anti_iso_lt])
            um.book([unit for d in jetFakesDS[ch_] for unit in nominals[args.era]['units'][ch_][d]], [same_sign])
            um.book([unit for d in signalsS for unit in nominals[args.era]['units'][ch_][d]])
            if args.skip_systematic_variations:
                pass
            else:
                um.book([unit for d in simulatedProcsDS[ch_] for unit in nominals[args.era]['units'][ch_][d]], [*jet_es, *met_unclustered, *btag_eff, *mistag_eff])
                um.book([unit for d in trueTauBkgS | leptonFakesS | signalsS - {"zl"} for unit in nominals[args.era]['units'][ch_][d]], [*tau_es_3prong, *tau_es_3prong1pizero, *tau_es_1prong, *tau_es_1prong1pizero])
                um.book([unit for d in {'ztt', 'zj', 'zl', 'w'} | signalsS for unit in nominals[args.era]['units'][ch_][d]], [*recoil_resolution, *recoil_response])
                um.book([unit for d in jetFakesDS[ch_] for unit in nominals[args.era]['units'][ch_][d]], [*jet_to_tau_fake])
                um.book([unit for unit in nominals[args.era]['units'][ch_]['zl']], [*mu_fake_es_1prong, *mu_fake_es_1prong1pizero])
                um.book([unit for d in ['ztt', 'zl', 'zj'] for unit in nominals[args.era]['units'][ch_][d]], [*zpt])
                um.book([unit for d in ['ttt', 'ttl', 'ttj'] for unit in nominals[args.era]['units'][ch_][d]], [*top_pt])
                um.book([unit for unit in nominals[args.era]['units'][ch_]['emb']], [*emb_tau_es_3prong, *emb_tau_es_3prong1pizero, *emb_tau_es_1prong, *emb_tau_es_1prong1pizero, *emb_decay_mode_eff, *emb_tau_id_eff_lt, *tau_id_eff_lt])
                um.book([unit for unit in nominals[args.era]['units'][ch_]['ggh']], [*ggh_acceptance])
                um.book([unit for unit in nominals[args.era]['units'][ch_]['qqh']], [*qqh_acceptance])
                um.book([unit for d in {'data', 'emb'} | leptonFakesS | trueTauBkgS for unit in nominals[args.era]['units'][ch_][d]], [*ff_variations_lt])
                # Booking of era dependent uncertainty shapes
                if "2016" in args.era:
                    um.book([unit for d in simulatedProcsDS[ch_] for unit in nominals[args.era]['units'][ch_][d]], [*prefiring, *lep_trigger_eff_mt_2016, *tau_trigger_eff_mt_2016])
                    um.book([unit for unit in nominals[args.era]['units'][ch_]['zl']], [*zll_mt_fake_rate_2016])
                    um.book([unit for unit in nominals[args.era]['units'][ch_]['emb']], [*lep_trigger_eff_mt_emb_2016, *tau_trigger_eff_mt_emb_2016])
                elif "2017" in args.era:
                    um.book([unit for d in simulatedProcsDS[ch_] for unit in nominals[args.era]['units'][ch_][d]], [*prefiring, *lep_trigger_eff_mt_2017, *tau_trigger_eff_mt_2017])
                    um.book([unit for unit in nominals[args.era]['units'][ch_]['zl']], [*zll_mt_fake_rate_2017])
                    um.book([unit for unit in nominals[args.era]['units'][ch_]['emb']], [*lep_trigger_eff_mt_emb_2017, *tau_trigger_eff_mt_emb_2017])
                elif "2018" in args.era:
                    um.book([unit for d in simulatedProcsDS[ch_] for unit in nominals[args.era]['units'][ch_][d]], [*lep_trigger_eff_mt_2018, *tau_trigger_eff_mt_2018])
                    um.book([unit for unit in nominals[args.era]['units'][ch_]['zl']], [*zll_mt_fake_rate_2018])
                    um.book([unit for unit in nominals[args.era]['units'][ch_]['emb']], [*lep_trigger_eff_mt_emb_2018, *tau_trigger_eff_mt_emb_2018])
        elif ch_ == 'et':
            um.book([unit for d in {'data', 'emb'} | trueTauBkgS | leptonFakesS for unit in nominals[args.era]['units'][ch_][d]], [same_sign, anti_iso_lt])
            um.book([unit for d in jetFakesDS[ch_] for unit in nominals[args.era]['units'][ch_][d]], [same_sign])
            um.book([unit for d in signalsS for unit in nominals[args.era]['units'][ch_][d]])
            if args.skip_systematic_variations:
                pass
            else:
                um.book([unit for d in simulatedProcsDS[ch_] for unit in nominals[args.era]['units'][ch_][d]], [*jet_es, *met_unclustered, *btag_eff, *mistag_eff])
                um.book([unit for d in trueTauBkgS | leptonFakesS | signalsS - {"zl"} for unit in nominals[args.era]['units'][ch_][d]], [*tau_es_3prong, *tau_es_3prong1pizero, *tau_es_1prong, *tau_es_1prong1pizero, *tau_id_eff_lt])
                um.book([unit for d in {'ztt', 'zj', 'zl', 'w'} | signalsS for unit in nominals[args.era]['units'][ch_][d]], [*recoil_resolution, *recoil_response])
                um.book([unit for d in jetFakesDS[ch_] for unit in nominals[args.era]['units'][ch_][d]], [*jet_to_tau_fake])
                um.book([unit for unit in nominals[args.era]['units'][ch_]['zl']], [*ele_fake_es_1prong, *ele_fake_es_1prong1pizero])
                um.book([unit for d in ['ztt', 'zl', 'zj'] for unit in nominals[args.era]['units'][ch_][d]], [*zpt])
                um.book([unit for d in ['ttt', 'ttl', 'ttj'] for unit in nominals[args.era]['units'][ch_][d]], [*top_pt])
                um.book([unit for unit in nominals[args.era]['units'][ch_]['emb']], [*emb_tau_es_3prong, *emb_tau_es_3prong1pizero, *emb_tau_es_1prong, *emb_tau_es_1prong1pizero, *emb_decay_mode_eff, *emb_tau_id_eff_lt, *tau_id_eff_lt, *emb_e_es])
                um.book([unit for unit in nominals[args.era]['units'][ch_]['ggh']], [*ggh_acceptance])
                um.book([unit for unit in nominals[args.era]['units'][ch_]['qqh']], [*qqh_acceptance])
                um.book([unit for d in {'data', 'emb'} | trueTauBkgS | leptonFakesS for unit in nominals[args.era]['units'][ch_][d]], [*ff_variations_lt])
                # Booking of era dependent uncertainty shapes
                if "2016" in args.era:
                    um.book([unit for d in simulatedProcsDS[ch_] for unit in nominals[args.era]['units'][ch_][d]], [*prefiring, *lep_trigger_eff_et_2016, *tau_trigger_eff_et_2016])
                    um.book([unit for unit in nominals[args.era]['units'][ch_]['zl']], [*zll_et_fake_rate_2016])
                    um.book([unit for unit in nominals[args.era]['units'][ch_]['emb']], [*lep_trigger_eff_et_emb_2016, *tau_trigger_eff_et_emb_2016])
                elif "2017" in args.era:
                    um.book([unit for d in simulatedProcsDS[ch_] for unit in nominals[args.era]['units'][ch_][d]], [*prefiring, *lep_trigger_eff_et_2017, *tau_trigger_eff_et_2017])
                    um.book([unit for unit in nominals[args.era]['units'][ch_]['zl']], [*zll_et_fake_rate_2017])
                    um.book([unit for unit in nominals[args.era]['units'][ch_]['emb']], [*lep_trigger_eff_et_emb_2017, *tau_trigger_eff_et_emb_2017])
                elif "2018" in args.era:
                    um.book([unit for d in simulatedProcsDS[ch_] for unit in nominals[args.era]['units'][ch_][d]], [*lep_trigger_eff_et_2018, *tau_trigger_eff_et_2018])
                    um.book([unit for unit in nominals[args.era]['units'][ch_]['zl']], [*zll_et_fake_rate_2018])
                    um.book([unit for unit in nominals[args.era]['units'][ch_]['emb']], [*lep_trigger_eff_et_emb_2018, *tau_trigger_eff_et_emb_2018])
        elif ch_ == 'tt':
            um.book([unit for d in {'data', 'emb'} | trueTauBkgS | leptonFakesS for unit in nominals[args.era]['units'][ch_][d]], [same_sign, anti_iso_tt])
            um.book([unit for d in jetFakesDS[ch_] for unit in nominals[args.era]['units'][ch_][d]], [same_sign])
            um.book([unit for d in signalsS for unit in nominals[args.era]['units'][ch_][d]])
            if args.skip_systematic_variations:
                pass
            else:
                um.book([unit for d in simulatedProcsDS[ch_] for unit in nominals[args.era]['units'][ch_][d]], [*jet_es, *met_unclustered, *tau_trigger_eff_tt, *btag_eff, *mistag_eff])
                um.book([unit for d in trueTauBkgS | leptonFakesS | signalsS for unit in nominals[args.era]['units'][ch_][d]], [*tau_es_3prong, *tau_es_3prong1pizero, *tau_es_1prong, *tau_es_1prong1pizero, *tau_id_eff_tt, *tau_trigger_eff_tt])
                um.book([unit for d in {'ztt', 'zj', 'zl', 'w'} | signalsS for unit in nominals[args.era]['units'][ch_][d]], [*recoil_resolution, *recoil_response])
                um.book([unit for d in jetFakesDS[ch_] for unit in nominals[args.era]['units'][ch_][d]], [*jet_to_tau_fake])
                um.book([unit for d in ['ztt', 'zl', 'zj'] for unit in nominals[args.era]['units'][ch_][d]], [*zpt])
                um.book([unit for d in ['ttt', 'ttl', 'ttj'] for unit in nominals[args.era]['units'][ch_][d]], [*top_pt])
                um.book([unit for unit in nominals[args.era]['units'][ch_]['emb']], [*emb_tau_es_3prong, *emb_tau_es_3prong1pizero, *emb_tau_es_1prong, *emb_tau_es_1prong1pizero, *emb_decay_mode_eff, *emb_tau_id_eff_tt, *tau_id_eff_tt])
                um.book([unit for unit in nominals[args.era]['units'][ch_]['ggh']], [*ggh_acceptance])
                um.book([unit for unit in nominals[args.era]['units'][ch_]['qqh']], [*qqh_acceptance])
                um.book([unit for d in {'data', 'emb'} | trueTauBkgS | leptonFakesS for unit in nominals[args.era]['units']['tt'][d]], [*ff_variations_tt])
                if "2016" in args.era:
                    um.book([unit for d in simulatedProcsDS[ch_] for unit in nominals[args.era]['units']['tt'][d]], [*prefiring])
                elif "2017" in args.era:
                    um.book([unit for d in simulatedProcsDS[ch_] for unit in nominals[args.era]['units']['tt'][d]], [*prefiring])
        elif channel == 'em':
            um.book([unit for d in {'data', 'emb'} | simulatedProcsDS[ch_] - signalsS for unit in nominals[args.era]['units'][ch_][d]], [same_sign_em])
            um.book([unit for d in signalsS for unit in nominals[args.era]['units'][ch_][d]])
            if args.skip_systematic_variations:
                pass
            else:
                um.book([unit for d in simulatedProcsDS[ch_] for unit in nominals[args.era]['units'][ch_][d]], [*prefiring, *jet_es, *met_unclustered, *btag_eff, *mistag_eff, *ele_es, *ele_res])
                um.book([unit for d in {'ztt', 'zj', 'zl', 'w'} | signalsS for unit in nominals[args.era]['units'][ch_][d]], [*recoil_resolution, *recoil_response])
                um.book([unit for d in ['ztt', 'zl'] for unit in nominals[args.era]['units'][ch_][d]], [*zpt])
                um.book([unit for d in ['ttt', 'ttl'] for unit in nominals[args.era]['units'][ch_][d]], [*top_pt])
                um.book([unit for unit in nominals[args.era]['units'][ch_]['emb']], [*emb_e_es])
                um.book([unit for unit in nominals[args.era]['units'][ch_]['ggh']], [*ggh_acceptance])
                um.book([unit for unit in nominals[args.era]['units'][ch_]['qqh']], [*qqh_acceptance])
                um.book([unit for d in {'data', 'emb'} | simulatedProcsDS[ch_] - signalsS for unit in nominals[args.era]['units'][ch_][d]], [*qcd_variations_em])
                if "2016" in args.era:
                    um.book([unit for d in simulatedProcsDS[ch_] for unit in nominals[args.era]['units']['em'][d]], [*prefiring])
                elif "2017" in args.era:
                    um.book([unit for d in simulatedProcsDS[ch_] for unit in nominals[args.era]['units']['em'][d]], [*prefiring])



    # Step 2: convert units to graphs and merge them
    g_manager = GraphManager(um.booked_units, True)
    g_manager.optimize(args.optimization_level)
    graphs = g_manager.graphs

    # Step 3: convert to RDataFrame and run the event loop
    r_manager = RunManager(graphs)
    r_manager.run_locally(output_file, args.num_processes, args.num_threads)
    return


if __name__ == "__main__":
    args = parse_arguments()
    setup_logging("produce_shapes.log", logging.INFO)
    main(args)