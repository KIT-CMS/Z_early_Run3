#!/usr/bin/env python
# -*- coding: utf-8 -*-

import Dumbledraw.dumbledraw as dd
import Dumbledraw.rootfile_parser_inputshapes as rootfile_parser
import Dumbledraw.styles as styles
import ROOT
import argparse
import copy
import yaml
import os, sys
from time import sleep
import gc

from config.common.variables import seperate_var, get_base_name, get_all_variables
from config.shapes.group_runs import get_subera_boundaries

import logging
logger = logging.getLogger("")
from multiprocessing import Pool
from multiprocessing import Process
import multiprocessing

# from plots_to_latex import plots_to_latex

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plot categories using Dumbledraw from shapes produced by shape-producer module."
    )
    parser.add_argument(
        "-l", "--linear", action="store_true", help="Enable linear x-axis"
    )
    parser.add_argument("-e", "--era", type=str, required=True, help="Era")
    parser.add_argument("--subera", type=str, default="Run2022", help="Experiment sub-era.")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="ROOT file with shapes of processes",
    )
    parser.add_argument(
        "--variables",
        type=str,
        default=None,
        help="Enable control plotting for given variable",
    )
    parser.add_argument(
        "--category-postfix",
        type=str,
        default=None,
        help="Enable control plotting for given category_postfix. Structure of a category: <variable>_<postfix>",
    )
    parser.add_argument(
        "--match-data", 
        default=0,
        help="When this it true it will normalize the data to the MC"
    )
    parser.add_argument(
        "--seperate-variables", 
        default=0,
        help="When this is true it will seperate variables based off suffixes you assign"
    )
    parser.add_argument(
        "--lumi-label", 
        type=str,
        default=None,
        help="Determines what the luminosity label reads when the runs are summed"
    )
    parser.add_argument(
        "--write-to-latex",
        default=0,
        help="If true this will automatically write all the plots to latex slides for presenting",
    )
    parser.add_argument(
        "--channels",
        type=str,
        default=None,
        help="Enable control plotting for given variable",
    )
    parser.add_argument(
        "--normalize-by-bin-width",
        action="store_true",
        help="Normalize plots by bin width",
    )
    parser.add_argument(
        "--fake-factor", action="store_true", help="Fake factor estimation method used"
    )
    parser.add_argument(
        "--embedding", action="store_true", help="Fake factor estimation method used"
    )
    parser.add_argument(
        "--syst",
        type=str,
        default="Nominal",
        help="Draw sys variation",
    )
    parser.add_argument(
        "--blinded", action="store_true", help="if true, no data is plottet"
    )
    parser.add_argument(
        "--tag", type=str, default=None, help="plots are stored in plots/tag/"
    )
    parser.add_argument(
        "--doQCD",
        action="store_true",
        help="run QCD templates",
    )
    parser.add_argument(
        "--doZpt",
        action="store_true",
        help="run over Zpt bins",
    )
    parser.add_argument(
        "--doLepCorrBins",
        action="store_true",
        help="run over lepton momentum correction bins",
    )
    parser.add_argument("--plot-postfix", type=str, default="", help="plot name postfix")

    return parser.parse_args()

# loads the plot names from the text file into a list
# plot_names = list()
# data_file_names = open("data_plot_names.txt", "r")
# plot_names = [line.strip() for line in data_file_names.readlines() if line]
# data_file_names.close()

def setup_logging(output_file, level=logging.DEBUG):
    logger.setLevel(level)
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    file_handler = logging.FileHandler(output_file, "w")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

def main(info):
    args = info["args"]
    subera = args.subera
    variable = info["variable"]
    channel = info["channel"]
    channel_base = channel.replace("_corr", "")
    draw_list = list()
    channel_dict = {
        "ee": "#font[42]{#scale[1.0]{ee}}",
        "emet": "#font[42]{#scale[1.0]{single-e}}",
        "mmet": "#font[42]{#scale[1.0]{single-#mu}}",
        "mm": "#font[42]{#scale[1.0]{#mu#mu}}",
    }
    doLogy = info["doLogy"]

    # HERE
    subera_boundaries = get_subera_boundaries()
    plot_names = [
        "%d-%d" % (subera_boundaries[args.subera][0][0], subera_boundaries[args.subera][0][1])

        # "355862-357900"

        # "data"

        # "355862-356446",
        # "356523-357482",
        # "357538-357732",
        # "357734-357900",

        # "355862-356578",
        # "356580-356969",
        # "356970-357271",
        # "357328-357442",
        # "357447-357696",
        # "357697-357778",
        # "357779-357900",

        # "355872-356075",
        # "356076-356323",
        # "356371-356378",
        # "356381",
        # "356383-356433",
        # "356434-356446"
    ]

    if args.linear == True:
        split_value = 0.1
    else:
        if args.normalize_by_bin_width:
            split_value = 10001
        else:
            split_value = 101

    split_dict = {c: split_value for c in ["mm", "mmet", "emet", "ee"]}

    category = "_".join([channel, variable])
    if args.category_postfix is not None:
        category += "_%s" % args.category_postfix
    rootfile = rootfile_parser.Rootfile_parser(args.input, variable)

    # create plot
    width = 600
    if args.linear == True:
        plot = dd.Plot([0.3, [0.3, 0.28]], "ModTDR", r=0.04, l=0.14, width=width)
    else:
        plot = dd.Plot([0.5, [0.3, 0.28]], "ModTDR", r=0.04, l=0.14, width=width)

    # get background histograms
    if args.syst is None:
        stype = "Nominal"
    else:
        stype = args.syst

    # HERE tmp
    # bkg_processes = ["DYtau", "VV", "ST", "TT", "DY"]
    bkg_processes = ["Wtau", "DYtau", "VV", "ST", "TT", "W", "DY"]

    is_first_total = True
    total_bkg = None
    
    is_first_ewk = True
    ewk_bkg = None

    ewk_processes = None
    if channel_base in ["mm", "ee"]:
        ewk_processes = ["Wtau", "DYtau", "VV", "ST", "W"]
    elif channel_base in ["mmet", "emet"]:
        ewk_processes = ["Wtau", "DYtau", "VV", "ST", "DY"]

    # HERE tmp
    # stack_processes = bkg_processes
    stack_processes = ["EWK"] + [p for p in bkg_processes if p not in ewk_processes]

    # print("stack_processes:", stack_processes)
    legend_bkg_processes = copy.deepcopy(stack_processes)
    legend_bkg_processes.reverse()

    for process in bkg_processes:
        if is_first_total:
            total_bkg = rootfile.get(
                channel, process, plot_names[0], shape_type=stype
            ).Clone()
            plot.add_hist(
                rootfile.get(channel, process, plot_names[0], shape_type=stype),
                process,
                "bkg",
            )
            is_first_total = False
        else:
            total_bkg.Add(
                rootfile.get(
                    channel, process, plot_names[0], shape_type=stype
                )
            )
            plot.add_hist(
                rootfile.get(
                    channel, process, plot_names[0], shape_type=stype
                ),
                process,
                "bkg",
            )

        if process in ewk_processes:
            if is_first_ewk:
                ewk_bkg = rootfile.get(
                    channel, process, plot_names[0], shape_type=stype
                ).Clone()
                is_first_ewk = False
            else:
                ewk_bkg.Add(
                    rootfile.get(
                        channel, process, plot_names[0], shape_type=stype
                    )
                )

        # if is_first_total:
        #     total_bkg = rootfile.get(
        #         channel, process, args.category_postfix, shape_type=stype
        #     ).Clone()
        #     plot.add_hist(
        #         rootfile.get(channel, process, args.category_postfix, shape_type=stype),
        #         process,
        #         "bkg",
        #     )
        #     is_first_total = False
        # else:
        #     total_bkg.Add(
        #         rootfile.get(
        #             channel, process, args.category_postfix, shape_type=stype
        #         )
        #     )
        #     plot.add_hist(
        #         rootfile.get(
        #             channel, process, args.category_postfix, shape_type=stype
        #         ),
        #         process,
        #         "bkg",
        #     )

        # if process in ewk_processes:
        #     if is_first_ewk:
        #         ewk_bkg = rootfile.get(
        #             channel, process, args.category_postfix, shape_type=stype
        #         ).Clone()
        #         is_first_ewk = False
        #     else:
        #         ewk_bkg.Add(
        #             rootfile.get(
        #                 channel, process, args.category_postfix, shape_type=stype
        #             )
        #         )

        plot.setGraphStyle(process, "hist", fillcolor=styles.color_dict[process])

    plot.add_hist(total_bkg, "total_bkg")
    plot.add_hist(ewk_bkg, "EWK", "bkg")
    plot.setGraphStyle("EWK", "hist", fillcolor=styles.color_dict["EWK"])

    for data_name in plot_names:
        if args.category_postfix != None and args.category_postfix != "None":
            plot.add_hist(
                rootfile.get(channel+"-"+args.category_postfix, "data", data_name, shape_type="Nominal"),
                data_name,
            )
        else:
            plot.add_hist(
                rootfile.get(channel, "data", data_name, shape_type="Nominal"),
                data_name,
            )

    if matchData:
        mc_norm = plot.subplot(0).get_hist("total_bkg").Integral()
        assert mc_norm > 0.
        plot.subplot(0).get_hist("total_bkg").Scale(1/mc_norm)
        plot.subplot(1).get_hist("total_bkg").Scale(1/mc_norm)
        plot.subplot(2).get_hist("total_bkg").Scale(1/mc_norm)
        for _proc in stack_processes:
            plot.subplot(0).get_hist(_proc).Scale(1/mc_norm)
            plot.subplot(1).get_hist(_proc).Scale(1/mc_norm)
            plot.subplot(2).get_hist(_proc).Scale(1/mc_norm)

    for data_name in plot_names:
        if matchData:
            data_norm = plot.subplot(0).get_hist(data_name).Integral()
            entry_number = plot.subplot(0).get_hist(data_name).GetEntries()
            if data_norm != 0:
                plot.subplot(0).get_hist(data_name).Scale(1/data_norm)
                plot.subplot(1).get_hist(data_name).Scale(1/data_norm)
                plot.subplot(2).get_hist(data_name).Scale(1/data_norm)
        else:
            plot.subplot(0).get_hist(data_name).GetXaxis().SetMaxDigits(4) 
            entry_number = plot.subplot(0).get_hist(data_name).GetEntries()
        if entry_number > 10000:
            draw_list.append(data_name)
        if args.blinded:
            plot.subplot(0).setGraphStyle(data_name, "e0", markercolor=styles.color_dict[data_name], markersize=0, linewidth=0)
            plot.subplot(0).setGraphStyle(data_name, "e0", markercolor=styles.color_dict[data_name], markersize=0, linewidth=0)
        else:
            plot.subplot(0).setGraphStyle(data_name, "e0", markercolor=styles.color_dict[data_name])
            plot.subplot(0).setGraphStyle(data_name, "e0", markercolor=styles.color_dict[data_name])
        if args.linear:
            pass
        else:
            if args.blinded:
                plot.subplot(1).setGraphStyle(data_name, "e0", markersize=0, linewidth=0)
            else:
                plot.subplot(1).setGraphStyle(data_name, "e0")

        plot.subplot(2).normalize([data_name], "total_bkg")
        plot.subplot(2).setGraphStyle(data_name, "e0", markercolor=styles.color_dict[data_name])

        # set axes limits and labels
        plot.subplot(0).setYlims(
            split_dict[channel_base],
            max(
                2 * plot.subplot(0).get_hist(data_name).GetMaximum(),
                split_dict[channel_base] * 2,
            ),
        )

    plot.subplot(2).normalize(["total_bkg"], "total_bkg")
    plot.setGraphStyle(
        "total_bkg", "e2", markersize=0, fillcolor=styles.color_dict["unc"], linecolor=0
    )

    # stack background processes
    plot.create_stack(stack_processes, "stack")

    # normalize stacks by bin-width
    if args.normalize_by_bin_width:
        plot.subplot(0).normalizeByBinWidth()
        plot.subplot(1).normalizeByBinWidth()

    plot.subplot(2).setYlims(0.7, 1.30)
    # plot.subplot(2).setYlims(0.01, 1.99)
    ymin_list = list()
    ymax_list = list()
    for i in plot_names:
        ymax_list.append(plot.subplot(0).get_hist(i).GetMaximum())
        ymin_list.append(plot.subplot(0).get_hist(i).GetMinimum(1e-5))
    ymax = max(ymax_list)
    ymin = min(ymin_list)
    ymin = min(ymin, 1e-2)
    plot.subplot(0).setYlims(0.5, ymax*1.7)

    if doLogy:
        plot.subplot(0).setLogY()
        if len(plot_names) == 1:
            plot.subplot(0).setYlims(0.5, ymax*1.e4)
        else:
            plot.subplot(0).setYlims(0.4*ymin, ymax*1.e4)
    if matchData:
        plot.subplot(0).setYlims(1e-5, 1e2)

    if args.linear != True:
        plot.subplot(1).setYlims(0.1, split_dict[channel_base])
        plot.subplot(1).setYlabel("")  # otherwise number labels are not drawn on axis
        #plot.subplot(1).setLogY()
    if variable != None:
        xLabelName = get_base_name(variable)
        if xLabelName in styles.x_label_dict[channel_base]:
            x_label = styles.x_label_dict[channel_base][xLabelName]
        else:
            x_label = variable
        plot.subplot(2).setXlabel(x_label)
    else:
        plot.subplot(2).setXlabel("NN output")
    if args.normalize_by_bin_width:
        plot.subplot(0).setYlabel("dN/d(NN output)")
    elif matchData:
        plot.subplot(0).setYlabel("a.u.")
    else:
        plot.subplot(0).setYlabel("# events")

    plot.subplot(2).setYlabel("Data/MC")
    plot.subplot(2).setGrid()
    plot.scaleYLabelSize(0.8)
    plot.scaleYTitleOffset(0.8)
    plot.scaleXTitleOffset(1.2)
    plot.scaleXTitleSize(0.8)

    # draw subplots. Argument contains names of objects to be drawn in corresponding order.
    procs_to_draw = ["stack", "total_bkg"]
    for drawing_element in draw_list:
        procs_to_draw.append(drawing_element)

    plot.subplot(0).Draw(procs_to_draw)
    if args.linear != True:
        plot.subplot(1).Draw(procs_to_draw)

    # path = "chi_square_data/data_outfile.csv"
    # for drawing_element in draw_list:
    #     #print chi square info in the terminal
    #     print("Run " + drawing_element + ":")
    #     plot.subplot(2).get_hist(drawing_element).Fit("pol0","0")

    #     #store chi square info to a text file that can be manually exported to excel
    #     chi2 = str(plot.subplot(2).get_hist(drawing_element).GetFunction("pol0").GetChisquare())
    #     num_dof = str(plot.subplot(2).get_hist(drawing_element).GetFunction("pol0").GetNDF())
    #     p_0 = str(plot.subplot(2).get_hist(drawing_element).GetFunction("pol0").GetParameter(0)) + " +/- " + str(plot.subplot(2).get_hist(drawing_element).GetFunction("pol0").GetParError(0))
    #     str_to_write = variable + "," + drawing_element + "," + channel + "," + chi2 + "," + num_dof +"," + p_0 + "\n"
    #     if os.path.exists(path): 
    #         with open(path , 'a') as f:
    #             f.write(str_to_write)
    #     else:
    #         with open(path , 'w') as f:
    #             f.write("variable,run,channel,chi2,dof,p_0\n")
    #             f.write(str_to_write)

    plot.subplot(2).Draw(procs_to_draw[1:])
            

    # create legends
    for i in range(2):
        plot.add_legend(width=0.5, height=0.15)
        for process in legend_bkg_processes:
            plot.legend(i).add_entry(
                0,
                process,
                styles.legend_label_dict[
                    process.replace("TTL", "TT").replace("VVL", "VV").replace("NLO", "")
                ],
                "f",
            )
        plot.legend(i).add_entry(0, "total_bkg", "Bkg. stat. unc.", "f")
        for run_num in draw_list:
            data_label_name = "Observed" if plot_names[0] == "data" else "Run " + run_num
            plot.legend(i).add_entry(0, run_num, data_label_name, "PE2L")
        plot.legend(i).setNColumns(2)
    plot.legend(0).Draw()
    plot.legend(1).setAlpha(0.0)
    plot.legend(1).Draw()

    if plot_names[0] == "data":
        for i in range(2):
            plot.add_legend(reference_subplot=2, pos=1, width=0.8, height=0.06)
            plot.legend(i + 2).add_entry(0, "data", "Observed", "PE2L")
            plot.legend(i + 2).add_entry(0, "total_bkg", "Bkg. stat. unc.", "f")
            plot.legend(i + 2).setNColumns(4)
        plot.legend(2).Draw()
        plot.legend(3).setAlpha(0.0)
        plot.legend(3).Draw()

    # draw additional labels
    plot.DrawCMS(variable, channel_base)
    if "2016" in args.era:
        plot.DrawLumi("35.9 fb^{-1} (2016, 13 TeV)")
    elif "2017" in args.era:
        plot.DrawLumi("41.5 fb^{-1} (2017, 13 TeV)")
    elif "2018" in args.era:
        plot.DrawLumi("59.8 fb^{-1} (2018, 13 TeV)")
    elif "2022" in args.era:
        if matchData or len(plot_names) != 1:
            plot.DrawLumi("(2022, 13.6 TeV)")
        else:
            lumiString = "{:.2f}".format(float(args.lumi_label))+" fb^{-1}"
            plot.DrawLumi(lumiString + " (2022, 13.6 TeV)")
    else:
        logger.critical("Era {} is not implemented.".format(args.era))
        raise Exception

    posChannelCategoryLabelLeft = None
    plot.DrawChannelCategoryLabel(
        "%s" % (channel_dict[channel_base]),  # "{cat}".format(cat=args.category_postfix)
        begin_left=posChannelCategoryLabelLeft,
    )

    if matchData:
        plot.DrawText(0.17, 0.73, "#bf{#font[42]{Normalized to 1}}", 0.03)
    else:
        if len(plot_names) != 1:
            plot.DrawText(0.17, 0.73, "#bf{#font[42]{Normalized to 1 pb^{-1}}}", 0.03)
        else:
            # plot.DrawText(0.17, 0.73, "#bf{#font[42]{Normalized to Z-counting lumi}}", 0.03)
            plot.DrawText(0.17, 0.73, "#bf{#font[42]{Normalized to online lumi}}", 0.03)

    # save plot
    # if not os.path.exists("plots/%s/%s/%s" % (args.tag, args.era, channel)):
    #     os.makedirs("plots/%s/%s/%s" % (args.tag, args.era, channel))

    out_name_base = "plots/%s/%s/%s/%s/%s_%s_%s_%s" % (
        args.tag,
        args.era,
        channel,
        args.syst,
        args.era+"_"+plot_names[0],
        channel,
        args.syst,
        variable,
    )
    if args.category_postfix != None and args.category_postfix != "None":
        out_name_base = out_name_base + "_" + args.category_postfix
    out_name_base = out_name_base+args.plot_postfix
    if doLogy:
        out_name_base = out_name_base + "_log"

    out_name_base = out_name_base.replace("-", "_")

    plot.save(
        "%s.%s"
        % (
            out_name_base,
            "pdf",
        )
    )
    # plot.save(
    #     "%s.%s"
    #     % (
    #         out_name_base,
    #         "png",
    #     )
    # )

if __name__ == "__main__":
    args = parse_arguments()
    
    # setup_logging("{}_plot_shapes.log".format(args.era), logging.DEBUG)
    setup_logging("{}_plot_shapes.log".format(args.era), logging.INFO)

    channels = args.channels.split(",")

    variable_dict = {}
    variables_base = args.variables.split(",") if args.variables is not None else []
    if len(variables_base) > 0:
        if bool(int(args.seperate_variables)):
            variable_dict = {
                "mm": seperate_var(variables_base, doZpt = args.doZpt, doQCD = args.doQCD, doLepCorrBins = args.doLepCorrBins),
                "mmet": seperate_var(variables_base, doZpt = args.doZpt, doQCD = args.doQCD, doLepCorrBins = args.doLepCorrBins),
                "ee": seperate_var(variables_base, doZpt = args.doZpt, doQCD = args.doQCD, doLepCorrBins = args.doLepCorrBins),
                "emet": seperate_var(variables_base, doZpt = args.doZpt, doQCD = args.doQCD, doLepCorrBins = args.doLepCorrBins),
            }
        else:
            variable_dict = {
                "mm": variables_base,
                "mmet": variables_base,
                "ee": variables_base,
                "emet": variables_base,
            }
    else:
        variable_dict_base = get_all_variables(args.doZpt, args.doQCD, args.doLepCorrBins)
        if bool(int(args.seperate_variables)):
            for _ch in channels:
                _ch_base = _ch.replace("_corr", "")
                variable_dict[_ch_base] = seperate_var(variable_dict_base[_ch_base], doZpt = args.doZpt, doQCD = args.doQCD, doLepCorrBins = args.doLepCorrBins)
        else:
            variable_dict = variable_dict_base

    writeLatex = bool(int(args.write_to_latex))
    matchData = bool(int(args.match_data))
    infolist = []

    for _ch in channels:
        _ch_base = _ch.replace("_corr", "")
        outdir = "plots/%s/%s/%s/%s" % (args.tag, args.era, _ch, args.syst)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        for v in variable_dict[_ch_base]:
            infolist.append({"args": args, "channel": _ch, "variable": v, "doLogy": True})
            infolist.append({"args": args, "channel": _ch, "variable": v, "doLogy": False})

    # nthread = min(len(infolist), 64)
    # pool = Pool(nthread)
    # pool.map(main, infolist)

    # if writeLatex:
    #     plots_to_latex(args.tag, args.era, args.channels)
    for info in infolist:
        # print(info["variable"])
        # if "_zpt" in info["variable"] and int(info["variable"].split("_zpt")[1]) < 19:
        #     continue
        main(info)
        sys.stdout.flush()
        # gc.collect()
        sleep(0.1)

    
