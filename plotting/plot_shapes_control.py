#!/usr/bin/env python
# -*- coding: utf-8 -*-

import Dumbledraw.dumbledraw as dd
import Dumbledraw.rootfile_parser_inputshapes as rootfile_parser
import Dumbledraw.styles as styles
import ROOT
import argparse
import copy
import yaml
import os

import logging

logger = logging.getLogger("")
from multiprocessing import Pool
from multiprocessing import Process
import multiprocessing

from plots_to_latex import plots_to_latex

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plot categories using Dumbledraw from shapes produced by shape-producer module."
    )
    parser.add_argument(
        "-l", "--linear", action="store_true", help="Enable linear x-axis"
    )
    parser.add_argument("-e", "--era", type=str, required=True, help="Era")
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
        default=1,
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
        "--draw-jet-fake-variation",
        type=str,
        default=None,
        help="Draw variation of jetFakes or QCD in derivation region.",
    )
    parser.add_argument(
        "--blinded", action="store_true", help="if true, no data is plottet"
    )
    parser.add_argument(
        "--tag", type=str, default=None, help="plots are stored in plots/tag/"
    )

    return parser.parse_args()

#loads the plot names from the text file into a list
plot_names = list()
data_file_names = open("data_plot_names.txt", "r")
plot_names = [line.strip() for line in data_file_names.readlines() if line]
data_file_names.close()

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
    variable = info["variable"]
    channel = info["channel"]
    draw_list = list()
    channel_dict = {
        "ee": "#font[42]{#scale[0.85]{ee}}",
        "emet": "#scale[0.85]{e} met",
        "mmet": "#mu met",
        "em": "#scale[0.85]{e}#mu",
        "et": "#font[42]{#scale[0.85]{e}}#tau_{#font[42]{h}}",
        "mm": "#mu#mu",
        "mt": "#mu#tau_{#font[42]{h}}",
        "tt": "#tau_{#font[42]{h}}#tau_{#font[42]{h}}",
    }
    if args.linear == True:
        split_value = 0.1
    else:
        if args.normalize_by_bin_width:
            split_value = 10001
        else:
            split_value = 101

    split_dict = {c: split_value for c in ["et", "mt", "tt", "em", "mm", "mmet", "emet", "ee"]}

    bkg_processes = [
        "VVL",
        "TTL",
        "ZL",
        "jetFakesEMB",
        "EMB",
        "HTT",
    ]  # "ggH125", "qqH125"

    if not args.fake_factor and args.embedding:
        bkg_processes = ["QCDEMB", "VVL", "VVJ", "W", "TTL", "TTJ", "ZJ", "ZL", "EMB"]
    if not args.embedding and args.fake_factor:
        bkg_processes = ["VVT", "VVL", "TTT", "TTL", "ZL", "jetFakes", "ZTT"]
    if not args.embedding and not args.fake_factor:
        bkg_processes = [
            # "VVL",
            "W",
            "TTL",
            "ZL",
        ]
    if args.draw_jet_fake_variation is not None:
        bkg_processes = ["VVL", "TTL", "ZL", "EMB"]
        if not args.fake_factor and args.embedding:
            bkg_processes = ["VVL", "VVJ", "W", "TTL", "TTJ", "ZJ", "ZL", "EMB"]
        if not args.embedding and args.fake_factor:
            bkg_processes = ["VVT", "VVL", "TTT", "TTL", "ZL", "ZTT"]
        if not args.embedding and not args.fake_factor:
            bkg_processes = [
                "VV",
                "W",
                "TTT",
                "TTL",
                "TTJ",
                "ZJ",
                "ZL",
                "ZTT",
            ]
    all_bkg_processes = [b for b in bkg_processes]
    legend_bkg_processes = copy.deepcopy(bkg_processes)
    legend_bkg_processes.reverse()
    if "2016" in args.era:
        era = "Run2016"
    elif "2017" in args.era:
        era = "Run2017"
    elif "2018" in args.era:
        era = "Run2018"
    else:
        logger.critical("Era {} is not implemented.".format(args.era))
        raise Exception

    category = "_".join([channel, variable])
    if args.category_postfix is not None:
        category += "_%s" % args.category_postfix
    rootfile = rootfile_parser.Rootfile_parser(args.input, variable)
    bkg_processes = [b for b in all_bkg_processes]
    if "em" in channel:
        if not args.embedding:
            bkg_processes = ["QCD", "VVT", "VVL", "W", "TTT", "TTL", "ZL", "ZTT"]
        if args.embedding:
            bkg_processes = ["QCDEMB", "VVL", "W", "TTL", "ZL", "EMB"]
        if args.draw_jet_fake_variation is not None:
            if not args.embedding:
                bkg_processes = ["VVT", "VVL", "W", "TTT", "TTL", "ZL", "ZTT"]
            if args.embedding:
                bkg_processes = ["VVL", "W", "TTL", "ZL", "EMB"]

    if "mmet" in channel:
        # bkg_processes = ["ZL", "TTL"]
        bkg_processes = ["ZL", "TTL", "W"]
        # bkg_processes = ["VVL", "W", "TTL", "ZL"]
    elif "emet" in channel:
        bkg_processes = ["ZL", "TTL", "W"]
        # bkg_processes = ["VVL", "W", "TTL", "ZL"]
    elif "mm" in channel:
        bkg_processes = ["TTL", "ZL"]
        # bkg_processes = ["W", "TTL", "ZL"]
        # bkg_processes = ["VVL", "W", "TTL", "ZL"]
    elif "ee" in channel:
        # bkg_processes = ["W", "TTL", "ZL"]
        bkg_processes = ["TTL", "ZL"]
        # bkg_processes = ["VVL", "W", "TTL", "ZL"]

    legend_bkg_processes = copy.deepcopy(bkg_processes)
    legend_bkg_processes.reverse()

    # create plot
    width = 600
    if args.linear == True:
        plot = dd.Plot([0.3, [0.3, 0.28]], "ModTDR", r=0.04, l=0.14, width=width)
    else:
        plot = dd.Plot([0.5, [0.3, 0.28]], "ModTDR", r=0.04, l=0.14, width=width)

    # get background histograms
    total_bkg = None
    if args.draw_jet_fake_variation is None:
        stype = "Nominal"
    else:
        stype = args.draw_jet_fake_variation
    for index, process in enumerate(bkg_processes):
        if index == 0:
            total_bkg = rootfile.get(
                channel, process, args.category_postfix, shape_type=stype
            ).Clone()
            plot.add_hist(
                rootfile.get(channel, process, args.category_postfix, shape_type=stype),
                process,
                "bkg",
            )
        else:
            if process == "HTT":
                HTT = rootfile.get(
                    channel, "ggH125", args.category_postfix, shape_type=stype
                ).Clone()
                HTT.Add(
                    rootfile.get(
                        channel, "qqH125", args.category_postfix, shape_type=stype
                    )
                )
                total_bkg.Add(HTT)
                plot.add_hist(HTT, "HTT", "bkg")
            else:
                total_bkg.Add(
                    rootfile.get(
                        channel, process, args.category_postfix, shape_type=stype
                    )
                )
                plot.add_hist(
                    rootfile.get(
                        channel, process, args.category_postfix, shape_type=stype
                    ),
                    process,
                    "bkg",
                )

        plot.setGraphStyle(process, "hist", fillcolor=styles.color_dict[process])

    plot.add_hist(total_bkg, "total_bkg")

    for index, category in enumerate(plot_names):
        print("PROCESS VALUES: ", category)
        plot.add_hist(
            rootfile.get(channel, "data", category, shape_type=stype),
            category,
        )

    for data_name in plot_names:
        if matchData:
            data_norm = plot.subplot(0).get_hist(data_name).Integral()
            entry_number = plot.subplot(0).get_hist(data_name).GetEntries()
            mc_norm = plot.subplot(0).get_hist("total_bkg").Integral()
            assert mc_norm > 0.
            plot.subplot(0).get_hist("total_bkg").Scale(1/mc_norm)
            plot.subplot(1).get_hist("total_bkg").Scale(1/mc_norm)
            plot.subplot(2).get_hist("total_bkg").Scale(1/mc_norm)
            if data_norm != 0:
                plot.subplot(0).get_hist(data_name).Scale(1/data_norm)
                plot.subplot(1).get_hist(data_name).Scale(1/data_norm)
                plot.subplot(2).get_hist(data_name).Scale(1/data_norm)
            for _proc in bkg_processes:
                plot.subplot(0).get_hist(_proc).Scale(1/mc_norm)
                plot.subplot(1).get_hist(_proc).Scale(1/mc_norm)
                plot.subplot(2).get_hist(_proc).Scale(1/mc_norm)
        else:
            plot.subplot(0).get_hist(data_name).GetXaxis().SetMaxDigits(4) 
            entry_number = plot.subplot(0).get_hist(data_name).GetEntries()
        if entry_number > 50:
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
        plot.subplot(2).setGraphStyle(data_name, "e0", markercolor=styles.color_dict[data_name])

        # set axes limits and labels
        plot.subplot(0).setYlims(
            split_dict[channel],
            max(
                2 * plot.subplot(0).get_hist(data_name).GetMaximum(),
                split_dict[channel] * 2,
            ),
        )

    plot.subplot(2).normalize(["total_bkg"], "total_bkg")
    plot.setGraphStyle(
        "total_bkg", "e2", markersize=0, fillcolor=styles.color_dict["unc"], linecolor=0
    )

    # stack background processes
    plot.create_stack(bkg_processes, "stack")

    # normalize stacks by bin-width
    if args.normalize_by_bin_width:
        plot.subplot(0).normalizeByBinWidth()
        plot.subplot(1).normalizeByBinWidth()

    # plot.subplot(2).setYlims(0.75, 1.55)
    plot.subplot(2).setYlims(0.00, 2)
    if channel == "mmet":
        plot.subplot(0).setLogY()
        #ymax = plot.subplot(0).get_hist("355206").GetMaximum()
        #plot.subplot(0).setYlims(0, ymax*1.5)
        ymax_list = list()
        for i in plot_names:
            ymax_list.append(plot.subplot(0).get_hist(i).GetMaximum())
        ymax = max(ymax_list)
        plot.subplot(0).setYlims(0, ymax*20000000)
        # plot.subplot(0).setXlims(50, 150)
        # plot.subplot(1).setXlims(50, 150)
        # plot.subplot(2).setXlims(50, 150)
    elif channel == "emet":
        plot.subplot(0).setLogY()
        ymax_list = list()
        for i in plot_names:
            ymax_list.append(plot.subplot(0).get_hist(i).GetMaximum())
        ymax = max(ymax_list)
        plot.subplot(0).setYlims(0, ymax*20000000)
        #plot.subplot(0).setYlims(0,1)
        # plot.subplot(0).setXlims(50, 150)
        # plot.subplot(1).setXlims(50, 150)
        # plot.subplot(2).setXlims(50, 150)
    elif channel == "mm":
        #plot.subplot(0).setLogY()
        ymax_list = list()
        for i in plot_names:
            ymax_list.append(plot.subplot(0).get_hist(i).GetMaximum())
        ymax = max(ymax_list)
        plot.subplot(0).setYlims(0, ymax*1.7)
        #plot.subplot(0).setYlims(0, 15000)
        # plot.subplot(0).setXlims(50, 150)
        # plot.subplot(1).setXlims(50, 150)
        # plot.subplot(2).setXlims(50, 150)
    elif channel == "ee":
        #plot.subplot(0).setLogY()
        ymax_list = list()
        for i in plot_names:
            ymax_list.append(plot.subplot(0).get_hist(i).GetMaximum())
        ymax = max(ymax_list)
        plot.subplot(0).setYlims(0, ymax*1.7)
        # plot.subplot(0).setXlims(50, 150)
        # plot.subplot(1).setXlims(50, 150)
        # plot.subplot(2).setXlims(50, 150)

    if args.linear != True:
        plot.subplot(1).setYlims(0.1, split_dict[channel])
        plot.subplot(1).setYlabel("")  # otherwise number labels are not drawn on axis
        #plot.subplot(1).setLogY()
    if variable != None:
        xLabelName = variable.replace("_barrel", "").replace("_endcap", "").replace("_pos", "").replace("_neg", "").replace("_nv015", "").replace("_nv1530", "").replace("_nv3045", "").replace("_isoSR", "").replace("_iso5", "").replace("_iso6", "").replace("_iso7", "").replace("_iso8", "").replace("_iso9", "").replace("_iso10", "").replace("_iso11", "").replace("_iso12", "").replace("_iso13", "")
        if xLabelName in styles.x_label_dict[channel]:
            x_label = styles.x_label_dict[channel][xLabelName]
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
        plot.subplot(0).setYlabel("N_{events}")

    plot.subplot(2).setYlabel("Data/MC")
    plot.subplot(2).setGrid()
    plot.scaleYLabelSize(0.8)
    plot.scaleYTitleOffset(0.8)
    plot.scaleXTitleOffset(1.5)
    plot.scaleXTitleSize(0.5)

    # draw subplots. Argument contains names of objects to be drawn in corresponding order.
    path = "chi_square_data/data_outfile.csv"
    procs_to_draw = ["stack", "total_bkg"]
    for drawing_element in draw_list:
        procs_to_draw.append(drawing_element)
        plot.subplot(0).Draw(procs_to_draw)
        if args.linear != True:
            plot.subplot(1).Draw(procs_to_draw)

        #print chi square info in the terminal
        print("\n Run " + drawing_element + ":")
        plot.subplot(2).get_hist(drawing_element).Fit("pol0","0")

        #store chi square info to a text file that can be manually exported to excel
        chi2 = str(plot.subplot(2).get_hist(drawing_element).GetFunction("pol0").GetChisquare())
        num_dof = str(plot.subplot(2).get_hist(drawing_element).GetFunction("pol0").GetNDF())
        p_0 = str(plot.subplot(2).get_hist(drawing_element).GetFunction("pol0").GetParameter(0)) + " +/- " + str(plot.subplot(2).get_hist(drawing_element).GetFunction("pol0").GetParError(0))
        str_to_write = variable + "," + drawing_element + "," + channel + "," + chi2 + "," + num_dof +"," + p_0 + "\n"
        if os.path.exists(path): 
            with open(path , 'a') as f:
                f.write(str_to_write)
        else:
            with open(path , 'w') as f:
                f.write("variable,run,channel,chi2,dof,p_0\n")
                f.write(str_to_write)
    plot.subplot(2).Draw(procs_to_draw[1:])
            

    # create legends
    suffix = ["", "_top"]
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
        for index in plot_names:
            data_label_name = "Observed" if plot_names[0] == "data" else "Run " + index
            plot.legend(i).add_entry(0, index, data_label_name, "PE2L")
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
    plot.DrawCMS(variable, channel)
    if "2016" in args.era:
        plot.DrawLumi("35.9 fb^{-1} (2016, 13 TeV)")
    elif "2017" in args.era:
        plot.DrawLumi("41.5 fb^{-1} (2017, 13 TeV)")
    elif "2018" in args.era and matchData:
        plot.DrawLumi("(2022, 13.6 TeV)")
    elif "2018" in args.era and not matchData:
        lumiString = "1 pb^{-1}" if len(plot_names) != 1 else args.lumi_label + " fb^{-1}"
        plot.DrawLumi(lumiString + "(2022, 13.6 TeV)")
    else:
        logger.critical("Era {} is not implemented.".format(args.era))
        raise Exception

    posChannelCategoryLabelLeft = None
    plot.DrawChannelCategoryLabel(
        "%s, %s" % (channel_dict[channel], "{cat}".format(cat=args.category_postfix)),
        begin_left=posChannelCategoryLabelLeft,
    )
    # save plot
    if not args.embedding and not args.fake_factor:
        postfix = "fully_classic"
    if args.embedding and not args.fake_factor:
        postfix = "emb_classic"
    if not args.embedding and args.fake_factor:
        postfix = "classic_ff"
    if args.embedding and args.fake_factor:
        postfix = "emb_ff"
    if args.draw_jet_fake_variation is not None:
        postfix = postfix + "_" + args.draw_jet_fake_variation

    if not os.path.exists("plots/%s" % (args.tag)):
        os.makedirs("plots/%s/%s_plots_%s/%s" % (args.tag, args.era, postfix, channel))

    plot.save(
        "plots/%s/%s_plots_%s/%s/%s_%s_%s_%s.%s"
        % (
            args.tag,
            args.era,
            postfix,
            channel,
            args.era,
            channel,
            variable,
            args.category_postfix,
            "pdf",
        )
    )
    plot.save(
        "plots/%s/%s_plots_%s/%s/%s_%s_%s_%s.%s"
        % (
            args.tag,
            args.era,
            postfix,
            channel,
            args.era,
            channel,
            variable,
            args.category_postfix,
            "png",
        )
    )

#function that seperates the variables just like in produce shapes
def seperateVariables(input_list):
    charge_list=["_pos", "_neg"]
    eta_list=["_barrel", "_endcap"]
    iso_list=["_isoSR","_iso5","_iso6","_iso7","_iso8","_iso9","_iso10","_iso11","_iso12","_iso13"]
    npv_list=["_npv015", "_npv1530", "_npv3045"]
    variable_list = list()
    for var in input_list:
        for charge in charge_list:
            for eta in eta_list:
                for iso in iso_list:
                    variable_list.append(var+charge+eta+iso)
    return variable_list


if __name__ == "__main__":
    args = parse_arguments()
    setup_logging("{}_plot_shapes.log".format(args.era), logging.DEBUG)
    variables = args.variables.split(",")
    seperateVariable = bool(int(args.seperate_variables))
    if seperateVariable: variables = seperateVariables(variables)
    channels = args.channels.split(",")
    writeLatex = bool(int(args.write_to_latex))
    matchData = bool(int(args.match_data))
    infolist = []

    if not args.embedding and not args.fake_factor:
        postfix = "fully_classic"
    if args.embedding and not args.fake_factor:
        postfix = "emb_classic"
    if not args.embedding and args.fake_factor:
        postfix = "classic_ff"
    if args.embedding and args.fake_factor:
        postfix = "emb_ff"
    for ch in channels:
        for v in variables:
            infolist.append({"args": args, "channel": ch, "variable": v})
    pool = Pool(1)
    pool.map(main, infolist)
    if writeLatex:
         plots_to_latex(args.tag, args.era, postfix, args.channels)
    # for info in infolist:
    #     main(info)

    