#!/usr/bin/env python
import argparse
import logging
import pickle
import re
import yaml
import os

from ntuple_processor import Histogram
from ntuple_processor import (
    dataset_from_crownoutput,
    dataset_from_artusoutput,
    Unit,
    UnitManager,
    GraphManager,
    RunManager,
)

from config.shapes.channel_selection import channel_selection
from config.shapes.run_selection import run_selection
from config.shapes.run_selection import run_group_selection
from config.shapes.group_runs import run_lumi_selection
from config.shapes.file_names import files
from config.shapes.process_selection import (
    DY_process_selection,
    W_process_selection,
    TT_process_selection,
    ST_process_selection,
    VV_process_selection,
    DYtau_process_selection,
    Wtau_process_selection,
)
from config.shapes.data_selection import (
    data_process_selection,
    data_group_process_selection,
)


# Muon ID weight uncertainties
from config.shapes.variations import (
    mu_id_weight,
)

from config.shapes.control_binning import get_control_binning, seperateVariables

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
        description="Produce shapes for the legacy NMSSM analysis."
    )
    parser.add_argument("--era", required=True, type=str, help="Experiment era.")
    parser.add_argument(
        "--channels",
        default=[],
        type=lambda channellist: [channel for channel in channellist.split(",")],
        help="Channels to be considered, seperated by a comma without space",
    )
    parser.add_argument(
        "--directory", required=True, type=str, help="Directory with Artus outputs."
    )
    parser.add_argument(
        "--et-friend-directory",
        type=str,
        default=[],
        nargs="+",
        help="Directories arranged as Artus output and containing a friend tree for et.",
    )
    parser.add_argument(
        "--mt-friend-directory",
        type=str,
        default=[],
        nargs="+",
        help="Directories arranged as Artus output and containing a friend tree for mt.",
    )
    parser.add_argument(
        "--tt-friend-directory",
        type=str,
        default=[],
        nargs="+",
        help="Directories arranged as Artus output and containing a friend tree for tt.",
    )
    parser.add_argument(
        "--mmet-friend-directory",
        type=str,
        default=[],
        nargs="+",
        help="Directories arranged as Artus output and containing a friend tree for mmet.",
    )
    parser.add_argument(
        "--emet-friend-directory",
        type=str,
        default=[],
        nargs="+",
        help="Directories arranged as Artus output and containing a friend tree for emet.",
    )
    parser.add_argument(
        "--mm-friend-directory",
        type=str,
        default=[],
        nargs="+",
        help="Directories arranged as Artus output and containing a friend tree for mm.",
    )
    parser.add_argument(
        "--ee-friend-directory",
        type=str,
        default=[],
        nargs="+",
        help="Directories arranged as Artus output and containing a friend tree for ee.",
    )
    parser.add_argument(
        "--em-friend-directory",
        type=str,
        default=[],
        nargs="+",
        help="Directories arranged as Artus output and containing a friend tree for em.",
    )
    parser.add_argument(
        "--optimization-level",
        default=2,
        type=int,
        help="Level of optimization for graph merging.",
    )
    parser.add_argument(
        "--num-processes", default=1, type=int, help="Number of processes to be used."
    )
    parser.add_argument(
        "--num-threads", default=1, type=int, help="Number of threads to be used."
    )
    parser.add_argument(
        "--skip-systematic-variations",
        action="store_true",
        help="Do not produce the systematic variations.",
    )
    parser.add_argument(
        "--output-file",
        required=True,
        type=str,
        help="ROOT file where shapes will be stored.",
    )
    parser.add_argument(
        "--control-plots",
        action="store_true",
        help="Produce shapes for control plots. Default is production of analysis shapes.",
    )
    parser.add_argument(
        "--run-plot",
        default=1,
        help="If true this will stack all runs on that same plot rather than separate run by run plots",
    )
    parser.add_argument(
        "--group-runs",
        default=0,
        help="If true this will group runs by luminosity which you can set below",
    )
    parser.add_argument(
        "--seperate-variables",
        default=0,
        help="If true this will seperate variables",
    )
    parser.add_argument(
        "--control-plot-set",
        default="pt_1",
        type=lambda varlist: [variable for variable in varlist.split(",")],
        help="Variables the shapes should be produced for.",
    )
    parser.add_argument(
        "--total-lumi",
        default=0.001,
        help="Luminosity value of the total runs summed",
    )
    parser.add_argument(
        "--run-list",
        default='',
        type=lambda runlist: [run_number for run_number in runlist.split(",")],
        help="List of run numbers being used",
    )
    parser.add_argument(
        "--lumi-list",
        default='',
        type=lambda lumilist: [lumi for lumi in lumilist.split(",")],
        help="List of luminosities which correspond to the run numbers",
    )
    parser.add_argument(
        "--only-create-graphs",
        action="store_true",
        help="Create and optimise graphs and create a pkl file containing the graphs to be processed.",
    )
    parser.add_argument(
        "--process-selection",
        default=None,
        type=lambda proclist: set([process for process in proclist.split(",")]),
        help="Subset of processes to be processed.",
    )
    parser.add_argument(
        "--graph-dir",
        default=None,
        type=str,
        help="Directory the graph file is written to.",
    )
    parser.add_argument(
        "--classdict",
        default=None,
        type=str,
        help="path to config file from NN training to extract the classes",
    )
    parser.add_argument(
        "--ntuple_type", default="artus", type=str, help="Options: crown or artus"
    )
    parser.add_argument(
        "--enable-booking-check",
        action="store_true",
        help="Enables check for double actions during booking. Takes long for all variations.",
    )
    return parser.parse_args()

def main(args):
    # Parse given arguments.
    friend_directories = {
        "mm": args.mm_friend_directory,
        "ee": args.ee_friend_directory,
        "mmet": args.mmet_friend_directory,
        "emet": args.emet_friend_directory,
    }

    if ".root" in args.output_file:
        output_file = args.output_file
        log_file = args.output_file.replace(".root", ".log")
    else:
        output_file = "{}.root".format(args.output_file)
        log_file = "{}.log".format(args.output_file)

    nominals = {}
    nominals[args.era] = {}
    nominals[args.era]["datasets"] = {}
    nominals[args.era]["units"] = {}
    runPlot = bool(int(args.run_plot))
    groupRuns = float(args.group_runs)
    totalLumi = float(args.total_lumi)

    def get_nominal_datasets(era, channel):
        datasets = dict()

        def filter_friends(dataset, friend):
            if re.match("(gg|qq|tt|w|z|v)h", dataset.lower()):
                if "FakeFactors" in friend or "EMQCDWeights" in friend:
                    return False
            elif re.match("data", dataset.lower()):
                if "crosssection" in friend:
                    return False
            return True

        for key, names in files[era][channel].items():
            datasets[key] = dataset_from_crownoutput(
                key,
                names,
                args.era,
                channel,
                channel + "_nominal",
                args.directory,
                [
                    fdir
                    for fdir in friend_directories[channel]
                    if filter_friends(key, fdir)
                ],
            )
        return datasets

    def get_control_units(channel, era, datasets, control_binning):
        #adds the MC samples first to the dictionary of histograms to be made
        data_dict = {
            "dy": [
                Unit(
                    datasets["DY"],
                    [
                        channel_selection(channel, era),
                        DY_process_selection(channel, era, runPlot, totalLumi)
                    ],
                    [
                        control_binning[channel][v]
                        for v in set(control_binning[channel].keys())
                        & set(variable_list)
                    ],
                )
            ],
            "w": [
                Unit(
                    datasets["W"],
                    [
                        channel_selection(channel, era),
                        W_process_selection(channel, era, runPlot, totalLumi),
                    ],
                    [
                        control_binning[channel][v]
                        for v in set(control_binning[channel].keys())
                        & set(variable_list)
                    ],
                )
            ],
            "tt": [
                Unit(
                    datasets["TT"],
                    [
                        channel_selection(channel, era),
                        TT_process_selection(channel, era, runPlot, totalLumi),
                    ],
                    [
                        control_binning[channel][v]
                        for v in set(control_binning[channel].keys())
                        & set(variable_list)
                    ],
                )
            ],
            "st": [
                Unit(
                    datasets["ST"],
                    [
                        channel_selection(channel, era),
                        ST_process_selection(channel, era, runPlot, totalLumi),
                    ],
                    [
                        control_binning[channel][v]
                        for v in set(control_binning[channel].keys())
                        & set(variable_list)
                    ],
                )
            ],
            "vv": [
                Unit(
                    datasets["VV"],
                    [
                        channel_selection(channel, era),
                        VV_process_selection(channel, era, runPlot, totalLumi),
                    ],
                    [
                        control_binning[channel][v]
                        for v in set(control_binning[channel].keys())
                        & set(variable_list)
                    ],
                )
            ],
            "ewktau": [
                Unit(
                    datasets["DYtau"],
                    [
                        channel_selection(channel, era),
                        DYtau_process_selection(channel, era, runPlot, totalLumi)
                    ],
                    [
                        control_binning[channel][v]
                        for v in set(control_binning[channel].keys())
                        & set(variable_list)
                    ],
                ),
                Unit(
                    datasets["Wtau"],
                    [
                        channel_selection(channel, era),
                        Wtau_process_selection(channel, era, runPlot, totalLumi),
                    ],
                    [
                        control_binning[channel][v]
                        for v in set(control_binning[channel].keys())
                        & set(variable_list)
                    ],
                )
            ],
        }

        # sums the runs and returns data as its own histogram
        if not runPlot:
            data_dict["data"] = [Unit(
                datasets["data"],
                [
                    channel_selection(channel, era),
                ],
                [
                    control_binning[channel][v]
                    for v in set(control_binning[channel].keys())
                    & set(variable_list)
                ],
            )]
        # seperates run by run and plots them seperatly
        else:
            run_list = args.run_list
            lumi_list = args.lumi_list
            run_lumi = {run_list[i]: float(lumi_list[i]) for i in range(len(run_list))}
            temp_prev_run_number = float(run_list[0])
            list_of_run_lists = list()
            if groupRuns:
                while temp_prev_run_number < float(run_list[-1]):
                    list_of_run_lists.append(run_lumi_selection(temp_prev_run_number, run_list, run_lumi, groupRuns)) 
                    temp_prev_run_number = float(list_of_run_lists[-1][-1])
                for run_list in list_of_run_lists:
                    if run_list[0] == run_list[-1]:
                        data_dict[run_list[0]] = [Unit(
                            datasets["data"],
                            [
                                channel_selection(channel, era),
                                run_selection(run_list[0]),
                                data_process_selection(channel, era, run_list[0], run_lumi),
                            ],
                            [
                                control_binning[channel][v]
                                for v in set(control_binning[channel].keys())
                                & set(variable_list)
                            ],
                        )]
                    else:
                        data_dict[run_list[0] + "-" + run_list[-1]] = [Unit(
                                datasets["data"],
                                [
                                    channel_selection(channel, era),
                                    run_group_selection(run_list[0], run_list[-1]),
                                    data_group_process_selection(channel, era, run_list, run_lumi)
                                ],
                                [
                                    control_binning[channel][v]
                                    for v in set(control_binning[channel].keys())
                                    & set(variable_list)
                                ],
                            )]
            else:
                for run_number in run_list:
                    data_dict[run_number] = [Unit(
                        datasets["data"],
                        [
                            channel_selection(channel, era),
                            run_selection(run_number),
                            data_process_selection(channel, era, run_number, run_lumi),
                        ],
                        [
                            control_binning[channel][v]
                            for v in set(control_binning[channel].keys())
                            & set(variable_list)
                        ],
                    )]
        return data_dict

    # Step 1: create units and book actions
    for channel in args.channels:
        control_binning = get_control_binning(channel, variable_list)
        print(control_binning)
        nominals[args.era]["datasets"][channel] = get_nominal_datasets(
            args.era, channel
        )
        if args.control_plots:
            nominals[args.era]["units"][channel] = get_control_units(
                channel, args.era, nominals[args.era]["datasets"][channel], control_binning
            )
    um = UnitManager()

    if args.process_selection is None:
        if runPlot:
            run_nums = set()
            for run_key in nominals[args.era]["units"][channel].keys():
                run_nums.add(run_key)

        if runPlot:
            procS = {
                *run_nums,
                "dy",
                "w",
                "tt",
                "st",
                "vv",
                "ewktau",
            }
        else:
            procS = {
                "data",
                "dy",
                "w",
                "tt",
                "st",
                "vv",
                "ewktau",
            }

    else:
        procS = args.process_selection

    print("Processes to be computed: ", procS)

    simulatedProcsDS = {
        "mm": {
            "dy",
            "w",
            "tt",
            "st",
            "vv",
            "ewktau",
        },
        "ee": {
            "dy",
            "w",
            "tt",
            "st",
            "vv",
            "ewktau",
        },
        "mmet": {
            "dy",
            "w",
            "tt",
            "st",
            "vv",
            "ewktau",
        },
        "emet": {
            "dy",
            "w",
            "tt",
            "st",
            "vv",
            "ewktau",
        }
    }

    dataS = set()
    if runPlot:
        for run_key in nominals[args.era]["units"][channel].keys():
            if run_key not in simulatedProcsDS[channel]:
                dataS.add(run_key)
    else: dataS = {"data"} & procS

    # exports the names of the data plots to a text file
    plot_names = ["{}\n".format(plot_name) for plot_name in dataS]
    plot_names.sort()
    with open(r'data_plot_names.txt', 'w') as fp:
        fp.writelines(plot_names)
        fp.close()

    for ch_ in args.channels:
        um.book(
            [
                unit
                for d in dataS | simulatedProcsDS[ch_]
                for unit in nominals[args.era]["units"][ch_][d]
            ],
            enable_check=args.enable_booking_check,
        )
        if args.skip_systematic_variations:
            pass
        else:
            # Book variations common to all channels
            um.book(
                [
                    unit
                    for d in simulatedProcsDS[ch_]
                    for unit in nominals[args.era]["units"][ch_][d]
                ],
                [
                    *mu_id_weight,
                ],
                enable_check=args.enable_booking_check,
            )

    # Step 2: convert units to graphs and merge them
    g_manager = GraphManager(um.booked_units, True)
    g_manager.optimize(args.optimization_level)
    graphs = g_manager.graphs
    for graph in graphs:
        print("%s" % graph)
    if args.only_create_graphs:
        if args.control_plots:
            graph_file_name = "control_unit_graphs-{}-{}-{}.pkl".format(
                args.era, ",".join(args.channels), ",".join(sorted(procS))
            )
        else:
            graph_file_name = "analysis_unit_graphs-{}-{}-{}-{}.pkl".format(
                args.tag, args.era, ",".join(args.channels), args.proc_arr
            )
        if args.graph_dir is not None:
            graph_file = os.path.join(args.graph_dir, graph_file_name)
        else:
            graph_file = graph_file_name
        logger.info("Writing created graphs to file %s.", graph_file)
        with open(graph_file, "wb") as f:
            pickle.dump(graphs, f)
    else:
        # Step 3: convert to RDataFrame and run the event loop
        print("GRAPHS:", graphs)
        r_manager = RunManager(graphs)
        r_manager.run_locally(output_file, args.num_processes, args.num_threads)
    return

if __name__ == "__main__":
    # from multiprocessing import set_start_method
    # set_start_method("spawn")
    args = parse_arguments()
    seperateVariable = bool(int(args.seperate_variables))
    if seperateVariable: variable_list = seperateVariables(args.control_plot_set)
    else: variable_list = args.control_plot_set
    if ".root" in args.output_file:
        log_file = args.output_file.replace(".root", ".log")
    else:
        log_file = "{}.log".format(args.output_file)
    setup_logging(log_file, logging.DEBUG)
    main(args)
