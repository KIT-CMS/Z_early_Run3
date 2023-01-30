#!/bin/sh

# shape or plot
type=$1

# enter which channel you want (mm/ee/mmet/emet)
channel=$2

nthreads=8
if [ ! -z "$3" ]
then
    nthreads=$3
fi

# if 0 then it sums the runs, if 1 you get a run by run plot
# runPlot=0
runPlot=1

# input the max luminosity you want runs to be grouped by (if 0 it will not group runs)
# groupRuns=0
groupRuns=1000

# luminosity of the summed runs in fb-1 (if you use runPlot it will not matter what number it here because it will be weighted to 1 pb-1)
# totalLuminosity=7.544816105842
totalLuminosity=4.919548635
# Run numbers and correspoding luminosities you want to plot
run_list="355862,355863,355870,355871,355872,355892,355912,355913,355921,355933,355942,355988,355989,356043,356071,356074,356075,356076,356077,356135,356309,356316,356321,356322,356323,356371,356375,356378,356381,356383,356385,356386,356426,356428,356433,356434,356435,356446,356523,356531,356563,356568,356569,356570,356576,356578,356580,356582,356614,356615,356619,356811,356812,356813,356814,356815,356824,356908,356919,356937,356946,356947,356948,356949,356951,356954,356955,356956,356968,356969,356970,356998,356999,357000,357001,357079,357080,357081,357101,357102,357104,357106,357112,357268,357271,357328,357329,357330,357331,357332,357333,357401,357406,357438,357440,357441,357442,357447,357472,357478,357479,357482" #,357538,357542,357550,357610,357611,357612,357613,357688,357696,357697,357698,357699,357700,357720,357732,357734,357735,357754,357756,357757,357758,357759,357766,357777,357778,357779,357781,357802,357803,357804,357805,357806,357807,357808,357809,357812,357813,357814,357815,357898,357899,357900"
lumi_list="0.622396724,0.732968640,1.561198678,0.146692485,21.029407135,1.976640736,9.007156928,5.865767709,22.873241336,0.724271840,1.371972162,2.780268121,1.265671750,5.567402786,15.248489527,2.999361674,13.902884694,17.133845292,49.787188261,2.612524012,10.957679862,9.292921771,7.684875420,2.060615508,67.033187162,1.569606629,6.515209807,41.464490983,131.203129037,2.867059944,2.467573112,8.706233768,2.563759568,45.822680251,47.969426483,1.724343482,0.441028997,92.536256058,159.668345814,9.002263709,51.798236133,3.193804083,45.207670585,19.147074536,25.929216259,154.074111931,7.723225698,14.350816146,7.836097276,229.874631464,23.957143767,12.758405557,29.774954537,14.427146840,155.063334801,19.045058546,9.448816473,7.223515810,20.963191799,27.999009333,34.170473236,89.708208698,19.762953445,21.181707938,58.164494098,67.344980476,60.931797755,14.862427888,46.499503522,65.343931988,93.623773165,1.536626252,16.473913043,10.743512440,50.102746465,7.271353976,188.379172811,174.168965073,16.780629770,33.998823482,1.227174925,16.692756409,131.700550888,20.398303272,405.808885633,15.868432136,209.375038734,40.740482691,5.338217312,93.824183818,41.319544454,161.737018459,32.380550655,60.683221298,112.206319937,24.488058518,307.175729738,2.360207763,7.243290668,1.576219850,323.123255284,49.434988240,5.338575692,81.695935046,11.679759817,44.545306255,91.469413147,153.895784138,45.152973037,128.258466980,104.182680415,12.552330421,18.430599972,9.129904784,155.235799381,17.516143488,40.211495629,107.254510168,315.110543965,32.728813604,145.822389695,2.085368937,24.880551721,20.209734650,35.896329417,22.849018795,85.681800505,16.547938028,2.080793211,51.652435163,51.556078769,7.542239346,29.014013547,15.299643768,76.420645302,4.233347072,11.826321268,14.290667592,77.920791535,50.094097961,191.247913699,102.645748876,176.490809327,109.830466091"

# if this value is not 0 it will normalize the data to the MC
matchData=0

# if 1 it will seperate variables by things you can define in control_binning.py (useful for QCD template)
seperateVariables=0

extra_args=""
# extra_args="--doLepCorrBins"
# extra_args="--doZpt"
# extra_args="--doQCD"

# variables by channel (this can be changed to fit what variables you want plotted)
# variable_list=""

# for ntuple_dir in /ceph/moh/CROWN_samples/Run3V03/ntuples_xsec_sf_lep_corr_final_v07*
# do
#     postfix=$(basename $ntuple_dir)
#     postfix=${postfix/ntuples_xsec_sf_/"_"}

    ntuple_dir="/ceph/moh/CROWN_samples/Run3V03/ntuples_xsec_sf_EraC_lep_corr_01_x0p60"
    # postfix="_NoSF"
    # postfix="_SF"
    # postfix="_met_corr"
    # postfix="_Zcounting"  # remove Z_mass_window cut from channel_selection, remove sf_trg from sfWeight
    # postfix="_toFit"
    # postfix="_lepuncorr"
    # postfix="_RecoilSyst"
    # postfix="_ZptCorr"
    postfix=""

    # set subera
    subera="C"
    #lumi_label=4.844307925632
    lumi_label=4.919548635
    # if [[ "$ntuple_dir" == *"_C"* ]]; then
    #     subera="C"
    #     lumi_label=4.844307925632
    # elif [[ "$ntuple_dir" == *"_D1"* ]]; then
    #     subera="D1"
    #     lumi_label=0.9192951682019999
    # elif [[ "$ntuple_dir" == *"_D2"* ]]; then
    #     subera="D2"
    #     lumi_label=1.7812130120119998
    # else
    #     subera="Run2022"
    #     lumi_label=$totalLuminosity
    # fi

    if [[ "$extra_args" == *"--doLepCorrBins"* ]]
    then
        postfix=$postfix'_LepCorrBins'
        seperateVariables=1
    fi
    if [[ "$extra_args" == *"--doZpt"* ]]
    then
        if [[ "$postfix" == *"_ZptCorr"* ]]
        then
            postfix=$postfix
        else
            postfix=$postfix'_Zpt'
        fi
        seperateVariables=1
    fi
    if [[ "$extra_args" == *"--doQCD"* ]]
    then
        postfix=$postfix'_QCD'
        seperateVariables=1
    fi

    echo ""
    echo ">>> Running $postfix"
    echo ">>>     ntuple_dir: $ntuple_dir"
    echo ">>>     postfix:    $postfix"
    echo ">>>     subera:     $subera"
    echo ">>>     extra_args: $extra_args"

    if [ $type == "shape" ]
    then
        if [[ ! -z "$postfix" ]]
        then
            extra_args="$extra_args --skip-systematic-variations"
        fi
        echo ">>>     extra_args: $extra_args"

        echo ">>> Producing shapes"
        python shapes/produce_shapes.py \
        --channels $channel --era 2022 \
        --output-file output/earlyRun3_2022$postfix \
        --directory $ntuple_dir \
        --num-processes 8 --num-threads $nthreads \
        --optimization-level 1 --control-plots \
        --run-plot $runPlot --group-runs $groupRuns \
        --run-list $run_list --lumi-list $lumi_list \
        --seperate-variables $seperateVariables \
        --total-lumi $totalLuminosity \
        $extra_args \
        --ntuple_type crown \
        --subera $subera \
        --applySF \
        --mm-friend-directory \
            /ceph/moh/CROWN_samples/Run3V03/friend_xsec_sf_EraC_lep_corr_01_x0p60_met_corr \
            /ceph/moh/CROWN_samples/Run3V03/friend_xsec_sf_EraC_lep_corr_01_x0p60_met_corr_double \
            /ceph/moh/CROWN_samples/Run3V03/friend_xsec_sf_EraC_lep_corr_01_x0p60_met_corr_sigOnly \
            /ceph/moh/CROWN_samples/Run3V03/friend_xsec_sf_EraC_lep_corr_01_x0p60_met_corr_zrap0 \
            /ceph/moh/CROWN_samples/Run3V03/friend_xsec_sf_EraC_lep_corr_01_x0p60_met_corr_zrap1 \
            /ceph/moh/CROWN_samples/Run3V03/friend_xsec_sf_EraC_lep_corr_01_x0p60_met_corr_zrap2 \
        --mmet-friend-directory \
            /ceph/moh/CROWN_samples/Run3V03/friend_xsec_sf_EraC_lep_corr_01_x0p60_met_corr \
            /ceph/moh/CROWN_samples/Run3V03/friend_xsec_sf_EraC_lep_corr_01_x0p60_met_corr_double \
            /ceph/moh/CROWN_samples/Run3V03/friend_xsec_sf_EraC_lep_corr_01_x0p60_met_corr_sigOnly \
            /ceph/moh/CROWN_samples/Run3V03/friend_xsec_sf_EraC_lep_corr_01_x0p60_met_corr_zrap0 \
            /ceph/moh/CROWN_samples/Run3V03/friend_xsec_sf_EraC_lep_corr_01_x0p60_met_corr_zrap1 \
            /ceph/moh/CROWN_samples/Run3V03/friend_xsec_sf_EraC_lep_corr_01_x0p60_met_corr_zrap2

    elif [ $type == "plot" ]
    then
        # for syst in Nominal SFTrkUp SFTrkDn SFStaUp SFStaDn SFIDUp SFIDDn SFIsoUp SFIsoDn SFTrgUp SFTrgDn
        for syst in Nominal
        do
            echo ">>> Producing plots"
            if [[ ! -f output/earlyRun3_2022$postfix.root ]]
            then
                echo ">>> output/earlyRun3_2022$postfix.root does not exist --> skip"
                continue
            fi
            python plotting/plot_shapes_control.py -l \
            --era Run2022 \
            --input output/earlyRun3_2022$postfix.root \
            --category-postfix None \
            --match-data $matchData \
            --seperate-variables $seperateVariables \
            --lumi-label $lumi_label \
            --channels $channel \
            $extra_args \
            --tag earlyRun3_2022$postfix \
            --subera $subera \
            --plot-postfix "$postfix" \
            --syst $syst
        done
    fi
# done

echo ""
echo ">>> finished"
echo ""
