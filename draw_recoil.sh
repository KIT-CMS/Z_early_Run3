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
#totalLuminosity=4.919548635
# totalLuminosity=5.054203038
totalLuminosity=5.035650234254
# Run numbers and correspoding luminosities you want to plot
#run_list="355862,355863,355870,355871,355872,355892,355912,355913,355921,355933,355942,355988,355989,355998,355999,356004,356005,356043,356071,356074,356075,356076,356077,356135,356309,356316,356322,356323,356371,356375,356378,356381,356383,356385,356386,356426,356428,356433,356434,356435,356446,356523,356531,356563,356568,356569,356570,356576,356578,356580,356582,356614,356615,356619,356810,356811,356812,356813,356814,356815,356824,356908,356919,356937,356946,356947,356948,356949,356951,356954,356955,356956,356968,356969,356970,356998,356999,357000,357001,357079,357080,357081,357101,357102,357104,357106,357112,357268,357271,357328,357329,357330,357331,357332,357333,357401,357406,357438,357440,357441,357442,357447,357472,357478,357479,357482" #,357538,357542,357550,357610,357611,357612,357613,357688,357696,357697,357698,357699,357700,357720,357732,357734,357735,357754,357756,357757,357758,357759,357766,357777,357778,357779,357781,357802,357803,357804,357805,357806,357807,357808,357809,357812,357813,357814,357815,357898,357899,357900"
run_list='355862,355863,355870,355871,355872,355892,355912,355913,355921,355933,355942,355988,355989,355998,355999,356004,356005,356043,356071,356074,356075,356076,356077,356135,356309,356316,356322,356323,356371,356375,356378,356381,356383,356385,356386,356426,356428,356433,356434,356435,356446,356523,356531,356563,356568,356569,356570,356576,356578,356580,356582,356614,356615,356619,356810,356811,356812,356813,356814,356815,356824,356908,356919,356937,356946,356947,356948,356949,356951,356954,356955,356956,356968,356969,356970,356998,356999,357000,357001,357079,357080,357081,357101,357102,357104,357106,357112,357268,357271,357328,357329,357330,357331,357332,357333,357401,357406,357438,357440,357441,357442,357447,357472,357478,357479,357482'
#lumi_list="0.622396724,0.732968640,1.561198678,0.146692485,21.029407135,1.976640736,9.007156928,5.865767709,22.873241336,0.724271840,1.371972162,2.780268121,1.265671750,5.567402786,15.248489527,2.999361674,13.902884694,17.133845292,49.787188261,2.612524012,10.957679862,9.292921771,7.684875420,2.060615508,67.033187162,1.569606629,6.515209807,41.464490983,131.203129037,2.867059944,2.467573112,8.706233768,2.563759568,45.822680251,47.969426483,1.724343482,0.441028997,92.536256058,159.668345814,9.002263709,51.798236133,3.193804083,45.207670585,19.147074536,25.929216259,154.074111931,7.723225698,14.350816146,7.836097276,229.874631464,23.957143767,12.758405557,29.774954537,14.427146840,155.063334801,19.045058546,9.448816473,7.223515810,20.963191799,27.999009333,34.170473236,89.708208698,19.762953445,21.181707938,58.164494098,67.344980476,60.931797755,14.862427888,46.499503522,65.343931988,93.623773165,1.536626252,16.473913043,10.743512440,50.102746465,7.271353976,188.379172811,174.168965073,16.780629770,33.998823482,1.227174925,16.692756409,131.700550888,20.398303272,405.808885633,15.868432136,209.375038734,40.740482691,5.338217312,93.824183818,41.319544454,161.737018459,32.380550655,60.683221298,112.206319937,24.488058518,307.175729738,2.360207763,7.243290668,1.576219850,323.123255284,49.434988240,5.338575692,81.695935046,11.679759817,44.545306255,91.469413147,153.895784138,45.152973037,128.258466980,104.182680415,12.552330421,18.430599972,9.129904784,155.235799381,17.516143488,40.211495629,107.254510168,315.110543965,32.728813604,145.822389695,2.085368937,24.880551721,20.209734650,35.896329417,22.849018795,85.681800505,16.547938028,2.080793211,51.652435163,51.556078769,7.542239346,29.014013547,15.299643768,76.420645302,4.233347072,11.826321268,14.290667592,77.920791535,50.094097961,191.247913699,102.645748876,176.490809327,109.830466091"
# lumi_list="0.63801,0.751527,1.593551,0.151426,70.984318,10.969796,9.169095,5.951526,23.231373,0.723036,1.391794,2.834268,1.291898,2.69918,0.658442,1.383661,12.394444,5.69176,15.582188,3.060384,14.177479,17.467281,50.71093,2.673444,11.150236,9.461482,2.102951,73.414434,1.58137,6.659355,42.348427,133.773893,2.920617,2.513286,8.86264,2.616358,46.652603,48.684804,1.747974,0.446594,94.481816,163.787204,9.130886,52.783124,3.255192,46.102188,19.429062,26.464137,157.21535,7.87151,14.514699,8.022597,234.972763,24.406282,25.888664,13.089688,30.513522,14.781274,158.691482,19.448532,9.63696,7.249653,21.521138,28.779681,35.052868,92.070014,20.269006,21.729444,59.684524,69.065375,64.33714,15.206617,47.913384,67.294466,96.40275,1.589575,17.031531,11.087986,51.66541,7.524719,194.391898,179.260552,17.409691,35.218879,1.27001,17.278594,136.203862,20.971251,415.579157,16.329014,216.355977,42.300421,5.541778,97.553233,42.867927,166.683657,33.18762,61.491152,113.144177,24.615518,308.296505,2.394476,7.149332,1.59943,326.429059,49.599818"
lumi_list="0.624516690,0.735540770,1.559658372,0.148169025,69.515629423,10.802601699,8.978575546,5.825279127,22.754448690,0.692376260,1.330926526,2.775116204,1.266494092,2.644759427,0.645082250,1.354422830,12.129229856,5.580313227,15.284372935,3.002265913,13.907312666,17.135853244,49.724475667,2.621015734,10.927730525,9.260169736,2.063053318,72.008391986,1.551548623,6.537815457,41.621399315,131.202161909,2.860354248,2.461457717,8.674933521,2.567762692,45.830078472,47.817308786,1.716021590,0.438410753,92.813938929,161.067921139,8.969844305,51.906559873,3.199248448,45.345399980,19.114854891,25.989430035,154.545768009,7.730211736,14.249111553,7.888687221,231.035636661,23.959490105,25.463143250,12.899831904,30.056960688,14.554483681,156.058618868,19.101244867,9.458727980,7.144502553,21.189077486,28.321158564,34.555747451,90.680232790,19.946081920,21.375902260,58.684501624,67.852169890,63.160285620,14.918061732,47.233298410,66.362559750,94.953445066,1.569041624,16.806636519,10.927795917,50.926111348,7.431664161,191.767464344,176.471880331,17.190981921,34.765081149,1.253125974,17.043522406,134.226009513,21.328085420,422.441519318,16.607171051,220.289160414,43.014505895,5.633690874,99.115875716,43.522343403,169.436618241,33.744821712,62.649565245,115.231693892,25.047635701,313.175603407,2.433494283,7.269330768,1.626275296,331.908025194,50.432329154"

# if this value is not 0 it will normalize the data to the MC
matchData=0

# if 1 it will seperate variables by things you can define in control_binning.py (useful for QCD template)
seperateVariables=0

extra_args=""
# extra_args="--doLepCorrBins"
# extra_args="--doZpt"
#extra_args="--doQCD"

# variables by channel (this can be changed to fit what variables you want plotted)
# variable_list=""

# for ntuple_dir in /ceph/moh/CROWN_samples/Run3V03/ntuples_xsec_sf_lep_corr_final_v07*
# do
#     postfix=$(basename $ntuple_dir)
#     postfix=${postfix/ntuples_xsec_sf_/"_"}

    # ntuple_dir="/work/jdriesch/CROWN_samples/Run3V04/ntuples_xsec_sf_scaleres_EraC"
    ntuple_dir="/ceph/jdriesch/CROWN_samples/Run3V07/ntuples"
    # postfix="_NoSF"
    # postfix="_SF"
    # postfix="_met_corr"
    # postfix="_Zcounting"  # remove Z_mass_window cut from channel_selection, remove sf_trg from sfWeight
    # postfix="_toFit"
    # postfix="_lepuncorr"
    postfix="_RecoilSyst"
    # postfix="_ZptCorr"
    # postfix=""

    # set subera
    subera="C"
    #lumi_label=4.844307925632
    #lumi_label=4.919548635
    # lumi_label=5.054203038
    lumi_label=5.035650234254
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
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/lepton \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/xsec \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/pu \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/xy \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/sf \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/met \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/met_double \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/met_sigOnly \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/met_zrap0 \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/met_zrap1 \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/met_zrap2 \
        --mmet-friend-directory \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/lepton \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/xsec \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/pu \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/xy \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/sf \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/met \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/met_double \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/met_sigOnly \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/met_zrap0 \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/met_zrap1 \
            /ceph/jdriesch/CROWN_samples/Run3V07/friends/met_zrap2 \

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
