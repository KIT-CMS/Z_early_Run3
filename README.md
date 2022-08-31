
Framework to produce Control Plots

Before cloning the early Run3 analysis framework start the ssh agent:
```bash
eval $(ssh-agent)
ssh-add
```
Àfterwards checkout the early Run3 analysis framework:
```bash
git clone --recursive git@github.com:KIT-CMS/Z_early_Run3.git
```
1. setup the environment
```bash
source utils/setup_root.sh
```
2. Go to 'make-histograms.sh' and make your changes according to your histogramming desires (comments to guide in the script)
    - To define new variables not in the ntuple go to control_binning.py
    - To change the binning/bounds go to control_binning.py
    - To change the cuts go to channel_selection.py
    - To change the MC weights go to process_selection.py
    - To change the data weights (such as lumi scaling) go to data_selection.py
    - To specify new file names for each sample go to file_names.py
    - Everything else should be able to be modified in the bash script but if you want to add new features a good place to start is produce_shapes.py

3. Run the histogramming code
```bash
. make-histograms.sh
```
* If you would prefer to type all of this in the terminal use here is the command:
```bash
   python shapes/produce_shapes.py --channels $channel --output-file output/earlyRun3_crown_2018_"$channel"\
   --directory /ceph/moh/CROWN_samples/EarlyRun3_V12/ntuples\
   --$channel-friend-directory /ceph/moh/CROWN_samples/EarlyRun3_V12/friends/crosssection\
   --era 2018 --num-processes 64 --num-threads 64 --optimization-level 1 --control-plots\
   --run-plot $runPlot --group-runs $groupRuns --seperate-variables $seperateVariable\
   --control-plot-set $variable_list --total-lumi $totalLuminosity --run-list $run_list --lumi-list $lumi_list --ntuple_type crown --skip-systematic-variations
```
Example:
```bash
python shapes/produce_shapes.py --channels mm --output-file output/earlyRun3_crown_2022_mm\
   --directory /ceph/moh/CROWN_samples/EarlyRun3_V11/ntuples\
   --mm-friend-directory /ceph/moh/CROWN_samples/EarlyRun3_V11/friends/crosssection\
   --era 2018 --num-processes 16 --num-threads 16 --optimization-level 1 --control-plots\
   --run-plot 1 --group-runs 5.2 seperate-variables 1\
   --control-plot-set pt_1,pt_2,eta_1,eta_2,m_vis --total-lumi 0.088 --run-list 355443,355444,355445\
   --lumi-list 1.971790835,1.243983634,1.537414347 --ntuple_type crown --skip-systematic-variations
```
- $channel takes in either a list of channels (mm,ee,mmet,emet) and histograms the same variables for all or you can histogram specific variables for each channel by going to make-histograms.sh and reading the comments, then run one channel at a time.
- $runPlot is 1 for a run by run plot and 0 to just sum all the runs
-$seperateVariables if 1 will allow you to seperate the variables by things like number of primary vertices, positive and negative boson etc. In order to choose which things it seperates by you must go to control_binning.py and change the code under seperateVariables accordingly (comments to guide)
- $groupRuns is 0 to not group runs by luminosity and any luminosity value to choose the max for grouping luminosities
- $variable_list takes in a list seperated by comments of variables you want (suggested to go to make-histograms.sh to store big lists of variables)
- $totalLuminosity can be anything if you are doing run by run histograms because it scales everything to 1 pb-1 but for the summed run control plots you should put your total luminosity here
- $run_list can be anything if you are not doing run by run histogramming and if runPlot is true you list your run numbers here seperated by commas
- $lumi_list has the corresponding luminosities to the run numbers (keep them in the same order respectively)

4. Open up a new terminal for the plotting

5. Go to 'make-plots.sh' and make your changes according to your plotting desires (comments to guide in the script)
    - To change your pdf slideshow go to presentation/control-plots-slides.tex
    - To add logos for your pdf slideshow go to /logos
    - To get your chi square data for each plot go to /chi_square_data

6. Run the plotting code
```bash
. make-plots.sh
```

* If you would prefer to type all of this in the terminal use here is the command:
```bash
source plotting/plot_shapes_control.sh 2018 output/earlyRun3_crown_2022_"$channel".root $variable_list None $matchData $seperateVariables $lumiLabel $latex $channel earlyRun3$channel
```
Example:
```bash
source plotting/plot_shapes_control.sh 2018 output/earlyRun3_crown_2022_mm.root pt_1,pt_2,eta_1,eta_2,m_vis None 1 1 0.088 1 mm earlyRun3mm
```
- $channel takes in either a list of channels (mm,ee,mmet,emet) and plots the same variables for all or you can plot specific variables for each channel by going to make-plots.sh and reading the comments, then run one channel at a time.
- $variable_list takes in a list seperated by comments of variables you want (suggested to go to make-plots.sh to store big lists of variables)
- $matchData is 1 to normalize data to MC and 0 to turn this feature off
-$seperateVariables should be the same as in the histogramming, but this time go to the bottom of plot_shapes_control.py to make sure you are seperating by the same variables
- $lumiLabel if what luminosity will be displayed on the plot (match $totalLuminosity or 0 if you are making a run plot
- $latex is 1 to write the plots to a latex file that will make a slideshow pdf and 0 to turn this feature off

If you did not use . make-plots.sh and used the terminal instead you must also run the following:
```bash
export LATEX=/ceph/mlink/texlive/2022
export PATH=$LATEX/bin/x86_64-linux:$PATH
export MANPATH=$LATEX/texmf-dist/doc/man:$MANPATH
export INFOPATH=$LATEX/texmf-dist/doc/info:$INFOPATH
pdflatex presentation/control-plots-slides.tex
```
