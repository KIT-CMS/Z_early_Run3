# !/bin/sh

# channel options include (mm/ee/mmet/emet)
channel=$1

# output postfix
postfix="DYWNLO"
# postfix="DYWLO"

# if 0 then it sums the runs, if 1 you get a run by run plot
runPlot=1

# if this value is not 0 it will normalize the data to the MC
matchData=0

# if this value is  1 it will take seperate the Variables just like in the histogramming code
seperateVariables=0

# if you summed the runs change this to match the total lumi for the label (in fb-1)
lumiLabel="656.885772534"

# if value is not 0 this will write all the plots you make to a latex file and create a slideshow with them
# Always adds on so you must go to (presentation/control-plots-slides.tex) if you want to delete/change the presentation
latex=0

# sets up a latex environment (thanks Marco)
export LATEX=/ceph/mlink/texlive/2022
export PATH=$LATEX/bin/x86_64-linux:$PATH
export MANPATH=$LATEX/texmf-dist/doc/man:$MANPATH
export INFOPATH=$LATEX/texmf-dist/doc/info:$INFOPATH

# variables by channel
if [ $channel == 'mm' ]
then
   variable_list='pt_1,eta_1,phi_1,m_vis,met_uncorrected,metphi_uncorrected'
   # variable_list='pt_1,pt_2,eta_1,eta_2,phi_1,phi_2,m_vis,pt_vis,mt_uncorrected,met_uncorrected,metphi_uncorrected,pfmet_uncorrected,pfmetphi_uncorrected,nTrackerLayers_1,nTrackerLayers_2,nStations_1,nStations_2'

elif [ $channel == 'mmet' ]
then
   variable_list='pt_1,eta_1,phi_1,met_uncorrected,metphi_uncorrected'
   # variable_list='pt_1,eta_1,phi_1,mt_uncorrected,met_uncorrected,metphi_uncorrected,pfmet_uncorrected,pfmetphi_uncorrected,nTrackerLayers_1,nStations_1'

elif [ $channel == 'ee' ]
then
   variable_list='pt_1,eta_1,phi_1,m_vis,met_uncorrected,metphi_uncorrected'
   # variable_list='pt_1,pt_2,eta_1,eta_2,etaSC_1,etaSC_2,phi_1,phi_2,m_vis,pt_vis,mt_uncorrected,met_uncorrected,metphi_uncorrected,pfmet_uncorrected,pfmetphi_uncorrected,deltaetaSC_1,deltaetaSC_2,eInvMinusPInv_1,eInvMinusPInv_2,hoe_1,hoe_2,scEtOverPt_1,scEtOverPt_2,sieie_1,sieie_2,lostHits_1,lostHits_2'
   
elif [ $channel == 'emet' ]
then
   variable_list='pt_1,eta_1,phi_1,met_uncorrected,metphi_uncorrected'
   # variable_list='pt_1,eta_1,etaSC_1,phi_1,mt_uncorrected,met_uncorrected,metphi_uncorrected,pfmet_uncorrected,pfmetphi_uncorrected,deltaetaSC_1,eInvMinusPInv_1,hoe_1,scEtOverPt_1,sieie_1,lostHits_1'
fi

# creates the plots
source plotting/plot_shapes_control.sh \
2022 \
output/earlyRun3_2022_"$channel"_runPlot"$runPlot"_"$postfix".root \
$variable_list \
None \
$matchData \
$seperateVariables \
$lumiLabel $latex $channel earlyRun3_2022_"$channel"_runPlot"$runPlot"_"$postfix"_matchData"$matchData"

# compiles control-plots-slides.tex into a pdf
if [ $latex != 0 ]
then
   pdflatex presentation/control-plots-slides.tex
fi
