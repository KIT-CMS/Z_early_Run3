#!/bin/bash



export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
export VO_CMS_SW_SETUPFILE=cmsset_default.sh
CVMFS_CMS_CHECK_CODE=`timeout 5 ls $VO_CMS_SW_DIR &> /dev/null; echo $?`

case "$CVMFS_CMS_CHECK_CODE" in
    0)
	if [ -f $VO_CMS_SW_DIR/$VO_CMS_SW_SETUPFILE ]; then
	    source $VO_CMS_SW_DIR/$VO_CMS_SW_SETUPFILE
	    cd /work/moh/EarlyRun3/SignalExtraction/CMSSW_10_6_0/src
	    cmsenv
	    cd -
	else
	    echo "folder '$VO_CMS_SW_DIR' exists, but '$VO_CMS_SW_SETUPFILE' is not found!?"
	fi
	;;
    2)
	echo "---> cvmfs folder missing"
	;;
    124)
	echo "---> timeout for cvmfs mount"
	;;
    *)
	echo "---> unknown cvmfs problem"
	;;
esac
