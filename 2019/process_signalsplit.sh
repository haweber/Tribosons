#!/bin/bash

# trap "kill 0" EXIT

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Help
usage()
{
    echo "ERROR - Usage:"
    echo
    echo "      sh $(basename $0) OPTIONSTRINGS ..."
    echo
    echo "Options with arguments:"
    echo "  -h    Help                   (Display this message)"
    echo "  -i    Input baby version     (e.g. -i WWW2017_v5.0.0)"
    echo "  -t    Job tag                (e.g. -t test1)"
    echo "  -u    Enable user study"
    echo "  -x    Skip cutflow histograms"
    echo "  -k    Do skim"
    echo "  -s    Do systematics"
    echo "  -r    Username for input     (e.g. -r mliu or -r phchang)"
    echo
    exit
}

# Default value
CUTFLOW="-C"

# Command-line opts
while getopts ":i:t:r:uxksh" OPTION; do
  case $OPTION in
    i) INPUT_BABY_VERSION=${OPTARG};;
    t) JOB_TAG=${OPTARG};;
    r) USERNAME=${OPTARG};;
    u) DO_USER_STUDY=" --user_study ";;
    x) CUTFLOW=" ";;
    k) DOSKIM=" -K";;
    s) SYSTEMATICS=" -S";;
    h) usage;;
    :) usage;;
  esac
done

if [ -z ${INPUT_BABY_VERSION} ]; then usage; fi
if [ -z ${JOB_TAG}  ]; then usage; fi
if [ -z ${USERNAME} ]; then USERNAME=mliu; fi
if [ -z ${DOSKIM} ]; then JOBS=" -L -H"; fi

# to shift away the parsed options
shift $(($OPTIND - 1))

# Verbose
date
echo "================================================"
echo "$(basename $0) $*"
echo "$(basename $0) $*" >> $DIR/.$(basename $0).history
echo "------------------------------------------------"
echo "INPUT_BABY_VERSION  : ${INPUT_BABY_VERSION}"
echo "JOB_TAG             : ${JOB_TAG}"
echo "================================================"

# DDFAKEDIR=data/
DDFAKEDIR=bkgdata/
#if [[ ${INPUT_BABY_VERSION} == *"2016"* ]]; then
#    DDFAKEDIR=data/
#fi

mkdir -p hists/${INPUT_BABY_VERSION}/${JOB_TAG}/

if [[ ${USERNAME} == *"mliu"* ]]; then
    SUBDIRECTORY=skim
else
    SUBDIRECTORY=""
fi

echo "Submitting parallel jobs ... ==>"
echo "(below are individual bash commands)"
# Split a couple of big jobs by a few sub jobs
(set -x ;./doAnalysis ${DOSKIM} ${SYSTEMATICS} ${CUTFLOW} ${JOBS} ${DO_USER_STUDY} -i "/nfs-7/userdata/${USERNAME}/WWW_babies/${INPUT_BABY_VERSION}/${SUBDIRECTORY}/grouped/vvv/www_*.root"    -o hists/${INPUT_BABY_VERSION}/${JOB_TAG}/www_onshell.root            -T t                   > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/www_onshell.log            2>&1) & 
(set -x ;./doAnalysis ${DOSKIM} ${SYSTEMATICS} ${CUTFLOW} ${JOBS} ${DO_USER_STUDY} -i "/nfs-7/userdata/${USERNAME}/WWW_babies/${INPUT_BABY_VERSION}/${SUBDIRECTORY}/grouped/vvv/wh_ww_*.root"    -o hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wh_www.root            -T t                   > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wh_www.log            2>&1) &
(set -x ;./doAnalysis ${DOSKIM} ${SYSTEMATICS} ${CUTFLOW} ${JOBS} ${DO_USER_STUDY} -i "/nfs-7/userdata/${USERNAME}/WWW_babies/${INPUT_BABY_VERSION}/${SUBDIRECTORY}/grouped/vvv/wwz_*.root"    -o hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wwz_onshell.root            -T t                   > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wwz_onshell.log            2>&1) & 
(set -x ;./doAnalysis ${DOSKIM} ${SYSTEMATICS} ${CUTFLOW} ${JOBS} ${DO_USER_STUDY} -i "/nfs-7/userdata/${USERNAME}/WWW_babies/${INPUT_BABY_VERSION}/${SUBDIRECTORY}/grouped/vvv/zh_ww_*.root"    -o hists/${INPUT_BABY_VERSION}/${JOB_TAG}/zh_wwz.root            -T t                   > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/zh_wwz.log            2>&1) &
(set -x ;./doAnalysis ${DOSKIM} ${SYSTEMATICS} ${CUTFLOW} ${JOBS} ${DO_USER_STUDY} -i "/nfs-7/userdata/${USERNAME}/WWW_babies/${INPUT_BABY_VERSION}/${SUBDIRECTORY}/grouped/vvv/wzz_*.root"    -o hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wzz_onshell.root            -T t                   > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wzz_onshell.log            2>&1) & 
(set -x ;./doAnalysis ${DOSKIM} ${SYSTEMATICS} ${CUTFLOW} ${JOBS} ${DO_USER_STUDY} -i "/nfs-7/userdata/${USERNAME}/WWW_babies/${INPUT_BABY_VERSION}/${SUBDIRECTORY}/grouped/vvv/wh_zz_*.root"    -o hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wh_wzz.root            -T t                   > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wh_wzz.log            2>&1) & 
(set -x ;./doAnalysis ${DOSKIM} ${SYSTEMATICS} ${CUTFLOW} ${JOBS} ${DO_USER_STUDY} -i "/nfs-7/userdata/${USERNAME}/WWW_babies/${INPUT_BABY_VERSION}/${SUBDIRECTORY}/grouped/vvv/zzz_*.root"    -o hists/${INPUT_BABY_VERSION}/${JOB_TAG}/zzz_onshell.root            -T t                   > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/zzz_onshell.log            2>&1) & 
(set -x ;./doAnalysis ${DOSKIM} ${SYSTEMATICS} ${CUTFLOW} ${JOBS} ${DO_USER_STUDY} -i "/nfs-7/userdata/${USERNAME}/WWW_babies/${INPUT_BABY_VERSION}/${SUBDIRECTORY}/grouped/vvv/zh_zz_*.root"    -o hists/${INPUT_BABY_VERSION}/${JOB_TAG}/zh_zzz.root            -T t                   > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/zh_zzz.log            2>&1) & 
sleep 1;
echo "<== Submitted parallel jobs ..."
wait

echo ""
echo "Jobs finished ..."
echo ""
echo "Printing how long it took for each job."

tail -n 3 hists/${INPUT_BABY_VERSION}/${JOB_TAG}/*.log

# Hadd the split jobs result
if [ -z ${DOSKIM} ]; then
    echo "Hadding some histogram outputs ... ==>"
    (set -x ; hadd hists/${INPUT_BABY_VERSION}/${JOB_TAG}/vvv_onshell.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/www_onshell.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wwz_onshell.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wzz_onshell.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/zzz_onshell.root > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/vvv_onshell.log) &
    (set -x ; hadd hists/${INPUT_BABY_VERSION}/${JOB_TAG}/vvv_higgs.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wh_www.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wh_wzz.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/zh_wwz.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/zh_zzz.root > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/vvv_higgs.log) &
    (set -x ; hadd hists/${INPUT_BABY_VERSION}/${JOB_TAG}/vvv_all.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/www_onshell.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wwz_onshell.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wzz_onshell.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/zzz_onshell.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wh_www.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wh_wzz.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/zh_wwz.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/zh_zzz.root > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/vvv_all.log) &
    (set -x ; hadd hists/${INPUT_BABY_VERSION}/${JOB_TAG}/www_all.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/www_onshell.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wh_www.root > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/www_all.log) &
    (set -x ; hadd hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wwz_all.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wwz_onshell.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/zh_wwz.root > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/www_all.log) &
    (set -x ; hadd hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wzz_all.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wzz_onshell.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/wh_wzz.root > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/www_all.log) &
    (set -x ; hadd hists/${INPUT_BABY_VERSION}/${JOB_TAG}/zzz_all.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/zzz_onshell.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/zh_zzz.root > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/www_all.log) &
fi
wait
echo "<== Done hadding histogram outputs!"

echo "histograms are at hists/${INPUT_BABY_VERSION}/${JOB_TAG}/*.root"
echo "You can find a print out of how long each subjobs took above for reference"
