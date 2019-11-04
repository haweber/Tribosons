#!/bin/bash

trap "kill 0" EXIT

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
#(set -x ;./doAnalysis ${DOSKIM} ${SYSTEMATICS} ${CUTFLOW} ${JOBS} ${DO_USER_STUDY} -i /nfs-7/userdata/${USERNAME}/WWW_babies/${INPUT_BABY_VERSION}/${SUBDIRECTORY}/grouped/bkg/          -o hists/${INPUT_BABY_VERSION}/${JOB_TAG}/fakes.root          -T t_fakes             > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/fakes.log          2>&1) &
(set -x ;./doAnalysis ${DOSKIM} ${SYSTEMATICS} ${CUTFLOW} ${JOBS} ${DO_USER_STUDY} -i "/nfs-7/userdata/${USERNAME}/WWW_babies/${INPUT_BABY_VERSION}/${SUBDIRECTORY}/grouped/sigofficial/www*root"          -o hists/${INPUT_BABY_VERSION}/${JOB_TAG}/www.root -T t_www               > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/signal_private.log 2>&1) & 
(set -x ;./doAnalysis ${DOSKIM} ${SYSTEMATICS} ${CUTFLOW} ${JOBS} ${DO_USER_STUDY} -i "/nfs-7/userdata/${USERNAME}/WWW_babies/${INPUT_BABY_VERSION}/${SUBDIRECTORY}/grouped/sigofficial/vh*.root"  -o hists/${INPUT_BABY_VERSION}/${JOB_TAG}/whtowww.root         -T t_www               > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/signal.log         2>&1) & # Official CMS sample
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
    (set -x ; hadd hists/${INPUT_BABY_VERSION}/${JOB_TAG}/lostlep.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/lostlep_*.root > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/lostlep.log) &
    (set -x ; hadd hists/${INPUT_BABY_VERSION}/${JOB_TAG}/ddfakes.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/ddfakes_*.root > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/ddfakes.log) &
    (set -x ; hadd hists/${INPUT_BABY_VERSION}/${JOB_TAG}/fitlostlep.root hists/${INPUT_BABY_VERSION}/${JOB_TAG}/fitlostlep_*.root > hists/${INPUT_BABY_VERSION}/${JOB_TAG}/fitlostlep.log) &
fi
wait
echo "<== Done hadding histogram outputs!"

echo "histograms are at hists/${INPUT_BABY_VERSION}/${JOB_TAG}/*.root"
echo "You can find a print out of how long each subjobs took above for reference"
