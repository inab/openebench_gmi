#!/bin/bash

scriptdir="$(dirname "$0")"
case "$scriptdir" in
	/*)
		true
		;;
	.)
		scriptdir="$(pwd)"
		;;
	*)
		scriptdir="$(pwd)/${scriptdir}"
		;;
esac

set -e

# Getting the workflow version from the nextflow.config
if [ $# = 0 ] ; then
	version_tag="$(grep -o "^ *version = '[^']*'" "${scriptdir}"/nextflow.config | cut -f 2 -d "'")"
else
	version_tag="$1"
fi


declare -A name2tag=(
	[checkFormat]=openebench_gmi/sample-checkresults
	[getQueryIds]=openebench_gmi/sample-getqueryids
	[compareIds]=openebench_gmi/sample-compareids
	[robinsonFouldsMetric]=openebench_gmi/sample-robinsonfoulds
	[calculateSnPrecision]=openebench_gmi/sample-calculatesnprecision
	[assessmentSnPrecision]=openebench_gmi/sample-assessment-snprecision
	[assessmentRfHeatmap]=openebench_gmi/sample-assessment-rfheatmap
)

for A in "${scriptdir}"/containers/* ; do
	name_tag="${name2tag[$(basename "$A")]}"
	if [ -n "${name_tag}" ] ; then
		tag="${name_tag}:${version_tag}"
		echo "[$(date -Is)] Building ${tag} from $A"
		docker build -t "${tag}" -f "${A}/Dockerfile" "${scriptdir}"
	fi
done
