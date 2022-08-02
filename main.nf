#!/usr/bin/env nextflow

/*
===============================================================================================================
                  O P E N E B E N C H  W O R K F L O W   O U T B R E A K  D E T E C T I O N
===============================================================================================================
 #### Homepage / Documentation
 https://github.com/BU-ISCIII/openebench_gmi
 https://github.com/javi-gv94/openebench_gmi
 https://github.com/inab/openebench_gmi
 @#### Authors
 Sara Monzon <smonzon@isciii.es>
 Javier Garrayo-Ventas <https://orcid.org/0000-0003-0015-1573>
 José Mª Fernández <https://orcid.org/0000-0002-4806-5140>
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Pipeline overview:
 - 1. : Validation & Check Results
    - 1.1. : GetIDsFromNewick
    - 1.2. : CheckNewickFormat
 - 2. : Metrics
 	- 2.1: Robin-Foulds calculation between participant result and golden dataset.
 	- 2.2: Metrics consolidation.
 - 3. : Visualization

 ----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
     BU-ISCIII/openebench_gmi : OpenEBench pipeline for Outbreak detection challenge v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run BU-ISCIII/openebench_gmi --input {test.newick.file} --goldstandard_dir {golden.folder.path} --assess_dir {assessment.path} --public_ref_dir {path.to.info.ref.dataset} --challenges_ids {event.id}

    Mandatory arguments:
      --input                   Path to input data (must be surrounded with quotes).
      --goldstandard_dir            Path to reference data. Golden datasets.
      --public_ref_dir				Path where public dataset info is stored for validation.
      --assess_dir					Path where benchmark data is stored.
      --challenges_ids                    Event identifier.
      --participant_id				Participant identifier.
      --tree_format					Format tree ["nexus","newick"].

    Other options:
	  --statsdir					The output directory with nextflow statistics
      --outdir                      The output directory where the results will be saved
	  --otherdir					The output directory where custom results will be saved (no directory inside)
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = workflow.manifest.version ?: '1.0'

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

/*
* DEFAULT AND CUSTOM VALUE FOR CONFIGURABLE VARIABLES
*/

if(params.input){
	input_file = file(params.input)
	if (!input_file.exists()) exit 1, "Input Newick file not found: ${params.input}"
}

if(params.goldstandard_dir){
	Channel
		.fromPath( params.goldstandard_dir, type: 'dir',checkIfExists:true )
		.into { goldstandard_dir_robinsonfoulds ; goldstandard_dir_snprecision }
	//if (!goldstandard_dir.exists()) exit 1, "Input Gold standard path not found: ${params.goldstandard_dir}"
}

if(params.public_ref_dir){
	ref_dir = Channel.fromPath( params.public_ref_dir, type: 'dir' ,checkIfExists:true)
	//if (!ref_dir.exists()) exit 1, "Input Reference dir path not found: ${params.ref_dir}"
}

if(params.assess_dir){
	Channel
		.fromPath( params.assess_dir, type: 'dir' ,checkIfExists:true)
		.into { asses_dir_rfheatmap ; asses_dir_snprecision ; assess_dir_robinsonfoulds  }
	//if (!assess_dir.exists()) exit 1, "Input Asses dir path not found: ${params.assess_dir}"
}

params.tree_format = "newick"
if ( ! (params.tree_format =~ /newick|nexus/) ) {
	exit 1, 'Please provide a valid --tree_format option [newick,nexus]'
}

/*
* CHECK MANDATORY INPUTS
*/

params.input = false
if(! params.input){
	exit 1, "Missing tree input file : $params.input. Specify path with --input"
}

params.challenges_ids = false
if(! params.challenges_ids){
	exit 1, "Missing Event identifier : $params.challenges_ids. Specify path with --challenges_ids"
}

params.participant_id = false
if(! params.participant_id){
	exit 1, "Missing Participant identifier : $params.participant_id. Specify path with --participant_id"
}

result_dir = file(params.outdir)


/*
* HEADER LOG INFO
*/
log.info "========================================="
log.info " BU-ISCIII/openebench_gmi : OpenEBench pipeline for Outbreak detection challenge v${version}"
log.info "========================================="
def summary = [:]
summary['Test tree input']   = params.input
summary['Goldstandard dir']  = params.goldstandard_dir
summary['Public ref dir']    = params.public_ref_dir
summary['Benchmark data dir']  = params.assess_dir
summary['Event ID']            = params.challenges_ids
summary['Participant ID']      = params.participant_id
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']        = "$HOME"
summary['Current user']        = "$USER"
summary['Current path']        = "$PWD"
summary['Working dir']         = workflow.workDir
summary['Stats dir']           = params.statsdir
summary['Output dir']          = params.outdir
summary['Other Output dir']    = params.otherdir
summary['Script dir']          = workflow.projectDir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "===================================="


// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

/*
===============================
PIPELINE
===============================
*/


/*
* The instance generated from this docker file has to check the syntax of the submitted results.
*/
process validateInputFormat {

  container "openebench_gmi/sample-checkresults:" + version

  publishDir path: "${params.otherdir}", mode: 'copy', overwrite: true

  input:
  file tree from input_file

  output:
  file "*.nwk" into canonical_getresultsids,canonical_robinsonfoulds,canonical_snprecision

  """
  checkTreeFormat.py --tree_file ${tree} --tree_format ${params.tree_format} --output ${params.participant_id}_canonical.nwk
  """

}

/*
* The instance generated from this docker file knows how to extract query ids from the query.
*/
process getQueryIds {

  container "openebench_gmi/sample-getqueryids:" + version

  publishDir path: "${params.otherdir}", mode: 'copy', overwrite: true

  input:
  file tree from canonical_getresultsids


  output:
  file "*.json" into query_ids_json

  """
  getLeavesFromNewick.py --challenges_ids ${params.challenges_ids} --tree_file $tree --tree_format ${params.tree_format} --output queryids.json
  """

}

/*
* The instance generated from this docker file knows how to extract results ids from the results canonical formats.
*/
process ValidateInputIds {

  container "openebench_gmi/sample-compareids:" + version

  publishDir path: "${params.otherdir}", mode: 'copy', overwrite: true

  input:
  file query_ids from query_ids_json
  file ref_dir

  output:
  val task.exitStatus into EXIT_STAT_ROBINSONFOULDS,EXIT_STAT_SNPRECISION

  """
  compareIds.py --ids1 $query_ids --ids2 $ref_dir/inputIDs.json
  """

}

/*
* The instance generated from this docker file compute robinson foulds metric between the participant and all the public participants.
*/
process RobinsonFouldsMetrics {

  container "openebench_gmi/sample-robinsonfoulds:" + version

  publishDir path: "${params.otherdir}", mode: 'copy', overwrite: true

  input:
  val file_validated from EXIT_STAT_ROBINSONFOULDS
  file tree1 from canonical_robinsonfoulds
  file benchmark_dir from assess_dir_robinsonfoulds

  output:
  file "*.json" into metrics_robinsonfoulds_json

  when:
  file_validated == 0

  """
  mkdir aggregated_data
  cp -r $benchmark_dir/* aggregated_data/
  calculateRobinsonFouldsMetric.py --tree_file1 $tree1 --benchmark_trees_path aggregated_data -e ${params.challenges_ids} -p ${params.participant_id}
  """

}

/*
* The instance generated from this docker file compute metrics based on the number of lines and words.
*/
process SnPrecisionMetrics {

  container "openebench_gmi/sample-calculatesnprecision:" + version

  publishDir path: "${params.otherdir}", mode: 'copy', overwrite: true

  input:
  val file_validated from EXIT_STAT_SNPRECISION
  file tree1 from canonical_snprecision
  file gold_dir from goldstandard_dir_snprecision

  output:
  file "*.json" into metrics_snprecision_json

  when:
  file_validated == 0

  """
  calculateSnPrecision.py --tree_file1 $tree1 --tree_file2 $gold_dir/SIM-Sbareilly.tre -e ${params.challenges_ids} -p ${params.participant_id} -o ${params.participant_id}"_snprecision.json"
  """

}

process manage_assessment_snprecision {
	container "openebench_gmi/sample-assessment-snprecision:" + version
	tag "Performing benchmark assessment and building plots"

  	// publishDir path: "${params.outdir}", mode: 'copy', overwrite: true

	input:
	file assess_dir from asses_dir_snprecision
	file participant_result from metrics_snprecision_json
  val result_dir
  file new_data from file(params.otherdir)

	// output:
	// file "$result_dir/*" into benchmark_snprecision_result

	"""
	mkdir aggregated_data
	cp -r $assess_dir/* aggregated_data/
	cp $new_data/participant_matrix.json $new_data/*.nwk aggregated_data/
	python /app/manage_assessment_data.py -b aggregated_data -p $participant_result -o results
	cp -pr results/* $result_dir
	"""

}

process manage_assessment_rfheatmap {
	container "openebench_gmi/sample-assessment-rfheatmap:" + version
	tag "Performing benchmark assessment and building plots"

  	publishDir path: "${params.otherdir}", mode: 'copy', overwrite: true

	input:
	file assess_dir from asses_dir_rfheatmap
	file rb_metrics from metrics_robinsonfoulds_json
	output:
	file "heatmap/*" into benchmark_rfheatmap_result

	"""
	manageAssessmentRfHeatmap.py --assess_dir $assess_dir --output heatmap
	"""

}

