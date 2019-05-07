#!/usr/bin/env python3

import sys
import os
import argparse
import json
import shutil
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("SVG")
plt.ioff()

#################
### FUNCTIONS ###
#################

def check_arg (args=None) :
    '''
    Description:
        Function collect arguments from command line using argparse
    Input:
        args # command line arguments
    Constant:
        None
    Variables
        parser
    Return
        parser.parse_args() # Parsed arguments

    '''
    parser = argparse.ArgumentParser(prog = 'calculateSnPrecision.py', formatter_class=argparse.RawDescriptionHelpFormatter, description= 'calculateSpPrecision.py calculates Specificity and Precision metric between two newick trees using RAxML.')

    parser.add_argument('--version','-v', action='version', version='%(prog)s 0.3.5')
    parser.add_argument('--assess_dir','-a' ,required= False,choices = ["newick","nexus"], help = 'Tree file format [newick,nexus]', default = "newick")
    parser.add_argument('--output' ,'-o',required= False, help = 'Path to result folder. Default = benchmark_result', default="benchmark_result")
    parser.add_argument('--event_id','-e' ,required= False, help = 'OpenEbench event identifier', default="default")
    parser.add_argument('--participant_id','-p' ,required= False, help = 'OpenEbench participant identifier', default="default")

    return parser.parse_args()

###################
### MAIN SCRIPT ###
###################

if __name__ == '__main__' :

    # Variables
    version = 'manageAssesmentRbHeatmap v1.0'  # Script version
    arguments = ""                        # Arguments from ArgParse

    # Grab arguments
    arguments = check_arg(sys.argv[1:])

    # Checking if file exists
    if (not os.path.exists(arguments.assess_dir)):
        print("Assess directory not found.")
        sys.exit(1)
    in_file = os.path.join(arguments.assess_dir,"participant_matrix.json")
    with open(in_file) as json_file:
        data = json.load(json_file)

    #Obtain participants and data matrix
    participants = data['participants']
    p_matrix_values = np.array(data['matrix']['values'])

    # Draw heatmap
    fig, ax = plt.subplots()
    im = ax.imshow(p_matrix_values)


    # We want to show all ticks...
    ax.set_xticks(np.arange(len(participants)))
    ax.set_yticks(np.arange(len(participants)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(participants)
    ax.set_yticklabels(participants)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    #for i in range(len(p_matrix_values)):
    #    for j in range(len(p_matrix_values)):
    #        text = ax.text(j, i, harvest[i, j],
    #                       ha="center", va="center", color="w")

    ax.set_title("Heat map for participants)")
    fig.tight_layout()
    outname = os.path.join(arguments.output,"benchmark_gmi_robinsonfoulds_heatmap.svg")
    fig.set_size_inches(18.5, 10.5)
    fig.savefig(outname, dpi=100)
