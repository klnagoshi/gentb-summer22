#!/bin/bash
#SBATCH -c 4                               # Request one core
#SBATCH -t 0-11:59                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=100G                         # Memory total in MiB (for all cores)
#SBATCH -o generate-mutation-db_7-29_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e generate-mutation-db_7-29_%j.err                 # File to which STDERR will be written, including job ID (%j)
                                           # You can change the filenames given with -o and -e to any filenames you'd like
#SBATCH --mail-type=ALL                    # ALL email notification type
#SBATCH --mail-user=klnagoshi@college.harvard.edu  # Email to which notifications will be sent


# You can change hostname to any command you would like to run
python3 generate-mutation-db.py