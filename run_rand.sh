#! /usr/bin/bash
./ericscript_tempus_rand.pl --simulator --db ../ericscript_db_homosapiens_ensembl73/ -o $1  --nsims 5 -rl 75 -v --background_1 /home/ec2-user/data/bio-5943-mojo-synthetic/background_data/20-A73997_RSQ1/20-A73997_RSQ1_1.fastq --background_2 /home/ec2-user/data/bio-5943-mojo-synthetic/background_data/20-A73997_RSQ1/20-A73997_RSQ1_3.fastq --nreads_background 200000 --min_cov 80 --max_cov 120 -ngenefusion 35
