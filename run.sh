#! /usr/bin/bash
dir=$1
name=$2
outfolder=$3
fusionRDS=$4
nsim=5
ngene=2
## do ericscript
./ericscript.pl --simulator --db ../ericscript_db_homosapiens_ensembl73/ -o $dir/$name --ngenefusion $ngene --nsims $nsim -rl 75 -v --background_1 /home/ec2-user/data/bio-5943-mojo-synthetic/background_data/20-A73997_RSQ1/20-A73997_RSQ1_1.fastq --background_2 /home/ec2-user/data/bio-5943-mojo-synthetic/background_data/20-A73997_RSQ1/20-A73997_RSQ1_3.fastq --nreads_background 200000 --min_cov 1000 --max_cov 1000 --fusionlistRDS $fusionRDS
folder="$dir"/"$name"/IE/reads/
## compress fastq
for ((sim=1;sim<=$nsim;sim++));
do
        sim1=$(printf "%05d" $sim)
        mkdir -p $outfolder/"$name"_"$sim1"/
        cp $folder/sim_$sim1/total.reads.1.fq $outfolder/"$name"_"$sim1"/"$name"_"$sim1"_1.fastq
        gzip $outfolder/"$name"_"$sim1"/"$name"_"$sim1"_1.fastq
        cp $folder/sim_$sim1/total.reads.2.fq $outfolder/"$name"_"$sim1"/"$name"_"$sim1"_3.fastq
        gzip $outfolder/"$name"_"$sim1"/"$name"_"$sim1"_3.fastq
	## add truth sheet to fastq folder
	Rscript ./print_truth_fusion.R $outfolder/"$name"_"$sim1"/
	tar -cvzf $outfolder/"$name"_"$sim1".tar.gz $outfolder/"$name"_"$sim1"
done
