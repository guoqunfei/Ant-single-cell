#

# Transfrom gtf to gff
perl bin/gtf2gff.pl Amel.gtf.gz > Amel.gff
perl bin/check_orf_for_gff.pl Amel.gff Amel.fa > Amel.gff.orf

python bin/cleanPEP.py Amel.pep.fa.gz Amel.geneid_proteinid.cor > Amel.pep.clean

perl bin/bestProteinFromEnsembl.pl Amel.pep.clean Amel.gff.orf 1 > Amel.pep.best

perl bin/select_gff.pl Amel.pep.best Amel.gff > Amel.best.gff

perl bin/getGene.pl Amel.best.gff Amel.fa > Amel.best.cds

perl bin/check_orf_for_cds.pl Amel.best.cds | awk '$4>0' > Amel.best.cds.preStop

perl bin/cds2aa.pl Amel.best.cds > Amel.best.pep

perl bin/check_orf_for_gff.pl Amel.best.gff Amel.fa > Amel.best.gff.orf

awk '$8==1 && $9==1' Amel.best.gff.orf > Amel.best.gff.orf.intact

python bin/pep_specie.py Amel.best.pep > Amel.best.pep.fa

cat Mpha.best.pep.fa Amel.best.pep.fa > 2species.best.pep.fa

balst/bin/formatdb -i 2species.best.pep.fa -n 2speices.best.pep.fa -t 2species.best.pep.fa -l 2species.best.pep.fa.log -p T

mkdir split
for i in split; do cd $i; perl bin/Split_fasta.pl -seq 2species.best.pep.fa -nf 50 -od ./; cd ..;done
for i in split; do cd $i; perl bin/call_blast.pl ./ 2species.best.pep.fa blastp ./ best.pep.fa 0.01 queue_list project_id; cd ..;done
cat split/2species_*.best.pep.fa.m8 > 2species.best.pep.fa.m8

perl bin/select_m8.pl 2species.best.pep.fa.m8 config ./
for i in *.m8; do perl bin/solar.pl -a prot2prot -f m8 $i > $i.solar;done
for i in *.solar; do perl bin/solar_add_readLen.pl $i 2species.best.pep.fa > $i.cor;done
for i in *.m8; do perl bin/solar_add_identity.pl --solar $i.solar.cor --m8 $i > $i.solar.cor.idAdd;done

less Mpha_Amel.m8.solar.cor.idAdd.ort | grep -v "NA" | awk '{print $1"\t"$7}' | sed 's/Mpha_//g' | sed 's/Amel_//g' > Mpha_Amel.ort.gene2protein
less Amel.gtf.gz | grep "protein_id" | awk -F '\t' '{print $NF}' | while read p;od echo "${p}" | sed 's/;/\n/g' | grep -E 'gene_id|protein_id' | xargs echo;done | awk '{print $2"\t"$4}' | sort |uniq >> Amel.geneid_proteinid.cor
less Mpha_Amel.ort.gene2protein | while read p;do echo "${p}" | awk '{print $2}'|while read s;do less Amel.geneid_proteinid.cor.v1 | grep "^${s}\s" | while read t;do echo -ne "${p}\t${t}\n";done;done;done | awk '{print $1"\t"$4}' > Mpha_Amel.ort.genes
less Mpha_Amel.ort.genes | while read p;do echo "${p}" | awk '{print $1}' | while read s;do less Mpha.genes.anno | grep "^${s}\s" | while read t; do echo -ne "${p}\t${t}\n";done;done;done | awk '{print $1"\t"$4"\t"$2}' >> Mpha_Amel.ort.genes.anno