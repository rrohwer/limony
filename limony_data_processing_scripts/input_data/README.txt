taxass_limonyIRD_silva138cyanos

RRR 2021/08/25

repeating TaxAss on limony with IRD QC pipeline
with the silva 138 database where I edited cyano names

Starting with the dada2 files created on 2019-07-24
	Ie.03_R1e2t247_R2e2t195_dmx_dada_seqtab_nochim.rds
	Ie.03_R1e2t247_R2e2t195_dmx_dada_nochim_tracking.rds

Re-Running TaxAss with cyano-edited silva138 files:
	FreshTrain15Jun2020silva138.fasta
	FreshTrain15Jun2020silva138.taxonomy
	silva_nr_v138_taxass.fasta
	silva_nr_v138_taxass_cyano_edits.taxonomy

Using tax-scripts that are TaxAss V2, latest commit:
	commit 1268707da8388328bdee2f35d93880d8b50b4af2 (HEAD -> master, origin/master)
	Author: Robin Rohwer <Robin.Rohwer@gmail.com>
	Date:   Wed Aug 25 12:41:44 2021 -0700
    fixes a bug when finding forcing diffs (database improvement plot). Fixes converting everything to unclassified (=0) for comparison purposes.
    seems like it was working just fine, but it gave a bunch of warning messages.

------------------------------------------------------------------------------------------

All Starting Files:
$ ls -lh
total 490312
-rw-r--r--@ 1 athena  staff    30K Aug 19 18:50 Compare_Edited_Cyano_and_Silva138_Names.xlsx
-rw-r--r--@ 1 athena  staff    19K Jul  7  2020 Compare_FreshTrain_and_Silva138_Names.xlsx
-rw-r--r--  1 athena  staff   1.8M Jul  1  2020 FreshTrain15Jun2020silva138.fasta
-rw-r--r--  1 athena  staff   106K Jul  1  2020 FreshTrain15Jun2020silva138.taxonomy
-rw-r--r--  1 athena  staff    22K Jul 24  2019 Ie.03_R1e2t247_R2e2t195_dmx_dada_nochim_tracking.rds
-rw-r--r--  1 athena  staff   3.9M Jul 24  2019 Ie.03_R1e2t247_R2e2t195_dmx_dada_seqtab_nochim.rds
-rw-r--r--@ 1 athena  staff   635K Jul  1  2020 README-138.html
-rw-r--r--@ 1 athena  staff   731K Aug 19 19:33 README-138_cyano_edits.html
-rw-r--r--  1 athena  staff   433B Sep  4  2018 README.md
-rw-r--r--@ 1 athena  staff   787B Aug 25 10:48 README.txt
-rwxr-xr-x@ 1 athena  staff   2.9K Jul  4  2020 RunStep_15.sh
-rwxr-xr-x@ 1 athena  staff   1.8K Jul  4  2020 RunStep_16.sh
-rwxr-xr-x@ 1 athena  staff   7.3K Jul  2  2020 RunSteps_1-14.sh
-rwxr-xr-x  1 athena  staff   4.3K Jul  2  2020 RunSteps_quickie.sh
-rw-r--r--@ 1 athena  staff   946K Jul 14  2020 TaxAss_Directions.html
-rwxrwxrwx@ 1 athena  staff   6.0K May  1  2017 calc_full_length_pident.R
-rwxrwxrwx  1 athena  staff   1.7K Feb  2  2016 create_fastas_given_seqIDs.py
-rwxrwxrwx  1 athena  staff   394B Oct  9  2017 deletemothurbarf.sh
-rwxrwxrwx  1 athena  staff   3.4K Aug 23  2016 filter_seqIDs_by_pident.R
-rwxr-xr-x@ 1 athena  staff    33K Aug 25 12:41 find_classification_disagreements.R
-rwxrwxrwx  1 athena  staff   2.0K Feb  2  2016 find_seqIDs_blast_removed.py
-rwxr-xr-x  1 athena  staff    13K Jul  2  2020 plot_blast_hit_stats.R
-rwxrwxrwx  1 athena  staff    60K Jul 14  2020 plot_classification_disagreements.R
-rwxrwxrwx  1 athena  staff    21K Jun  4  2018 plot_classification_improvement.R
-rwxrwxrwx@ 1 athena  staff   4.2K Jul 24  2019 reformat_dada2_seqtabs.R
-rw-r--r--  1 athena  staff   1.6K Jul  2  2020 reformat_fasta.R
-rwxrwxrwx@ 1 athena  staff   3.8K May  1  2017 reformat_mothur_OTU_tables.R
-rwxr-xr-x@ 1 athena  staff    15K Nov  9  2020 reformat_taxonomy_nomenclature.R
-rw-r--r--  1 athena  staff   210M Jul  1  2020 silva_nr_v138_taxass.fasta
-rw-r--r--@ 1 athena  staff    18M Aug 19 17:59 silva_nr_v138_taxass_cyano_edits.taxonomy

List of commands run:
$ mv silva_nr_v138_taxass.fasta silva_nr_v138_taxass_cyano_edits.fasta  # "general" names need to match
$ Rscript reformat_dada2_seqtabs.R Ie.03_R1e2t247_R2e2t195_dmx_dada_seqtab_nochim.rds IRD.fasta IRD.abund IRD.count
$ ./RunSteps_1-14.sh IRD FreshTrain15Jun2020silva138 silva_nr_v138_taxass_cyano_edits "100 99 98 97 96 95" 80 80 2 > termout.txt
$ ./RunStep_15.sh IRD FreshTrain15Jun2020silva138 silva_nr_v138_taxass_cyano_edits 98 80 80 2 >> termout.txt
$ ./RunStep_16.sh IRD FreshTrain15Jun2020silva138 silva_nr_v138_taxass_cyano_edits


