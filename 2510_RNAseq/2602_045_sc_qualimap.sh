


qualimap rnaseq \
-outdir results/qualimap/Mov10_oe_1 \
-a proportional \
-bam results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam \
-p strand-specific-reverse \
-gtf /n/groups/hbctraining/RNA_seq_part_1/reference_data/Homo_sapiens.GRCh38.92.1.gtf \
--java-mem-size=8G