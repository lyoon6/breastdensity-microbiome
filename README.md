# The association between breast density and gut microbiota composition at 2 years post menarche: A cross-sectional study of adolescents in Santiago, Chile

### Author: Lara Yoon (lyoon6@ucla.edu) 


### Use these R scripts to perform diversity and maaslin analyses 




Note: The QIIME2 PICRUSt2 plugin was used to produce metagenome predictions. The code below can be implemented in QIIME2 and was adapted from the q2-picrust2 tutorial: https://github.com/picrust/picrust2/wiki/q2-picrust2-Tutorial 

 qiime picrust2 full-pipeline --i-table /home/qiime2/q2-filtering/meta-taxa-samp-filt-table.qza --i-seq /media/sf_qiime2sf/representative_sequences.qza --output-dir q2-picrust2_aim3 --p-placement-tool sepp --p-threads 1 --p-hsp-method pic --p-max-nsti 2 –verbose
