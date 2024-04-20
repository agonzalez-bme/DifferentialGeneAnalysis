% Read the healthy and cancer RNA sequences into tables
clc; clear;
healthy_1 = readtable("BreastHealthy_1.tsv", "FileType","text",'Delimiter', '\t');
healthy_2 = readtable("BreastHealthy_2.tsv", "FileType","text",'Delimiter', '\t');
healthy_3 = readtable("BreastHealthy_3.tsv", "FileType","text",'Delimiter', '\t');
healthy_4 = readtable("BreastHealthy_4.tsv", "FileType","text",'Delimiter', '\t');
healthy_5 = readtable("BreastHealthy_5.tsv", "FileType","text",'Delimiter', '\t');
healthy_6 = readtable("BreastHealthy_6.tsv", "FileType","text",'Delimiter', '\t');

cancer_1 = readtable("BreastCancer_1.tsv", "FileType","text",'Delimiter', '\t');
cancer_2 = readtable("BreastCancer_2.tsv", "FileType","text",'Delimiter', '\t');
cancer_3 = readtable("BreastCancer_3.tsv", "FileType","text",'Delimiter', '\t');
cancer_4 = readtable("BreastCancer_4.tsv", "FileType","text",'Delimiter', '\t');
cancer_5 = readtable("BreastCancer_5.tsv", "FileType","text",'Delimiter', '\t');
cancer_6 = readtable("BreastCancer_6.tsv", "FileType","text",'Delimiter', '\t');

% Extract the name of genes
gene_names = healthy_1.gene_name;
gene_names(1:4,:) = [];

% We are putting all the normalized gene expression next to eachother basically
% One sample per column, one gene per row
tpm_c = [cancer_1.tpm_unstranded,cancer_2.tpm_unstranded,cancer_3.tpm_unstranded,cancer_4.tpm_unstranded,cancer_5.tpm_unstranded,cancer_6.tpm_unstranded];
tpm_c(1:4,:) = [];
tpm_c = log2(tpm_c);

tpm_h = [healthy_1.tpm_unstranded,healthy_2.tpm_unstranded,healthy_3.tpm_unstranded,healthy_4.tpm_unstranded,healthy_5.tpm_unstranded,healthy_6.tpm_unstranded];
tpm_h(1:4,:) = [];
tpm_h = log2(tpm_h);

h = zeros(60660,1);
p = zeros(60660,1);
for ii = 1:length(tpm_h)
    [h(ii), p(ii)] = ttest2(tpm_h(ii,:), tpm_c(ii,:));
end

mavolcanoplot(tpm_c,tpm_h,p,'Labels',gene_names)

