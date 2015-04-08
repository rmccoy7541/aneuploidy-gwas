%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: glm_gwas.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Rajiv McCoy
% This file performs the genome-wide association
% test, reading in the genotype and phenotype
% data and fitting a logistic GLM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = ('/sandbox/test.results'); % write results to this file
fid = fopen(file, 'a');

% then read in the phenotype matrix, which is Y x 2 dimensions where Y is the number
% of individuals with non-missing phenotype...the order must correspond to
% the SNP vector
pheno = csvread('/sandbox/world_ped_files/phenotype_files/mitotic_f.matlab.response', 0, 1);

for i = 1:23
% then read in SNPs one chromosome at a time, since reading in the entire
% file requires too much RAM
    gtid = strcat('/sandbox/world_ped_files/genotype_files/gwas/chr', num2str(i), '.matlab.csv');
    snps = csvread(gtid, 1, 6);
    snpid = strcat('/sandbox/world_ped_files/genotype_files/gwas/chr', num2str(i), '.rsid.csv');
    snpid = fopen(snpid);
    rsid = textscan(snpid, '%s', 'delimiter', ',');
    locid = strcat('/sandbox/world_ped_files/genotype_files/gwas/chr', num2str(i), '.map.csv');
    loc = csvread(locid, 0, 3);
    for j = 1:numel(rsid{:})    
    
% fit a binomial model with a logit link; include covariates if desired
% this one does not include any covariates
    [a,b,c] = glmfit(snps(:,j), [pheno(:,2) pheno(:,1) + pheno(:,2)], 'binomial', 'link', 'logit', 'estdisp', 'on');
    
% this test controls for age with a quadratic term
%    [a,b,c] = glmfit([snps(:,j) pheno(:,3) pheno(:,3).^2], [pheno(:,2) pheno(:,1) + pheno(:,2)], 'binomial', 'link', 'logit', 'estdisp', 'on');
% this test controls for age and sample type
%    [a,b,c] = glmfit([snps(:,j) pheno(:,3) pheno(:,4)], [pheno(:,2) pheno(:,1) + pheno(:,2)], 'binomial', 'link', 'logit', 'estdisp', 'on');
% this test controls for age with a quadratic term and ten principal components        
%    [a,b,c] = glmfit([snps(:,j) pheno(:,3) pheno(:,3).^2 pheno(:,4:13)], [pheno(:,2) pheno(:,1) + pheno(:,2)], 'binomial', 'link', 'logit', 'estdisp', 'on');
% this test controls for genotype at rs2305957 (i.e., second-order analysis) - requires that the printed
% output statistics also be changed to grab the third elements from results
% vectors
%    [a,b,c] = glmfit([pheno(:,3) snps(:,j)], [pheno(:,2) pheno(:,1) + pheno(:,2)], 'binomial', 'link', 'logit', 'estdisp', 'on');

% write the SNP ID, beta, SE, t, and p-value to disk.
        format_rs = rsid{:}{j};
        fprintf(fid, '%s\t', format_rs);
        fprintf(fid, '%d\t%d\t%f\t%f\t%f\t%e\n', [i loc(j) c.beta(2) c.se(2) c.t(2) c.p(2)]);
%        fprintf(fid, '%d\t%d\t%f\t%f\t%f\t%e\n', [i loc(j) c.beta(3) c.se(3) c.t(3) c.p(3)]);

    end 
    fclose(snpid);
end

fclose(fid);
