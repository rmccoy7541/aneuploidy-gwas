%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: calculate_pca_scores.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Rajiv McCoy
% This file performs PCA on HapMap data then 
% loops through client genotype data to compute 
% principal component scores on those predefined 
% axes. This is for the purpose of defining sub-
% samples for stratified analysis. The program
% KING was used to calculate PCA scores within
% the sample for use as GWAS covariates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%client file directory
files = dir('/sandbox/HapMap3/client_gt_input/*.hapmap');

disp('Loading data matrix...');
A = csvread('hapmap3.ref_recode.sorted.pca', 0, 3);

% since PCA does not allow for missing data, need to impute missing data; in this case, we use k nearest neighbor imputation
disp('Data matrix loaded. Imputing missing data...');
B = knnimpute(A);
clear A;

% select 20K random SNPs for PCA    
disp('Imputation complete. Randomly downsampling matrix...');
sample=randsample(numel(B(:,1)), 20000);
SUBSET = B(sample, :);
clear B;
disp('Downsampling complete. Computing principle components...');
[COEFF, SCORE, VARIANCE] = princomp(SUBSET');


% define PC ranges for different population groupings -- only EUR and ASN groups has sufficient sample sizes to present in paper
RANGE_pc1_EUR = [min([SCORE(88:252,1) ; SCORE(1093:1194,1)]) max([SCORE(88:252,1) ; SCORE(1093:1194,1)])];
RANGE_pc2_EUR = [min([SCORE(88:252,2) ; SCORE(1093:1194,2)]) max([SCORE(88:252,2) ; SCORE(1093:1194,2)])];
RANGE_pc3_EUR = [min([SCORE(88:252,3) ; SCORE(1093:1194,3)]) max([SCORE(88:252,3) ; SCORE(1093:1194,3)])];

RANGE_pc1_ASIA = [min([SCORE(253:498,1) ; SCORE(600:712,1)]) max([SCORE(253:498,1) ; SCORE(600:712,1)])];
RANGE_pc2_ASIA = [min([SCORE(253:498,2) ; SCORE(600:712,2)]) max([SCORE(253:498,2) ; SCORE(600:712,2)])];
RANGE_pc3_ASIA = [min([SCORE(253:498,3) ; SCORE(600:712,3)]) max([SCORE(253:498,3) ; SCORE(600:712,3)])];

RANGE_pc1_AFR = [min([SCORE(1:87,1) ; SCORE(713:822,1) ; SCORE(909:1092,1) ; SCORE(1195:1397,1)]) max([SCORE(1:87,1) ; SCORE(713:822,1) ; SCORE(909:1092,1) ; SCORE(1195:1397,1)])];
RANGE_pc2_AFR = [min([SCORE(1:87,2) ; SCORE(713:822,2) ; SCORE(909:1092,2) ; SCORE(1195:1397,2)]) max([SCORE(1:87,2) ; SCORE(713:822,2) ; SCORE(909:1092,2) ; SCORE(1195:1397,2)])];
RANGE_pc3_AFR = [min([SCORE(1:87,3) ; SCORE(713:822,3) ; SCORE(909:1092,3) ; SCORE(1195:1397,3)]) max([SCORE(1:87,3) ; SCORE(713:822,3) ; SCORE(909:1092,3) ; SCORE(1195:1397,3)])];

RANGE_pc1_ASW = [min(SCORE(1:87,1)) max(SCORE(1:87,1))]; 
RANGE_pc2_ASW = [min(SCORE(1:87,2)) max(SCORE(1:87,2))]; 
RANGE_pc3_ASW = [min(SCORE(1:87,3)) max(SCORE(1:87,3))]; 

RANGE_pc1_CEU = [min(SCORE(88:252,1)) max(SCORE(88:252,1))]; 
RANGE_pc2_CEU = [min(SCORE(88:252,2)) max(SCORE(88:252,2))]; 
RANGE_pc3_CEU = [min(SCORE(88:252,3)) max(SCORE(88:252,3))]; 

RANGE_pc1_CHB = [min(SCORE(253:389,1)) max(SCORE(253:389,1))]; 
RANGE_pc2_CHB = [min(SCORE(253:389,2)) max(SCORE(253:389,2))]; 
RANGE_pc3_CHB = [min(SCORE(253:389,3)) max(SCORE(253:389,3))]; 

RANGE_pc1_CHD = [min(SCORE(390:498,1)) max(SCORE(390:498,1))]; 
RANGE_pc2_CHD = [min(SCORE(390:498,2)) max(SCORE(390:498,2))]; 
RANGE_pc3_CHD = [min(SCORE(390:498,3)) max(SCORE(390:498,3))]; 

RANGE_pc1_GIH = [min(SCORE(499:599,1)) max(SCORE(499:599,1))]; 
RANGE_pc2_GIH = [min(SCORE(499:599,2)) max(SCORE(499:599,2))]; 
RANGE_pc3_GIH = [min(SCORE(499:599,3)) max(SCORE(499:599,3))]; 

RANGE_pc1_JPT = [min(SCORE(600:712,1)) max(SCORE(600:712,1))]; 
RANGE_pc2_JPT = [min(SCORE(600:712,2)) max(SCORE(600:712,2))]; 
RANGE_pc3_JPT = [min(SCORE(600:712,3)) max(SCORE(600:712,3))]; 

RANGE_pc1_LWK = [min(SCORE(713:822,1)) max(SCORE(713:822,1))]; 
RANGE_pc2_LWK = [min(SCORE(713:822,2)) max(SCORE(713:822,2))]; 
RANGE_pc3_LWK = [min(SCORE(713:822,3)) max(SCORE(713:822,3))]; 

RANGE_pc1_MEX = [min(SCORE(823:908,1)) max(SCORE(823:908,1))]; 
RANGE_pc2_MEX = [min(SCORE(823:908,2)) max(SCORE(823:908,2))]; 
RANGE_pc3_MEX = [min(SCORE(823:908,3)) max(SCORE(823:908,3))]; 

RANGE_pc1_MKK = [min(SCORE(909:1092,1)) max(SCORE(909:1092,1))]; 
RANGE_pc2_MKK = [min(SCORE(909:1092,2)) max(SCORE(909:1092,2))]; 
RANGE_pc3_MKK = [min(SCORE(909:1092,3)) max(SCORE(909:1092,3))]; 

RANGE_pc1_TSI = [min(SCORE(1093:1194,1)) max(SCORE(1093:1194,1))]; 
RANGE_pc2_TSI = [min(SCORE(1093:1194,2)) max(SCORE(1093:1194,2))]; 
RANGE_pc3_TSI = [min(SCORE(1093:1194,3)) max(SCORE(1093:1194,3))]; 

RANGE_pc1_YRI = [min(SCORE(1195:1397,1)) max(SCORE(1195:1397,1))]; 
RANGE_pc2_YRI = [min(SCORE(1195:1397,2)) max(SCORE(1195:1397,2))]; 
RANGE_pc3_YRI = [min(SCORE(1195:1397,3)) max(SCORE(1195:1397,3))]; 


% now calculate client PC scores on those pre-defined axes
disp('PCA complete. Looping through client data...');

CLIENT_MATRIX = [];
for m = 1:numel(files)   
    CLIENT = load(files(m).name);
    CLIENT_MATRIX = [CLIENT_MATRIX CLIENT(sample, :)];
end


disp('Client data loaded. Imputing missing data...');
CLIENT_MATRIX = [SUBSET CLIENT_MATRIX];
CLIENT_MATRIX = knnimpute(CLIENT_MATRIX);
CLIENT_IMPUTE = CLIENT_MATRIX(:,1398:(1398 + (numel(files) * 2) - 1));
clear CLIENT_MATRIX;

pc1_list = [];
pc2_list = [];
pc3_list = [];
fid = fopen('/sandbox/HapMap3/inferred_ancestry.20000.1.txt', 'a'); % choose output location
l = 1;
one_two = 1;
one_two_num = 1;
for m = 1:numel(CLIENT_IMPUTE(1,:))
     if mod(m,2) == 0,
         %father is indicated by a "one" in output file
         one_two=' father';
         one_two_num=1;
     else
         %mother is indicated by a "two", consistent with PLINK
         one_two=' mother';
         one_two_num=2;
     end    
     family = files(l).name;
     family = regexp(family, '\.', 'split');
     family = family{:};
     
     
     % matlab centers the data (subtracts the mean) before computing scores, so
     % need to center the new observation vector before converting to PCs
     MEAN = mean(SUBSET,2);
     % (SUBSET(:,1)-MEAN(:,1))' * COEFF(:,1)
     % SCORE(1,1)
 
     CLIENT_SUBSET = CLIENT_IMPUTE(:, m);
     pc1 = (CLIENT_SUBSET(:,1)-MEAN(:,1))'*COEFF(:,1);
     pc2 = (CLIENT_SUBSET(:,1)-MEAN(:,1))'*COEFF(:,2);
     pc3 = (CLIENT_SUBSET(:,1)-MEAN(:,1))'*COEFF(:,3);
 
     pc1_list = [pc1_list, pc1];
     pc2_list = [pc2_list, pc2];
     pc3_list = [pc3_list, pc3];
 
     if (pc1 > RANGE_pc1_EUR(1)) && (pc1 < RANGE_pc1_EUR(2)) && (pc2 > RANGE_pc2_EUR(1)) && (pc2 < RANGE_pc2_EUR(2)) && (pc3 > RANGE_pc3_EUR(1)) && (pc3 < RANGE_pc3_EUR(2))
         disp(['Family ', num2str(l), num2str(one_two), ' within bounds of European principle components.']);
         EUR_ANC = 1;
     else
         disp(['Family ', num2str(l), num2str(one_two), ' lies outside the bounds of European principle components.']);
         EUR_ANC = 0;
     end
     
     
     if (pc1 > RANGE_pc1_ASIA(1)) && (pc1 < RANGE_pc1_ASIA(2)) && (pc2 > RANGE_pc2_ASIA(1)) && (pc2 < RANGE_pc2_ASIA(2)) && (pc3 > RANGE_pc3_ASIA(1)) && (pc3 < RANGE_pc3_ASIA(2))
         disp(['Family ', num2str(l), num2str(one_two), ' within bounds of Asian principle components.']);
         ASIA_ANC = 1;
     else
         disp(['Family ', num2str(l), num2str(one_two), ' lies outside the bounds of Asian principle components.']);
         ASIA_ANC = 0;
     end
     
     
     if (pc1 > RANGE_pc1_AFR(1)) && (pc1 < RANGE_pc1_AFR(2)) && (pc2 > RANGE_pc2_AFR(1)) && (pc2 < RANGE_pc2_AFR(2)) && (pc3 > RANGE_pc3_AFR(1)) && (pc3 < RANGE_pc3_AFR(2))
         disp(['Family ', num2str(l), num2str(one_two), ' within bounds of African principle components.']);
         AFR_ANC = 1;
     else
         disp(['Family ', num2str(l), num2str(one_two), ' lies outside the bounds of African principle components.']);
         AFR_ANC = 0;
     end     
     
     
     if (pc1 > RANGE_pc1_ASW(1)) && (pc1 < RANGE_pc1_ASW(2)) && (pc2 > RANGE_pc2_ASW(1)) && (pc2 < RANGE_pc2_ASW(2)) && (pc3 > RANGE_pc3_ASW(1)) && (pc3 < RANGE_pc3_ASW(2))
         disp(['Family ', num2str(l), num2str(one_two), ' within bounds of ASW principle components.']);
         ASW_ANC = 1;
     else
         disp(['Family ', num2str(l), num2str(one_two), ' lies outside the bounds of ASW principle components.']);
         ASW_ANC = 0;
     end     

     if (pc1 > RANGE_pc1_CEU(1)) && (pc1 < RANGE_pc1_CEU(2)) && (pc2 > RANGE_pc2_CEU(1)) && (pc2 < RANGE_pc2_CEU(2)) && (pc3 > RANGE_pc3_CEU(1)) && (pc3 < RANGE_pc3_CEU(2))
         disp(['Family ', num2str(l), num2str(one_two), ' within bounds of CEU principle components.']);
         CEU_ANC = 1;
     else
         disp(['Family ', num2str(l), num2str(one_two), ' lies outside the bounds of CEU principle components.']);
         CEU_ANC = 0;
     end     
     
     if (pc1 > RANGE_pc1_CHB(1)) && (pc1 < RANGE_pc1_CHB(2)) && (pc2 > RANGE_pc2_CHB(1)) && (pc2 < RANGE_pc2_CHB(2)) && (pc3 > RANGE_pc3_CHB(1)) && (pc3 < RANGE_pc3_CHB(2))
         disp(['Family ', num2str(l), num2str(one_two), ' within bounds of CHB principle components.']);
         CHB_ANC = 1;
     else
         disp(['Family ', num2str(l), num2str(one_two), ' lies outside the bounds of CHB principle components.']);
         CHB_ANC = 0;
     end     

     if (pc1 > RANGE_pc1_CHD(1)) && (pc1 < RANGE_pc1_CHD(2)) && (pc2 > RANGE_pc2_CHD(1)) && (pc2 < RANGE_pc2_CHD(2)) && (pc3 > RANGE_pc3_CHD(1)) && (pc3 < RANGE_pc3_CHD(2))
         disp(['Family ', num2str(l), num2str(one_two), ' within bounds of CHD principle components.']);
         CHD_ANC = 1;
     else
         disp(['Family ', num2str(l), num2str(one_two), ' lies outside the bounds of CHD principle components.']);
         CHD_ANC = 0;
     end     
     
     if (pc1 > RANGE_pc1_GIH(1)) && (pc1 < RANGE_pc1_GIH(2)) && (pc2 > RANGE_pc2_GIH(1)) && (pc2 < RANGE_pc2_GIH(2)) && (pc3 > RANGE_pc3_GIH(1)) && (pc3 < RANGE_pc3_GIH(2))
         disp(['Family ', num2str(l), num2str(one_two), ' within bounds of GIH principle components.']);
         GIH_ANC = 1;
     else
         disp(['Family ', num2str(l), num2str(one_two), ' lies outside the bounds of GIH principle components.']);
         GIH_ANC = 0;
     end     
     
     if (pc1 > RANGE_pc1_JPT(1)) && (pc1 < RANGE_pc1_JPT(2)) && (pc2 > RANGE_pc2_JPT(1)) && (pc2 < RANGE_pc2_JPT(2)) && (pc3 > RANGE_pc3_JPT(1)) && (pc3 < RANGE_pc3_JPT(2))
         disp(['Family ', num2str(l), num2str(one_two), ' within bounds of JPT principle components.']);
         JPT_ANC = 1;
     else
         disp(['Family ', num2str(l), num2str(one_two), ' lies outside the bounds of JPT principle components.']);
         JPT_ANC = 0;
     end     
     
     if (pc1 > RANGE_pc1_LWK(1)) && (pc1 < RANGE_pc1_LWK(2)) && (pc2 > RANGE_pc2_LWK(1)) && (pc2 < RANGE_pc2_LWK(2)) && (pc3 > RANGE_pc3_LWK(1)) && (pc3 < RANGE_pc3_LWK(2))
         disp(['Family ', num2str(l), num2str(one_two), ' within bounds of LWK principle components.']);
         LWK_ANC = 1;
     else
         disp(['Family ', num2str(l), num2str(one_two), ' lies outside the bounds of LWK principle components.']);
         LWK_ANC = 0;
     end     
     
     if (pc1 > RANGE_pc1_MEX(1)) && (pc1 < RANGE_pc1_MEX(2)) && (pc2 > RANGE_pc2_MEX(1)) && (pc2 < RANGE_pc2_MEX(2)) && (pc3 > RANGE_pc3_MEX(1)) && (pc3 < RANGE_pc3_MEX(2))
         disp(['Family ', num2str(l), num2str(one_two), ' within bounds of MEX principle components.']);
         MEX_ANC = 1;
     else
         disp(['Family ', num2str(l), num2str(one_two), ' lies outside the bounds of MEX principle components.']);
         MEX_ANC = 0;
     end     
     
     if (pc1 > RANGE_pc1_MKK(1)) && (pc1 < RANGE_pc1_MKK(2)) && (pc2 > RANGE_pc2_MKK(1)) && (pc2 < RANGE_pc2_MKK(2)) && (pc3 > RANGE_pc3_MKK(1)) && (pc3 < RANGE_pc3_MKK(2))
         disp(['Family ', num2str(l), num2str(one_two), ' within bounds of MKK principle components.']);
         MKK_ANC = 1;
     else
         disp(['Family ', num2str(l), num2str(one_two), ' lies outside the bounds of MKK principle components.']);
         MKK_ANC = 0;
     end     
     
     if (pc1 > RANGE_pc1_TSI(1)) && (pc1 < RANGE_pc1_TSI(2)) && (pc2 > RANGE_pc2_TSI(1)) && (pc2 < RANGE_pc2_TSI(2)) && (pc3 > RANGE_pc3_TSI(1)) && (pc3 < RANGE_pc3_TSI(2))
         disp(['Family ', num2str(l), num2str(one_two), ' within bounds of TSI principle components.']);
         TSI_ANC = 1;
     else
         disp(['Family ', num2str(l), num2str(one_two), ' lies outside the bounds of TSI principle components.']);
         TSI_ANC = 0;
     end     
     
     if (pc1 > RANGE_pc1_YRI(1)) && (pc1 < RANGE_pc1_YRI(2)) && (pc2 > RANGE_pc2_YRI(1)) && (pc2 < RANGE_pc2_YRI(2)) && (pc3 > RANGE_pc3_YRI(1)) && (pc3 < RANGE_pc3_YRI(2))
         disp(['Family ', num2str(l), num2str(one_two), ' within bounds of YRI principle components.']);
         YRI_ANC = 1;
     else
         disp(['Family ', num2str(l), num2str(one_two), ' lies outside the bounds of YRI principle components.']);
         YRI_ANC = 0;
     end          
     
     
     %write family ID, parent sex, ancestry (1=European, 0=other), and first 3 PC scores to file
     RESULTS = sprintf('%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f', [str2num(family) one_two_num EUR_ANC ASIA_ANC AFR_ANC ASW_ANC CEU_ANC CHB_ANC CHD_ANC GIH_ANC JPT_ANC LWK_ANC MEX_ANC MKK_ANC TSI_ANC YRI_ANC pc1 pc2 pc3]);
     fprintf(fid, '%s\n', RESULTS);  
     if mod(m,2)==0,
         l=l+1;
     end
end
clear SUBSET; clear CLIENT_SUBSET; clear CLIENT_IMPUTE;

%following code is for plotting

scatter(SCORE(:,1), SCORE(:,2), 2, 'black'); hold on; scatter(pc1_list(:), pc2_list(:), 2, 'green');

% scatter(SCORE(1:87,1), SCORE(1:87,2), 2, 'black'); %%ASW Africans in SW USA 6
% hold on;
% scatter(SCORE(88:252,1), SCORE(88:252,2), 2, 'blue'); %%CEU N/W Europeans in Utah 7
% hold on;
% scatter(SCORE(253:389,1), SCORE(253:389,2), 2, 'black'); %%CHB Han Chinese Beijing 8
% hold on;
% scatter(SCORE(390:498,1), SCORE(390:498,2), 2, 'black'); %%CHD Chinese in Denver 9
% hold on;
% scatter(SCORE(499:599,1), SCORE(499:599,2), 2, 'black'); %GIH Gujarati in Houston 10
% hold on;
% scatter(SCORE(600:712,1), SCORE(600:712,2), 2, 'black'); %%JPT Japanese in Tokyo 11
% hold on;
% scatter(SCORE(713:822,1), SCORE(713:822,2), 2, 'black'); %%LWK Luhya (Bantu) in Kenya 12
% hold on;
% scatter(SCORE(823:908,1), SCORE(823:908,2), 2, 'black'); %MEX Mexican in LA 13
% hold on;
% scatter(SCORE(909:1092,1), SCORE(909:1092,2), 2, 'black'); %%MKK Masaii in Kenya 14
% hold on;
% scatter(SCORE(1093:1194,1), SCORE(1093:1194,2), 2, 'black'); %%TSI Toscans in Italy 15
% hold on;
% scatter(SCORE(1195:1397,1), SCORE(1195:1397,2), 2, 'black'); %%YRI Yoruba in Nigeria 16
% hold on;
% scatter(pc1_list(:), pc2_list(:), 2, 'green' )
% hold on;
% rectangle('Position',[RANGE_pc1_AFR(1),RANGE_pc2_AFR(1),(RANGE_pc1_AFR(2)-RANGE_pc1_AFR(1)),(RANGE_pc2_AFR(2)-RANGE_pc2_AFR(1))]);

% figure;
% plot(SCORE(1:87,2), SCORE(1:87,3),'.', 'Color', 'black'); %ASW Africans in SW USA
% hold on;
% plot(SCORE(88:252,2), SCORE(88:252,3),'.', 'Color', 'blue'); %CEU N/W Europeans in Utah
% hold on;
% plot(SCORE(253:389,2), SCORE(253:389,3),'.', 'Color', 'black'); %CHB Han Chinese Beijing
% hold on;
% plot(SCORE(390:498,2), SCORE(390:498,3),'.', 'Color', 'black'); %CHD Chinese in Denver
% hold on;
% plot(SCORE(499:599,2), SCORE(499:599,3),'.', 'Color', 'black'); %GIH Gujarati in Houston
% hold on;
% plot(SCORE(600:712,2), SCORE(600:712,3),'.', 'Color', 'black'); %JPT Japanese in Tokyo
% hold on;
% plot(SCORE(713:822,2), SCORE(713:822,3),'.', 'Color', 'black'); %LWK Luhya (Bantu) in Kenya
% hold on;
% plot(SCORE(823:908,2), SCORE(823:908,3),'.', 'Color', 'black'); %MEX Mexican in LA
% hold on;
% plot(SCORE(909:1092,2), SCORE(909:1092,3),'.', 'Color', 'black'); %MKK Masaii in Kenya
% hold on;
% plot(SCORE(1093:1194,2), SCORE(1093:1194,3),'.', 'Color', 'blue'); %TSI Toscans in Italy
% hold on;
% plot(SCORE(1195:1397,2), SCORE(1195:1397,3),'.', 'Color', 'black'); %YRI Yoruba in Nigeria
% hold on;
% plot(pc2_list(:), pc3_list(:), '.', 'Color', 'green')
% hold on;
% rectangle('Position',[RANGE_pc2_EUR(1),RANGE_pc3_EUR(1),(RANGE_pc2_EUR(2)-RANGE_pc2_EUR(1)),(RANGE_pc3_EUR(2)-RANGE_pc3_EUR(1))]);
% 
% figure;
% scatter3(SCORE(1:87,1), SCORE(1:87,2), SCORE(1:87,3), 2, 'black'); %ASW Africans in SW USA
% hold on;
% scatter3(SCORE(88:252,1), SCORE(88:252,2), SCORE(88:252,3), 2, 'blue'); %CEU N/W Europeans in Utah
% hold on;
% scatter3(SCORE(253:389,1), SCORE(253:389,2), SCORE(253:389,3), 2, 'black'); %CHB Han Chinese Beijing
% hold on;
% scatter3(SCORE(390:498,1), SCORE(390:498,2), SCORE(390:498,3), 2, 'black'); %CHD Chinese in Denver
% hold on;
% scatter3(SCORE(499:599,1), SCORE(499:599,2), SCORE(499:599,3), 2, 'black'); %GIH Gujarati in Houston
% hold on;
% scatter3(SCORE(600:712,1), SCORE(600:712,2), SCORE(600:712,3), 2, 'black'); %JPT Japanese in Tokyo
% hold on;
% scatter3(SCORE(713:822,1), SCORE(713:822,2), SCORE(713:822,3), 2, 'black'); %%LWK Luhya (Bantu) in Kenya
% hold on;
% scatter3(SCORE(823:908,1), SCORE(823:908,2), SCORE(823:908,3), 2, 'black'); %MEX Mexican in LA
% hold on;
% scatter3(SCORE(909:1092,1), SCORE(909:1092,2), SCORE(909:1092,3), 2, 'black'); %MKK Masaii in Kenya
% hold on;
% scatter3(SCORE(1093:1194,1), SCORE(1093:1194,2), SCORE(1093:1194,3), 2, 'blue'); %TSI Toscans in Italy
% hold on;
% scatter3(SCORE(1195:1397,2), SCORE(1195:1397,2), SCORE(1195:1397,3), 2, 'black'); %YRI Yoruba in Nigeria
% hold on;
% scatter3(pc1_list(:), pc2_list(:), pc3_list(:), 2, 'green')
% hold on;

