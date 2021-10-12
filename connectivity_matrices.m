%% 
%Author: Marisa C. Ross, PhD
%Date: April 2021

addpath(genpath('/Volumes/Vol2/cisler/DOP/network/'))
addpath(genpath('/Volumes/Vol2/cisler/matlab_toolboxes/'))
addpath(genpath('/Volumes/Vol2/cisler/network/functions/'))
addpath(genpath('/Volumes/Vol2/cisler/Marisa/enigma_rs_networks/final_analyses/BrainNetViewer_20191031/'))
clear

%%
%define initial variables  

sub_ids = load('/Volumes/Vol2/cisler/Marisa/enigma_rs_networks/site_testing/sub_ids_all_col2');


%% Load in and Define Clinical Data
clinical_data=load(['/Volumes/Vol2/cisler/Marisa/enigma_rs_networks/site_testing/clinical_data_all_col2']);
keep_clinical=find(ismember(clinical_data(:,1),sub_ids)==1);
clinical_data=clinical_data(keep_clinical,:);

age = clinical_data(:,3); 
ptsd_severity=clinical_data(:,7); 
sex = clinical_data(:,4); %female = -.5, male = 0.5
ptsd_dx = clinical_data(:,5); %control = -0.5, ptsd = 0.5
group = clinical_data(:,6); %-1 = NTC, 0 = TEC, 1 = PTSD
index_trauma = clinical_data(:,8);
site = clinical_data(:,2);
subs = clinical_data(:,1);

%% Load in network results
load result_step_one_all_col2.mat
good_ROIs = result_step_one_surprise_250.good_ROIs;
mask_values = result_step_one_surprise_250.mask_values;

load result_step_two_all_col2.mat
group_r = result_step_two.big_r_resid;
all_sub_r = result_step_two.harmonized;

load result_step_four_all_col2.mat
community_structure = result_step_four.community_structure;
community_sort_index = result_step_four.community_sort_index;
communities_group = numel(unique(community_structure));
surprise = result_step_four.surprise;

load result_step_five_full_sample.mat
pos_par = result_step_five_surprise_sample_mask.all_graph_theory_indices(:,:,4);
mean_pos_par = mean(pos_par,2);

load all_whole_brain_strength 
load all_cen_strength
load all_dmn_strength
load mean_pos_par_dmn 
load mean_pos_par_cen
load minus_dmn_pos_par
load minus_dmn_strength
load minus_cen_pos_par
load minus_cen_strength

%% make connectivity matrices
imagesc(group_r(community_sort_index,community_sort_index),[0.4 1]);colorbar
set(gcf,'Color',[1 1 1])
set(gca,'YGrid','off')
set(gca,'Box','off')
title('Group-Level Connectivity Matrix')
print('group_level_connectivity_matrix','-dtiff')

%PTSD only
imagesc(squeeze(mean(all_sub_r(community_sort_index,community_sort_index,group==1),3)),[0.4 1]);colorbar
set(gcf,'Color',[1 1 1])
set(gca,'YGrid','off')
set(gca,'Box','off')
title('PTSD Group Connectivity Matrix')
print('ptsd_group_connectivity_matrix','-dtiff')

%TEC only
figure
imagesc(squeeze(mean(all_sub_r(community_sort_index,community_sort_index,group==0),3)),[0.4 1]);colorbar
set(gcf,'Color',[1 1 1])
set(gca,'YGrid','off')
set(gca,'Box','off')
title('TEC Group Connectivity Matrix')
print('tec_group_connectivity_matrix','-dtiff')

%NTC only
figure
imagesc(squeeze(mean(all_sub_r(community_sort_index,community_sort_index,group==-1),3)),[0.4 1]);colorbar
set(gcf,'Color',[1 1 1])
set(gca,'YGrid','off')
set(gca,'Box','off')
title('NTC Group Connectivity Matrix')
print('ntc_group_connectivity_matrix','-dtiff')

