%%
%Author: Marisa C. Ross, PhD
%Date: March 2021

addpath(genpath('/Volumes/Vol2/cisler/DOP/network/'))
addpath(genpath('/Volumes/Vol2/cisler/matlab_toolboxes/'))
addpath(genpath('/Volumes/Vol2/cisler/network/functions/'))
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
all_graph_indices = result_step_five_surprise_sample_mask.all_graph_theory_indices;
load all_cen_strength.mat

%% ROI-by-ROI break down of graph-theory indices

% First, make .nii of good_ROIs so number of nodes is accurate

ROIs = load_nii('/Volumes/Vol2/cisler/Marisa/enigma_rs_networks/final_analyses/good_ROIs_221.nii');
ROI_size = size(ROIs.img);
in_brain = find(ROIs.img > 0);

store_all_pos_par=zeros(numel(sub_ids),numel(in_brain));

for sn=1:numel(sub_ids)
    if sub_ids(sn)<10
       sub_num=['000' num2str(sub_ids(sn))];
    elseif sub_ids(sn)>=10&&sub_ids(sn)<100
       sub_num=['00' num2str(sub_ids(sn))];
    elseif sub_ids(sn)>=100&&sub_ids(sn)<1000
       sub_num=['0' num2str(sub_ids(sn))];
    elseif sub_ids(sn)>1000
        sub_num=num2str(sub_ids(sn));
    end
    disp(sub_num)
    data=load_nii(['/Volumes/Vol2/cisler/Marisa/enigma_rs_networks/final_analyses/sub_maps/all_subs_' sub_num '.nii']);
    sub_brain=double(reshape(squeeze(data.img),ROI_size(1)*ROI_size(2)*ROI_size(3),size(squeeze(data.img),4)));
    select_briks = size(sub_brain,3):size(sub_brain,2);
    sub_briks = sub_brain(in_brain,select_briks);
    store_all_briks(sn,:,:,:) = sub_briks;
    sub_brain_pos_par=sub_brain(in_brain,4); % this saves the brik containing the pos_par for every subject, for every node, in a vector
    store_all_pos_par(sn,:)=sub_brain_pos_par;
end

save store_all_pos_par store_all_pos_par
%to use whole_brain_regression function
%make design matrix
X=[ptsd_dx sex age site];
filename='ptsd_v_controls_t_';
whole_brain_data=store_all_pos_par;%data you care about in participant x brain space
mask=in_brain;
%data is an original nifti file loaded in and untouched
message = whole_brain_regression(whole_brain_data,X,filename,mask,data);

% Node analysis at the group level to determine critical t-values
% load in all relevant inputs

all_graph_theory_indices_resid=zeros(numel(sub_ids),numel(mask_values),7);
for l=1:numel(mask_values)
    find_roi=find(ROIs.img(in_brain)==mask_values(l));
    temp=squeeze(store_all_briks(:,find_roi,:));
    all_graph_theory_indices_resid(:,l,:) = squeeze(temp(:,1,:));
end
save all_graph_theory_indices_mask all_graph_theory_indices_resid

atlas_location = '/Volumes/Vol2/cisler/Marisa/enigma_rs_networks/final_analyses/good_ROIs_221.nii';
covariates = [age sex site];
clinical_variables = ptsd_dx;
clinical_labels = {'ptsd'};
group_spatial_map_output_prefix = '/Volumes/Vol2/cisler/Marisa/enigma_rs_networks/final_analyses/ptsd_v_controls_map_';
p_crit = 0.05;
stat_choice = 2; %stat_choice 1 = partial correlation (spearman), stat_choice 2 = regression (tstat)

result = group_level_node_analysis(all_graph_theory_indices_resid,atlas_location,covariates,clinical_variables,clinical_labels,mask_values,group_spatial_map_output_prefix,p_crit,stat_choice);

save ptsd_v_controls_node_result result


%% Node-level group analysis

load ptsd_v_controls_node_result.mat
critcal_values = result.save_crit_values;
t_values = result.save_stat;
p_values = result.save_ps;



