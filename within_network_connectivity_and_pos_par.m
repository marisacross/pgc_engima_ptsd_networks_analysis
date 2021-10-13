%%
%Author: Marisa C. Ross, PhD
%Date: May 2021

% Use this script to calculate graph network measures for individual networks of choice (this example uses the DMN and CEN)
clear

%%
%define initial variables  

sub_ids = load('/path/to/subject_list');


%% Load in and Define Clinical Data
clinical_data=load(['/path/to/clinical_data']);
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
load result_step_one.mat
good_ROIs = result_step_one.good_ROIs;
mask_values = result_step_one.mask_values;

load result_step_two.mat
group_r = result_step_two.big_r_resid;
all_sub_r = result_step_two.harmonized;

load result_step_four.mat
community_structure = result_step_four.community_structure;
community_sort_index = result_step_four.community_sort_index;
communities_group = numel(unique(community_structure));
surprise = result_step_four.surprise;

load result_step_five.mat
pos_par = result_step_five.all_graph_theory_indices(:,:,4);
mean_pos_par = mean(pos_par,2);
%% use group-level corr matrices to isolate upper triangle and calculate mean correlation strength

M2 = group_r(triu(true(size(group_r)),1)); % to get rid of zeros along diagonal, make into a vector where each element
% represents connectivity of that node to the rest of the brain
% size = 24310 x 1
whole_brain_strength = mean(M2); %average connectivity strength for sample across the brain


%% subject-level correlation matrices, average strength
% 
num_rois=numel(community_sort_index);
all_whole_brain_strength=zeros(numel(subs),1);
for sn = 1:size(all_sub_r,3)
    C = all_sub_r(:,:,sn);
    all_sub_r_v = C(triu(true(size(C)),1)); %isolate upper triangle for each sub, put in vector 
    % where each element represents connectivity strength of that node to the rest of the brain
    % size = 24310 x 1063
    all_whole_brain_strength(sn,1)=mean(all_sub_r_v,1); 
    %strength for each subject; size = 1063 x 1                    
end
    
save all_whole_brain_strength all_whole_brain_strength
% 
%% within-networks

%define the canonical networks based on sample spatial map
base=load_nii(['group_level_structure.nii']);
brain=double(base.img);
brain_size=size(brain);
max_modules=max(max(max(brain)));

DMN=[2];
CEN =[6];

% this finds the voxels (not ROIs) associated with each network
find_dmn = find(ismember(brain,DMN)==1);
find_cen = find(ismember(brain,CEN)==1);
all_minus_dmn = find(ismember(brain,DMN)~=1);
all_minus_cen = find(ismember(brain,CEN)~=1);

atlas=load_nii(['good_ROIs_221.nii']);
dmn_rois_in_atlas=unique(atlas.img(find_dmn));
cen_rois_in_atlas=unique(atlas.img(find_cen));

minus_dmn_rois_in_atlas=unique(atlas.img(all_minus_dmn)); 
minus_cen_rois_in_atlas=unique(atlas.img(all_minus_cen));

%get rid of the zero value
non_zero_dmn = find(minus_dmn_rois_in_atlas>0);
no_zero_cen = find(minus_cen_rois_in_atlas>0);
minus_dmn_rois_in_atlas = minus_dmn_rois_in_atlas(non_zero_dmn);
minus_cen_rois_in_atlas = minus_cen_rois_in_atlas(no_zero_cen);
% 
%find in your list of ROIs that are used where in that
%list these R0Is are. That index is what is used in
%all_graph_theory_indices
% mask_values contains the labels of the good_ROIs in
% all_graph_theory_indices

dmn_index = find(ismember(mask_values,dmn_rois_in_atlas)==1);
cen_index = find(ismember(mask_values,cen_rois_in_atlas)==1);
all_minus_dmn_index = find(ismember(mask_values,minus_dmn_rois_in_atlas)==1);
all_minus_cen_index = find(ismember(mask_values,minus_cen_rois_in_atlas)==1);
% 
dmn_strength = all_sub_r(dmn_index,dmn_index,:); %size = 28 x 28 x1063
cen_strength = all_sub_r(cen_index,cen_index,:); %size = 23x23x1063
all_minus_dmn_strength = all_sub_r(all_minus_dmn_index,all_minus_dmn_index,:); %size = 193 x 193 x1063
all_minus_cen_strength = all_sub_r(all_minus_cen_index,all_minus_cen_index,:); %size = 198 x 198 x1063

all_cen_strength=zeros(numel(subs),1);
all_dmn_strength=zeros(numel(subs),1);
minus_dmn_strength=zeros(numel(subs),1);
minus_cen_strength=zeros(numel(subs),1);
for sn = 1:size(all_sub_r,3)
    C = dmn_strength(:,:,sn);
    all_dmn_r = C(triu(true(size(C)),1));
    % size = 378x1063; represents the connectivity of each node in the DMN
    % to other nodes in the DMN
    all_dmn_strength(sn,1)=mean(all_dmn_r,1);  
    % size = 1063x1
    D = cen_strength(:,:,sn);
    all_cen_r = D(triu(true(size(D)),1));
    % size = 253x1063
    all_cen_strength(sn,1)=mean(all_cen_r,1); 
    E = all_minus_dmn_strength(:,:,sn);
    all_minus_dmn_r = E(triu(true(size(E)),1));
    minus_dmn_strength(sn,1)=mean(all_minus_dmn_r,1);  
    F = all_minus_cen_strength(:,:,sn);
    all_minus_cen_r = F(triu(true(size(F)),1));
    minus_cen_strength(sn,1)=mean(all_minus_cen_r,1); 
end

save all_dmn_strength all_dmn_strength
save all_cen_strength all_cen_strength
save minus_dmn_strength minus_dmn_strength
save minus_cen_strength minus_cen_strength

%% Determine average participation coefficient for each network

all_pos_par = result_step_five_surprise_sample_mask.all_graph_theory_indices(:,:,4);

% Next, find the pos_par values for each ROI within each canonical network
DMN_pos_par = all_pos_par(:,dmn_index);%this is a 1063x28 matrix with the pos_par of each ROI within the DMN for each subject
mean_pos_par_dmn = mean(DMN_pos_par,2);
minus_dmn_pos_par = all_pos_par(:,all_minus_dmn_index);
minus_dmn_pos_par = mean(minus_dmn_pos_par,2);

save mean_pos_par_dmn mean_pos_par_dmn
save minus_dmn_pos_par minus_dmn_pos_par

CEN_pos_par = all_pos_par(:,cen_index); %this is a 1063x23 matrix with the pos_par of each ROI within the CEN for each subject
mean_pos_par_cen = mean(CEN_pos_par,2);
minus_cen_pos_par = all_pos_par(:,all_minus_cen_index);
minus_cen_pos_par = mean(minus_cen_pos_par,2);

save mean_pos_par_cen mean_pos_par_cen
save minus_cen_pos_par minus_cen_pos_par

