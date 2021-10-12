function result_step_two_combat = stack_sub_matrices_function_combat_enigma(sub_ids,good_ROIs, resting_filename, atlas_location,mask_values,site,covariates)

%The purpose of this script is to create large matrices of subject data do
%eventually do community detection. This function allows for the inclusion
%of negative weights by not zeroing out negative correlations in the
%matricies
%
%Inputs:
%   sub_ids = subject IDs
%   good_ROIs = matrix of subject ROIs following removal of noise/bad ROIs
%   resting_filename = path to subject's resting state file 
%   atlas_location = path to atlas of choice
%   site = vector representing scanner site
%
%Outputs: result_step_two, a structure containing:
%   stacked_sub_matrices = large matrix of individual subject time course
%   data, concatinated across ROIs
%   big_r = group-level correlation matrix with site variance regressed out
%   all_sub_r = saving all subject correlation matrices, site variance
%   regressed out
%   save all_TRs = sum of TRs loaded in from dataset 
%   save_all_FD =  saves the mean of framewise displacement across the
%   collected motion measures
%
%Author: Marisa C. Ross, PhD Feb. 2021

%load in atlas
original_atlas=load_nii(atlas_location);
brain_size=size(original_atlas.img);
original_atlas=double(reshape(original_atlas.img,size(original_atlas.img,1)*size(original_atlas.img,2)*size(original_atlas.img,3),1));
in_brain=find(original_atlas~=0); %creates variable storing only the voxels with non-zero values (showing brain) from original_atlas 

stacked_sns=[];
stacked_sns_resid = [];
stacked_sub_matrices_resid=[];
stacked_sub_matrices=[];

for sn=1:numel(sub_ids)
    disp(sub_ids(sn));
    eval('!rm temp_c.nii')      %removes file from foldersave good_ROIs_DOP good_ROIs, if it exists, that is
    
    eval(['!3dAFNItoNIFTI -prefix ./temp_c ' char(resting_filename(sn))]);        %converts AFNI file to 3d NIFTI file, opens file
  
    original=load_nii(['./temp_c.nii']);       %loads NIFTI dataset
    eval('!rm temp_c.nii')      %removes file from folder
    sub_data=reshape(original.img,brain_size(1)*brain_size(2)*brain_size(3),size(original.img,4));  %reshapes the 4-D (3 spatial dimensions plus time) into 2-D matrix (space x time)
    sub_data=double(sub_data(in_brain,:)'); %converts the matrix into a "double precision" array, and restricts to just voxels in brain
    clear original
     
    %account for fact that fourier transform results in first two
    %timepoints being identical
    sub_data=sub_data(2:end,:);
      
    
    roi_data=zeros(size(sub_data,1),numel(mask_values));    %creates a 2D array of zeros the size of sub_datax1xmask_values, preallocating a variable
    roi_data_no_censor=zeros(size(sub_data,1),numel(mask_values));
    for rois=1:numel(mask_values)       %loops through all ROIs on mask
        this_roi=find(original_atlas(in_brain)==mask_values(rois));     %looping through each roi in the mask and finding the voxels in that ROI in 1d space
        %find voxels whos sum is not 0, which means there is actually 
        %data in the voxel. this is excluding voxels that are outside the participant's brain
        find_non_zero=find(sum(sub_data(:,this_roi))~=0);

        if numel(find_non_zero)>10%only calculate mean 

            roi_data(:,rois)=zscore(mean(sub_data(:,this_roi(find_non_zero)),2)); %creates a coefficient 
            %using SD of population flagged as '2'
            roi_data_no_censor(:,rois)=zscore(mean(sub_data(:,this_roi(find_non_zero)),2));
        end

    end
    
    stacked_sub_matrices=[stacked_sub_matrices;roi_data]; 
    stacked_sns=[stacked_sns;ones(size(roi_data,1),1).*sn]; %concatenates subject data across pp

    r=corrcoef(roi_data); %determines correlation coefficient of roi data across time
    find_nans=find(isnan(r)==1); %finds diagonals (correlation of 1)
    r(find_nans)=0; %sets diagonal value to zero
    
    for rd=1:size(r,1)
        r(rd,rd)=0;
    end
    r=.5*log((1+r)./(1-r));
    all_sub_r(:,:,sn)=r;   %saving all correlation matrices across subjects
    result_step_two_combat.all_sub_ROI_data{sn}=roi_data;
    result_step_two_combat.all_sub_ROI_data_no_censor{sn}=roi_data_no_censor;
end   
save all_sub_r all_sub_r

big_r=corrcoef(stacked_sub_matrices);
for rd=1:size(big_r,1)
    big_r(rd,rd)=0;
end

big_r=.5*log((1+big_r)./(1-big_r));
    

%get ROI-ROI FC to put into a matrix to regress out site using ComBat
    
    %create covariates in matrix, include age, sex, group,
    %etc. in a variable called all_x
    mod = zscore(covariates);
    batch = [site]; 
    
    for sub = 1:numel(sub_ids)
        C = all_sub_r(:,:,sub);
        dat(:,sub) = C(triu(true(size(C)),1)); %isolate upper triangle, put into feature x n matrix
    end
 
    bayesdata = combat(dat, batch, mod, 1); %run combat on upper triangle
    save bayesdata bayesdata

for sn = 1:numel(sub_ids) %reshape vector of upper triangle back into square matrix
    p = bayesdata(:,sn)';
    a = triu(ones(size(all_sub_r,1)),1); %creates upper triangle in the size of the square matrix, puts ones in upper, zeros in lower
    a(a>0)=p; %fills in 1s with appropriate values from bayesdata
    all_bayes(:,:,sn) = (a + a') ./ (eye(size(all_sub_r,1))+1); %sets the zeros on the diagonal 

end
    
    big_r_resid = mean(all_bayes,3); 


    
result_step_two_combat.big_r=big_r;
result_step_two_combat.big_r_resid=big_r_resid;
result_step_two_combat.all_sub_r=all_sub_r;
result_step_two_combat.harmonized = all_bayes;
result_step_two_combat.stacked_sub_matrices=stacked_sub_matrices;
result_step_two_combat.stacked_sub_matrices_resid=stacked_sub_matrices_resid;
result_step_two_combat.stacked_sub_nums=stacked_sns;
result_step_two_combat.stacked_sub_nums_resid=stacked_sns_resid;












