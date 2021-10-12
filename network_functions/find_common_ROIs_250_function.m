function result_step_one =find_common_ROIs_250_function(sub_ids,atlas_location,resting_filename)

%The purpose of this function is to use a pre-determined atlas and subject
%data to find common ROIs across subjects and remove uncommon or junk ROIs.
%
%Inputs: 
%    sub_ids = subject IDs
%    atlas_location = path to atlas of choice
%    resting_filename = path to subject's resting state file 
%
%Outputs: result_step_one, a structure containing:
%    good_ROIs = matrix of subject ROIs following removal of noise/bad ROIs
%    mask_values = unique ROIs in the subjects' brains from within the
%    original atlas
%
% "Happiness can be found even in the darkest of times if one simply remembers to
% turn on the light" - Albus Dumbledore
%Author: Marisa C. Ross July 2019

%load in atlas
original_atlas=load_nii(atlas_location);
brain_size=size(original_atlas.img);
original_atlas=double(reshape(original_atlas.img,size(original_atlas.img,1)*size(original_atlas.img,2)*size(original_atlas.img,3),1));
in_brain=find(original_atlas~=0);
mask_values=unique(original_atlas(in_brain)); % mask values = unique ROIs in the brain within the original atlas
initial_ROIs=1:numel(mask_values);

%preallocating variables
track_bad_ROIs=[];
all_disconnected=[];
all_nans=[];
%looping through each number of subjects, and calling the looping variable
%sn
for sn=1:numel(sub_ids)
    
    disp(sub_ids(sn));
    eval('!rm temp_c.nii')
    
    eval(['!3dAFNItoNIFTI -prefix ./temp_c ' char(resting_filename(sn))]);
  
    original=load_nii(['./temp_c.nii']);
    eval('!rm temp_c.nii')
    sub_data=reshape(original.img,brain_size(1)*brain_size(2)*brain_size(3),size(original.img,4));
    sub_data=double(sub_data(in_brain,:)');
    clear original
    
    %account for fact that fourier transform results in first two
    %timepoints being identical
    sub_data=sub_data(2:end,:);
   
    roi_data=zeros(size(sub_data,1),numel(mask_values));
    for rois=1:numel(mask_values)
        this_roi=find(original_atlas(in_brain)==mask_values(rois));
        find_non_zero=find(sum(sub_data(:,this_roi))~=0);

        if numel(find_non_zero)>10

            roi_data(:,rois)=zscore(mean(sub_data(:,this_roi(find_non_zero)),2));
        end

    end
    r=corrcoef(roi_data);
    

    not_invariant=find(var(roi_data)~=0);
    invariant=setdiff(initial_ROIs, not_invariant);
    bad_ROIs(sn,1) = numel(invariant);
    track_bad_ROIs=[track_bad_ROIs setdiff(initial_ROIs,initial_ROIs(not_invariant))];
end
track_bad_ROIs=unique(track_bad_ROIs);
good_ROIs=setdiff(1:numel(mask_values),track_bad_ROIs);
mask_values=mask_values(good_ROIs);

result_step_one.good_ROIs = good_ROIs;
result_step_one.mask_values = mask_values;
result_step_one.bad_ROIs_sub = bad_ROIs;



