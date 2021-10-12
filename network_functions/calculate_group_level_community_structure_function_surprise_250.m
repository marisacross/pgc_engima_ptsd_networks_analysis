function [message result_step_four] = calculate_group_level_community_structure_function_surprise_250(atlas_location, num_iter, big_r, spatial_map_filename,mask_values)

%Author: Marisa C. Ross, PhD
%Date: Jan. 2021

disp('doing community detection')

[module_community,qual]=paco_mx(big_r,'quality',2,'nrep',num_iter);

disp(['quality of partition  ='  num2str(qual)])

%fix weird labelign that paco produces
modules=unique(module_community);
new_structure=module_community.*0;
for node=1:numel(modules)
    find_mod=find(module_community==modules(node));
    new_structure(find_mod)=node;
end
module_community=new_structure;

%loading in template mask
original=load_nii(atlas_location);
original.img = original.img(:,:,:,23);
brain_size=size(original.img);
mask=double(reshape(original.img,size(original.img,1)*size(original.img,2)*size(original.img,3),1));


new_mask_values=zeros(size(mask,1),1);

for roi=1:numel(mask_values)
    this_roi_f=find(mask==mask_values(roi));
    new_mask_values(this_roi_f,1)=module_community(roi);
end

[temp community_sort_index]=sort(module_community);

community_structure=module_community;

original.img=reshape(new_mask_values,brain_size(1),brain_size(2),brain_size(3),1);
original.hdr.dime.dim=[4 brain_size(1) brain_size(2) brain_size(3) 1 1 1 1];
original.hdr.dime.datatype=16;
original.hdr.dime.bitpix=16;
%make the filename of this an input to the function
save_nii(original,spatial_map_filename);
result_step_four.community_structure = community_structure;
result_step_four.community_sort_index = community_sort_index;
result_step_four.surprise = qual;
message=['your 3D network spatial map is located here: ' spatial_map_filename];


