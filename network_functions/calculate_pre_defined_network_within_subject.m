function result = calculate_pre_defined_network_within_subject(atlas_location,sub_num, pre_defined_network_location,subject_spatial_map_output_prefix,community_structure,mask_values,all_sub_r)
% Use this function to generate the relevant community structure for each
% subject, given predefined atlas 
%
%Inputs: 
%   atlas_location = path to atlas of choice
%   sub_num = subject IDs
%   subject_spatial_map_output_prefix = user-defined filename prefix for
%   individual spatial map 
%   community_structure = best structure of network organization as
%   defied by the predefined network
%   mask_values = unique ROIs in the subjects' brains from within the
%   original atlas
%   all_sub_r = saving all subject correlation matrices

% Outputs: 
%   result_predefined_within, a structure containing:
%       all_graph_theory_indices = saved matrix of graph theory metrics,
%       including module_community, roi_strength, roi_clust, roi_participation, roi_connector_hubs, 
%       roi_within_z and roi_betweeness
%   
% "Words, in my not-so-humble opinion, are our most inexhaustible source of
% magic" - Albus Dumbledore
%-Marisa Ross & Josh Cisler - 7-25-2019


for sn=1:numel(sub_num)   
    %loading in atlas
    original_atlas=load_nii(atlas_location);
    brain_size=size(original_atlas.img);
    original_atlas=double(reshape(original_atlas.img,brain_size(1)*brain_size(2)*brain_size(3),1));

    %load in pre-defined network
    pre_defined_network=load_nii(pre_defined_network_location);
    pre_defined_network=reshape(pre_defined_network.img,brain_size(1)*brain_size(2)*brain_size(3),1);

    community_structure=zeros(numel(mask_values),1);
    for roi=1:numel(mask_values)
        this_roi_f=find(original_atlas==mask_values(roi));
        find_not_zero=find(pre_defined_network(this_roi_f)>0);
        community_structure(roi,1)=mode(pre_defined_network(this_roi_f(find_not_zero)));
    end


    [temp community_sort_index]=sort(community_structure);


    result.community_structure = community_structure;
    result.community_sort_index = community_sort_index;
    %the yeo network structure is actually smaller (less GM) than the 1000
    %atlas, so this can result in NANs in the community structure. account for
    %that here and create a reduced mask_values and good_ROIs if so
%     if sum(isnan(community_structure))>0
%         f_good_ROIs=find(isnan(community_structure)==0);
%         result.new_mask_values=mask_values(f_good_ROIs);
%         result.new_good_ROIs=good_ROIs(f_good_ROIs);
%         %restrict community structure to remaining good ROIs
%         result.community_structure=community_structure(f_good_ROIs);
%         %now have to resort with the truncated ROIs
%         [temp community_sort_index]=sort(community_structure(f_good_ROIs));
% 
%         result.community_sort_index=community_sort_index;
%     end

        %%
    %doing other graph theory stuff
    
    %calculating node statistics for strength, clustering coefficient,
    %participation, identifying 'connector' hubs, and within module
    %strength
    disp('doing graph theory calcs')
    roi_strength=strengths_und(squeeze(all_sub_r(:,:,sn)));
    roi_clust=clustering_coef_wu_sign(squeeze(all_sub_r(:,:,sn)),3);
    [pos_par neg_par]=participation_coef_sign(squeeze(all_sub_r(:,:,sn)),community_structure);
    roi_within_z=module_degree_zscore(squeeze(all_sub_r(:,:,sn)),community_structure);
    %scale z to range 0-1, then take inverse in order to create length
    %matrix
   % z_scaled=(z-min((min(z))))/(max(max(z))-min(min(z)));
    %z_length=1-z_scaled;
    %roi_betweeness=betweenness_wei(z_length);
    roi_betweeness=betweenness_wei(squeeze(all_sub_r(:,:,sn)));
    
    %stack indices together for an easier loop later
    all_graph_theory_indices=[community_structure roi_strength' roi_clust pos_par neg_par...
        roi_within_z roi_betweeness];
    
    %put back in original MNI atlas space
    original=load_nii(atlas_location);
    brain_size=size(original.img);
    mask=double(reshape(original.img,size(original.img,1)*size(original.img,2)*size(original.img,3),1));
    new_mask_values=zeros(size(mask,1),size(all_graph_theory_indices,2));
    for roi=1:numel(mask_values)
        this_roi_f=find(mask==mask_values(roi));
        for indices=1:size(all_graph_theory_indices,2);
            new_mask_values(this_roi_f,indices)=all_graph_theory_indices(roi,indices);
        end
    end

    original.img=reshape(new_mask_values,size(original.img,1),size(original.img,2),size(original.img,3),size(all_graph_theory_indices,2));
    original.hdr.dime.dim=[4 brain_size(1) brain_size(2) brain_size(3) size(all_graph_theory_indices,2) 1 1 1];
    original.hdr.dime.datatype=16;
    original.hdr.dime.bitpix=16;
    save_nii(original,[subject_spatial_map_output_prefix sub_num{sn} '.nii']);
    
    graph_theory_indices(sn,:,:) = all_graph_theory_indices;
    result.all_graph_theory_indices = graph_theory_indices;

end



