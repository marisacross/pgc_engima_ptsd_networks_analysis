function result = group_level_node_analysis(all_graph_theory_indices,atlas_location,covariates,clinical_variables,clinical_labels,mask_values,group_spatial_map_output_prefix,p_crit,stat_choice)

% The purpose of this function is:
% 1) Do permutation testing to determine critical threshold test value (either
%    t-value from regression or rho value from partial correlation) and 
% 2) Perform a between-group statistical test (either regression or partial correlation) on
%    each node/ROI within a functional brain map for each graph theory
%    index. 
% 
% Inputs:
%   all_graph_theory_indices: sub_num x roi x graph theory index matrix
%   (where your data live)
%   atlas_location: original atlas of group-level ROIs
%   covariates: vector of covariates for stats testing
%   clinical_labels: string containing the name of your clinical variable
%   mask_values: vector containing node assignments for each ROI
%   spatial_map_output_prefix: name for eventual group-level spatial map
%   p_crit: critical p-value threshold (i.e. 0.05)
%   stat_choice: takes values of 1 or 2. 1 = partial correlation using
%   spearman; 2 = linear regression using regstats
% Output is a result structure containing:
%   save_crit_values: critical values for each graph theory index
%   save_stat: matrix of group-level rho or t value for each ROI for each graph theory index 
%   save_ps: matrix of group-level p values for each ROI for each graph theory index
%   spatial map: nii file containing group-level spatial map with t or rho
%   values for each ROI
% Authors: Marisa C. Ross, PhD & Josh M. Cisler, PhD July 2019


%% permutation testing to define critical test value level
p=gcp('nocreate');
if isempty(p)==1
    parpool
end
iters=5000;
pctRunOnAll warning('off','all')
stat_crit_values=zeros(size(clinical_variables,2),size(all_graph_theory_indices,3));
for clin=1:size(clinical_variables,2)
    
    for index=1:size(all_graph_theory_indices,3)
        disp(['permutation testing for index ' num2str(index)])
        crit_value = perm_group_testing(all_graph_theory_indices(:,:,index),clinical_variables(:,clin),covariates,iters,p_crit,stat_choice);
        stat_crit_values(clin,index)=crit_value;
    end
end
result.save_crit_values=stat_crit_values;

save save_crit_values stat_crit_values

disp('finished permutation testing, doing actual testing....')


%% do actual testing and spit out maps
save_stat=zeros(size(all_graph_theory_indices,2),size(all_graph_theory_indices,3),size(clinical_variables,2));
save_ps=save_stat;

%original atlas for indicing purposes
original=load_nii(atlas_location);
warning off
for clin=1:size(clinical_variables,2)    
    for roi=1:size(all_graph_theory_indices,2)
        for index=1:size(all_graph_theory_indices,3)
            if stat_choice == 1
                [rho p]=partialcorr([clinical_variables(:,clin) all_graph_theory_indices(:,roi,index)],covariates,'type','spearman');
                save_stat(roi,index,clin)=rho(1,2);
                save_ps(roi,index,clin)=p(1,2);
            elseif stat_choice == 2
                [stats] = regstats(all_graph_theory_indices(:,roi,index),[clinical_variables(:,clin) covariates]);
                save_stat(roi,index,clin)=stats.tstat.t(2);
                save_ps(roi,index,clin)=stats.tstat.pval(2);
            end
        end
    end
 
    brain_size=size(original.img);
    mask=double(reshape(original.img,size(original.img,1)*size(original.img,2)*size(original.img,3),1));
    new_mask_values=zeros(size(mask,1),size(all_graph_theory_indices,3));
    for roi=1:numel(mask_values)
        this_roi_f=find(mask==mask_values(roi));
        for indices=1:size(all_graph_theory_indices,3);
            new_mask_values(this_roi_f,indices)=save_stat(roi,indices);

        end
    end

    original.img=reshape(new_mask_values,size(original.img,1),size(original.img,2),size(original.img,3),size(all_graph_theory_indices,3));
    original.hdr.dime.dim=[4 brain_size(1) brain_size(2) brain_size(3) size(all_graph_theory_indices,3) 1 1 1];
    original.hdr.dime.datatype=16;
    original.hdr.dime.bitpix=16;
    save_nii(original,[group_spatial_map_output_prefix clinical_labels{clin} '.nii']);
end
result.save_stat=save_stat;
result.save_ps=save_ps;