function thresh=calculate_max_connected_threshold(M)
%%
%Author: Marisa C. Ross, PhD
%Date: Feb. 2020

sort_all_zs=sort(M(:));
sort_all_zs=sort_all_zs(sort_all_zs>0);
max_z=find(sort_all_zs==max(sort_all_zs),1);
min_z=find(sort_all_zs==min(sort_all_zs),1);
thresh=sort_all_zs(round((max_z-min_z)./2));

iter=0;
done=0;
while done==0
    iter=iter+1;
    new_M=M;
    new_M(new_M<thresh)=0;
    all_strengths=strengths_und(new_M);
    find_zeros=find(all_strengths==0);
    
    if numel(find_zeros)==1
        
        max_unconnected_voxel=max(M(:,find_zeros));
        find_max_unconnect=find(sort_all_zs==max_unconnected_voxel,1);
        thresh=sort_all_zs(find_max_unconnect-1);
        done=1;
    elseif numel(find_zeros)>1
        max_z=find(sort_all_zs==thresh,1)-1;
        find_this_thresh=find(sort_all_zs==thresh,1);
        thresh=sort_all_zs(find_this_thresh-round((max_z-min_z)./2));
    elseif isempty(find_zeros)==1
        min_z=find(sort_all_zs==thresh,1)+1;
        find_this_thresh=find(sort_all_zs==thresh,1);
        thresh=sort_all_zs(find_this_thresh+ceil((max_z-min_z)./2));
    end
    if iter > 100
       done = 1;
    end
    disp(['iter = ' num2str(iter)])
    disp(['thresh = ' num2str(thresh)])
end