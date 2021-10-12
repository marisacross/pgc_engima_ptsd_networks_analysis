function crit_value = perm_group_testing_stat_choice(data,clinical_variable,covariates,iters,p_crit,stat_choice)
%%
%Author: Marisa C. Ross, PhD
%Date: Feb. 2020    
   
feature_size=size(data,2);
max_t=zeros(iters,1);
    parfor i=1:iters
        r=randperm(numel(clinical_variable));
        perm_variable=clinical_variable(r);
        test_values=zeros(feature_size,1);
        for roi=1:feature_size
            if stat_choice == 1
                rho=partialcorr([perm_variable data(:,roi)],covariates,'type','spearman');
                test_values(roi)=rho(1,2);
            elseif stat_choice == 2
                [stats] = regstats(perm_variable, [data(:,roi) covariates]);
                test_values(roi) = stats.tstat.t(2)
            end
        end
        max_t(i,1)=max(test_values);
    end
    sort_all_values=sort(max_t,'descend');
    crit_value=sort_all_values(round(numel(sort_all_values).*p_crit));
end