function [effect_position_mat,unique_effects]=effect_position_finder(effect_list)

    %effect_list={'a','b','d','a','d','d'};
    %effect_list = [1 1 3 5 4 5 5 2];
%    [unique_effects_sorted,unique_idx,col_position_idx]=unique(effect_list,'first');
    [unique_effects,~,column_positions]=unique(effect_list,'first');
%    [~,unique_order]=sort(unique_idx);
    
%    unique_effects = unique_effects_sorted(unique_order);
%    unique_effect_positions=unique_order(unique_order);
%    column_positions=unique_idx(col_position_idx);

    effect_num=length(unique_effects);
    row_number=length(effect_list);

%    unique_effect_mat=repmat(unique_effects,[row_number,1]);
%    effect_list_mat=repmat(effect_list',[1,effect_num]);
    
    effect_position_mat=sparse(1:row_number,column_positions,...
        ones([1 row_number]),row_number,effect_num);

    %effect_position_mat=(unique_effect_mat==effect_list_mat)
%    if iscell(unique_effect_mat)
%        effect_position_mat=sparse(strcmp(unique_effect_mat,effect_list_mat));
%    else
%        effect_position_mat=sparse((unique_effect_mat==effect_list_mat));
%    end

end