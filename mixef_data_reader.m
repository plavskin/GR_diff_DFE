function [response_vector,fixef_ID_mat,ranef_corr_struct,unique_fixefs,...
    unique_ranef_categories,order_vector,block_start_positions] = ...
    mixef_data_reader(data_file,response_name,...
    fixef_name,random_effect_names,block_effect_name,prop_data_to_use)

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_random_effects = length(random_effect_names);
        % number of random effects including epsilon

    %data_read_string = strcat('%f %q',repmat(' %q',[1 num_random_effects]),'\n');

    %fid=fopen(data_file);
    %data_read = textscan(fid, ...
    %data_read_string,'Delimiter',',');
    %fclose(fid);

    data_table = readtable(data_file);
    indices_to_use = logical(binornd(1,prop_data_to_use,[1 size(data_table,1)]));
    data_table = data_table(indices_to_use,:);
    
%    response_vector = data_read{1}';
    response_vector = data_table{:,response_name};

    % make correlation structure for fixed effect and identify position of intercept effect, if such exists
%    fixef_vector=data_read{2}';
    fixef_vector = data_table{:,fixef_name};
    [fixef_ID_mat,unique_fixefs] = effect_position_finder(fixef_vector);
    
    unique_ranef_categories = zeros(1,num_random_effects);
    ranef_corr_struct = struct;

    if isnan(block_effect_name)
        order_vector = NaN;
        block_start_positions = NaN;
    else
        order_vector = NaN([1,size(data_table,1)]);
    end

    for current_ranef_idx = 1:num_random_effects
        current_ranef_name = random_effect_names{current_ranef_idx};
%        current_ranef_vector = data_read{2+current_ranef_idx}';
        current_ranef_vector = data_table{:,current_ranef_name};
        [current_ID_mat,~] = effect_position_finder(current_ranef_vector);
        ranef_corr_struct.(current_ranef_name) = current_ID_mat;
        unique_ranef_categories(current_ranef_idx) = size(current_ID_mat,2);

        % if current ranef corresponds to the known block-diagonal
            % structure, get info on block-diagonal structure from it
        if strcmp(current_ranef_name,block_effect_name)
            % for each instantiation of the current random effect, make
                % list of indices which correspond to that instantiation
                % (i.e. which come from experimental plate
                % #current_column); concatenate these in order to create
                % order_vector, which sorts the data by the current ranef
            current_ranef_instantiations = size(current_ID_mat,2);
            block_start_positions = NaN([1,current_ranef_instantiations]);
            next_start_position = 1;
            for current_column = 1:current_ranef_instantiations
                current_start_position = next_start_position;
                current_order_vector = find(current_ID_mat(:,current_column));
                current_element_num = length(current_order_vector);
                next_start_position = current_element_num+current_start_position;
                order_vector(current_start_position:(next_start_position-1)) = current_order_vector;
                block_start_positions(current_column) = current_start_position;
            end

        end

    end

%    if ~(strcmp(epsilon_effect_name,'NA') || strcmp(epsilon_effect_name,'NaN'))
%        ranef_corr_struct.(epsilon_effect_name) = spdiags(ones([number_datapoints,1]),0,number_datapoints,number_datapoints);
%        unique_ranef_categories = [unique_ranef_categories,number_datapoints];
%    end

end