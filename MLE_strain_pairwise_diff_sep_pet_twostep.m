function MLE_strain_pairwise_diff_sep_pet_twostep(key_list, value_list)

    % EP 18-08-03

    % Master function that runs MLE_strain_pairwise_diff_sep_pet, but
    	% does this in two steps: first with high tolerance values
    	% and using a small subset of the data, then with the full data
    	% and higher tolerance, but starting from the local solutions
    	% of the previous run

    parameter_dict = containers.Map(key_list,value_list);
    datafile_path = parameter_dict('datafile_path');
    output_id_parameter = parameter_dict('output_id_parameter');
    external_counter = str2double(parameter_dict('external_counter'))

    ms_starting_point_file = fullfile(datafile_path,strcat('ms_starting_points-',...
        output_id_parameter,'_',int2str(external_counter),'.csv'))

    initial_data_fraction = 0.2;
    tolx_val = 10^-3;
    tolfun_val = 10^-2;
    write_solutions_output = true;

    % run initial mle with subset of data
    key_list_prerun = [key_list {'ms_starting_point_file','initial_data_fraction','tolx_val','tolfun_val','write_solutions_output'}];
    value_list_prerun = [value_list {ms_starting_point_file,initial_data_fraction,tolx_val,tolfun_val,write_solutions_output}];

    MLE_strain_pairwise_diff_sep_pet(key_list_prerun, value_list_prerun);

    % run final mle with default parameters
	MLE_strain_pairwise_diff_sep_pet(key_list, value_list);    

