% Input script
data_QBS_LHC;

% version
version = 2;

variable_list = {
        'Cell_list',
        'Type_list',
        'Sector_list',
        'EH84x_list',
        'TT84x_list',
        'CV94x_list',
        'PT961_list',
        'PT991_list',
	'TT94x_list',
        'TT961_list',
        'R_list',
        'Qs_list',
        'Kv_list',
        'nc_list',
        'L_list'
        };

filename = sprintf('data_qbs_lhc_%i.csv', version);
output_file = fopen(filename, 'w');

for i = 1:length(variable_list)
	if i == 1
		fprintf(output_file,'%s', variable_list{i});
	else
		fprintf(output_file,'\t%s',variable_list{i});
	end
end
fprintf(output_file,'\n');

for i = 1:length(Cell_list)
	fprintf(output_file, '%s', Cell_list{i});
	fprintf(output_file, '\t%s', Type_list{i});
	fprintf(output_file, '\t%s', Sector_list{i});
	fprintf(output_file, '\t%s', EH84x_list{i});
	fprintf(output_file, '\t%s', TT84x_list{i});
	fprintf(output_file, '\t%s', CV94x_list{i});
	fprintf(output_file, '\t%s', PT961_list{i});
	fprintf(output_file, '\t%s', PT991_list{i});
	fprintf(output_file, '\t%s', TT94x_list{i});
	fprintf(output_file, '\t%s', TT961_list{i});
	fprintf(output_file, '\t%e', R_list(i));
	fprintf(output_file, '\t%e', Qs_list(i));
	fprintf(output_file, '\t%e', Kv_list(i));
	fprintf(output_file, '\t%e', nc_list(i));
	fprintf(output_file, '\t%e', L_list(i));
	fprintf(output_file,'\n');
end

fclose(output_file);
