function invoke_matrix_outputting(matrix, matrixname, rowLabels, columnLabels, some_comment)

base_path = pwd;

current_time = datestr(now);
ct1 = strrep(current_time,':','_'); % substitute : with _
ct2 = strrep(ct1,' ','_'); % substitute spaces with _
some_comment = strrep(some_comment,':','_');
some_comment = strrep(some_comment,' ','_');

final_name = [matrixname '_' some_comment '_' ct2 '.tex'];

if exist([base_path '\Latex'], 'dir')
   addpath([base_path '\Latex']) 
   path.out = ['Latex\' final_name]; %for Microsoft
else
    addpath([base_path '/Latex'])
    path.out = ['Latex/' final_name]; %for Mac    
end


matrix2latex_black(matrix, path.out, 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.1f', 'size', 'small'); 
