function matlab2vtk (filename,title, msh, elementType, scalar_point, vector_point, scalar_cell)

output_unit = fopen(filename,'w+');

switch lower(elementType)

case 'line'

element_order = 2;
cell_type = 3;

case 'quad'

element_order = 4;
cell_type = 9;

case 'hex'
element_order = 8;
cell_type = 12;

case 'tet'

element_order = 4; % Linear
cell_type = 10;

end

fprintf ( output_unit, '# vtk DataFile Version 2.0\n' );
fprintf ( output_unit, '%s\n', title );
fprintf ( output_unit, 'ASCII\n' );
fprintf ( output_unit, '\n' );
fprintf ( output_unit, 'DATASET UNSTRUCTURED_GRID\n' );
fprintf ( output_unit, 'POINTS %d double\n', msh.nnode );

for i = 1 : msh.nnode
fprintf ( output_unit, '  %f  %f  %f\n', msh.coords(i,:) );
end

fprintf ( output_unit, '\n' );
fprintf ( output_unit, 'CELLS  %d  %d\n', msh.nelem, (element_order+1)*msh.nelem );
for ie = 1 : msh.nelem
fprintf ( output_unit, '  %d', element_order );
for j = 1 : element_order
fprintf ( output_unit, '  %d', msh.elements{ie}.connectivity(j) - 1 ); % '-1' due to 0 node numbering
end
fprintf ( output_unit, '\n' );
end

fprintf ( output_unit, '\n' );
fprintf ( output_unit, 'CELL_TYPES %d\n', msh.nelem );


for i = 1 : msh.nelem
fprintf ( output_unit, '%d\n', cell_type);
end


if (isempty(scalar_point)==0) || (isempty(vector_point)==0)
% POINT_DATA
fprintf ( output_unit, '\n' );
fprintf ( output_unit, 'POINT_DATA %d\n', msh.nnode );

if isempty(scalar_point)==0
        
% SCALAR
fprintf ( output_unit, 'SCALARS %s double\n', scalar_point.name );
fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
for i = 1 : msh.nnode
fprintf ( output_unit, '  %f\n', scalar_point.data(i) );
end
end

      if isempty(vector_point)==0
          % VECTOR
          fprintf ( output_unit, 'VECTORS %s double\n', vector_point.name);
          for i = 1 : msh.nnode
            fprintf ( output_unit, '  %f  %f  %f\n',  vector_point.data(i,:));
          end
    end
end

% CELL_DATA
if isempty(scalar_cell) == 0
    
    for i = 1 : length(scalar_cell)
     
        
        if(i == 1)    
            fprintf ( output_unit, '\n');
            fprintf ( output_unit, 'CELL_DATA %i\n', msh.nelem );
            fprintf ( output_unit, 'SCALARS %s double\n', scalar_cell{i}.name);
            fprintf ( output_unit, 'LOOKUP_TABLE default\n');
            for j = 1 : msh.nelem
                fprintf ( output_unit, ' %d\n', scalar_cell{i}.data(j) );
            end
        
        else
            
            if(i == 2)
                fprintf ( output_unit, '\n');
                fprintf ( output_unit, 'FIELD FieldData %i\n', length(scalar_cell) - 1 );
            end
                fprintf ( output_unit, '%s ', scalar_cell{i}.name);
                fprintf ( output_unit, '%i %i double\n', [1, msh.nelem]);
               % fprintf ( output_unit, 'LOOKUP_TABLE default\n');
                for j = 1 : msh.nelem
                    fprintf ( output_unit, ' %d\n', scalar_cell{i}.data(j) );
                end
        end
     end
            
        
end



return
end

