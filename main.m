 
addpath('include/FEM/');
addpath('Mesh/');

plotSnapShot = true;

fileName = 'layercake.msh';

msh = constructMesh(fileName);

msh.rhs = find(msh.coords(:,1) < min(msh.coords(:,1)) + 1e-6);

msh.bnd = [msh.rhs; msh.rhs + msh.nnode; msh.rhs + 2 * msh.nnode];

numSamples = 10;

theta = pi * rand(numSamples, 5);

SS = @(x) makeSnapShot(msh,x);

for i = 1 : numSamples
    
    disp(strcat('Building Snapshot \t',int2str(i)));
    
    [U{i},K{i},F{i}] = SS(theta(i,:));
    
    if(plotSnapShot)
    
        vector_data.name = 'displacement';

        vector_data.data = zeros(msh.nnode,3);

        vector_data.data(:,1) = U{i}(1:msh.nnode);
        vector_data.data(:,2) = U{i}(msh.nnode +1 :2 *msh.nnode);
        vector_data.data(:,3) = U{i}(2 * msh.nnode + 1 : 3 *msh.nnode);

        name = strcat('snapShot_',int2str(i),'.vtk');

        matlab2vtk(name,'', msh, 'hex', [], vector_data, []);
        
    end
      
end


