function msh = constructMesh(fileName)


msh = readMesh(fileName,'HEXAS');

msh.rhs = find(msh.coords(:,1) < 1e-6);

msh.bnd = [msh.rhs, msh.rhs + msh.nnode, msh.rhs + 2 * msh.nnode];

msh.free = 1 : msh.tdof; msh.free(msh.bnd) = [];

msh.top = find(msh.coords(:,3) > max(msh.coords(:,3)) - 1e-6);


% Data will be x, y, z displacements at these nodes from DIC camera

msh.M = sparse(3 * length(msh.top), msh.tdof);

for i = 1 : length(msh.top)

    msh.M(i,msh.top(i)) = 1;
    msh.M(i + msh.nnode,msh.top(i) + msh.nnode) = 1;
    msh.M(i + 2 * msh.nnode,msh.top(i) + 2 * msh.nnode) = 1;

end



end
