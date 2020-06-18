


function [Ui,Ki,Fi] = makeSnapShot(msh,theta)

for i=1:size(theta,1)

    [Ki,Fi] = elasticity3D(theta, msh);

    Ui = zeros(msh.tdof,1);

    Ui(msh.free) = Ki(msh.free,msh.free) \ (Fi(msh.free) - Ki(msh.free,msh.bnd) * Ui(msh.bnd));
 
    if (mod(i,10)==0)
        fprintf('snapshot %d done\n', i);
    end
    
end    

end