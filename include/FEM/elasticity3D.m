function [K, F] = elasticity3D(theta, msh)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Model Setup 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numLayers = 9;
ss = [theta(1),-1,theta(2),-1,theta(3),-1,theta(4),-1,theta(5)];

E_R = 4.5; %    GPa
nu_R = 0.35;

E1 = 135; %    GPa
E2 = 8.5; %    GPa
E3 = 8.5; %    GPa

nu_21 = 0.022; 
nu_31 = 0.022;
nu_32 = 0.5;

G_12 = 5;   %   GPa
G_13 = 5;   %   GPa
G_23 = 5; %   GPa

rho = 0.0001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[isotropic,composite] = makeMaterials(E_R,nu_R,E1,E2,E3,nu_21,nu_31,nu_32,G_12,G_13,G_23);

msh = defineIPs(msh);

[msh.N,msh.dN] = ShapeFunctions(msh);

% For Each Element

indx_j = repmat(1:24,24,1); indx_i = indx_j';
Kindx.i = msh.e2g(:,indx_i(:)); Kindx.j = msh.e2g(:,indx_j(:));

Ke_all = zeros(24^2,msh.nelem); 

F = zeros(msh.tdof,1);

for ie = 1:msh.nelem
    
    if ss(msh.elements{ie}.region) >= 0
        C = composite;
        ang = ss(msh.elements{ie}.region);
    else
        C = isotropic;
        ang = 0;
    end
        
    [Ke_all(:,ie),F(msh.e2g(ie,:))] = elementStiffness(msh.coords(msh.elements{ie}.connectivity,:),msh.elements{ie},msh.N,msh.dN,msh.nip,C,ang,msh.ip,rho); 
    
end

K = sparse(Kindx.i',Kindx.j',Ke_all);


end

function U = Convert2Output(msh,V)

U = zeros(msh.nnode,3);

U(:,1) = V(1:msh.nnode);
U(:,2) = V(msh.nnode+1:2*msh.nnode);
U(:,3) = V(2*msh.nnode+1:end);

end

function [isotropic,composite] = makeMaterials(E_R,nu_R,E1,E2,E3,nu_21,nu_31,nu_32,G_12,G_13,G_23)

% Isotropic Resin
lambda = E_R*nu_R/((1+nu_R)*(1-2*nu_R)); mu = E_R/(2*(1+nu_R));
isotropic = zeros(6);
isotropic(1:3,1:3) = lambda*ones(3);
isotropic = isotropic + diag(mu*[2,2,2,1,1,1]);

% Orthotropic Composite
S = zeros(6);
S(1,1) = 1/E1; S(1,2) = -nu_21/E2; S(1,3) = -nu_31/E3;
S(2,1) = S(1,2);    S(2,2) = 1/E2; S(2,3) = -nu_32/E3;
S(3,1) = S(1,3);    S(3,2) = S(2,3);    S(3,3) = 1/E3;
S(4,4) = 1/G_23;
S(5,5) = 1/G_13;
S(6,6) = 1/G_12;

composite = inv(S);


end

function msh = defineIPs(msh)

    msh.nip = 8;
    
    P   =   0.577350269189626;
    
    msh.ip.coords(1,:) = [-P,-P,-P];
    msh.ip.coords(2,:) = [P,-P,-P];
    msh.ip.coords(3,:) = [P,-P,P];
    msh.ip.coords(4,:) = [P,-P,P];
    
    msh.ip.coords(5,:) = [-P,P,-P];
    msh.ip.coords(6,:) = [P,P,-P];
    msh.ip.coords(7,:) = [P,P,P];
    msh.ip.coords(8,:) = [P,P,P];

    msh.ip.wgts(1:8) = ones(8,1);

end


function [N,dN] = ShapeFunctions(msh)


	N = cell(msh.nip,1); dN = cell(msh.nip,1);

	for ip = 1:msh.nip

		N{ip} = zeros(8,1); dN{ip} = zeros(8,3);

		xi = msh.ip.coords(ip,1);	eta = msh.ip.coords(ip,2);	mu = msh.ip.coords(ip,3);

		% Displacement Shape Functions
		
		N{ip}(1) = 0.125 * (1 - xi) * (1 - eta) * (1 - mu);
		N{ip}(2) = 0.125 * (1 + xi) * (1 - eta) * (1 - mu);
		N{ip}(3) = 0.125 * (1 + xi) * (1 + eta) * (1 - mu);
		N{ip}(4) = 0.125 * (1 - xi) * (1 + eta) * (1 - mu);
		N{ip}(5) = 0.125 * (1 - xi) * (1 - eta) * (1 + mu);
		N{ip}(6) = 0.125 * (1 + xi) * (1 - eta) * (1 + mu);
		N{ip}(7) = 0.125 * (1 + xi) * (1 + eta) * (1 + mu);
		N{ip}(8) = 0.125 * (1 - xi) * (1 + eta) * (1 + mu);

		% Derivative Shape Functions

		dN{ip}(1,1) = -0.125 * (1 - eta) * (1 - mu);
		dN{ip}(2,1) =  0.125 * (1 - eta) * (1 - mu);
		dN{ip}(3,1) =  0.125 * (1 + eta) * (1 - mu);
		dN{ip}(4,1) = -0.125 * (1 + eta) * (1 - mu);
		dN{ip}(5,1) = -0.125 * (1 - eta) * (1 + mu);
		dN{ip}(6,1) =  0.125 * (1 - eta) * (1 + mu);
		dN{ip}(7,1) =  0.125 * (1 + eta) * (1 + mu);
		dN{ip}(8,1) = -0.125 * (1 + eta) * (1 + mu);


		dN{ip}(1,2) = -0.125 * (1 - xi) * (1 - mu);
		dN{ip}(2,2) = -0.125 * (1 + xi) * (1 - mu);
		dN{ip}(3,2) =  0.125 * (1 + xi) * (1 - mu);
		dN{ip}(4,2) =  0.125 * (1 - xi) * (1 - mu);
		dN{ip}(5,2) = -0.125 * (1 - xi) * (1 + mu);
		dN{ip}(6,2) = -0.125 * (1 + xi) * (1 + mu);
		dN{ip}(7,2) =  0.125 * (1 + xi) * (1 + mu);
		dN{ip}(8,2) =  0.125 * (1 - xi) * (1 + mu);


		dN{ip}(1,3) = -0.125 * (1 - xi) * (1 - eta);
		dN{ip}(2,3) = -0.125 * (1 + xi) * (1 - eta);
		dN{ip}(3,3) = -0.125 * (1 + xi) * (1 + eta);
		dN{ip}(4,3) = -0.125 * (1 - xi) * (1 + eta);
		dN{ip}(5,3) =  0.125 * (1 - xi) * (1 - eta);
		dN{ip}(6,3) =  0.125 * (1 + xi) * (1 - eta);
		dN{ip}(7,3) =  0.125 * (1 + xi) * (1 + eta);
		dN{ip}(8,3) =  0.125 * (1 - xi) * (1 + eta);

	end	%	end	N and dN

end	% 3DShapeFunctions