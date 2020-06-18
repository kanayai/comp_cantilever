function [Ke, F] = elementStiffness(coords,elem,N,dN,nip,C,ang,intpoints,rho)

    Ke = zeros(24);
    
    F = zeros(24,1);

    
    g = 9.81;
    
    for ip = 1 : nip

        J	=	coords'*dN{ip};
        dNdX	=	dN{ip}*inv(J); %#ok<MINV>

        % For Composite Elements

        C = Rotation_Matrices(C,[0,0,ang]);

        % Local Coordinates to Global Coordinates

        th = zeros(3,1);

        C = Rotation_Matrices(C,th);
        
        % Calculate Strain Matrix
        
        B = zeros(6,24);
        
       

        B(1,1:8) = dNdX(:,1); % e_11 = u_1,1
        B(2,9:16) = dNdX(:,2); % e_22 = u_2,2 
        B(3,17:24) = dNdX(:,3); % e_33 = u_3,3

        B(4,9:16) = dNdX(:,3);	B(4,17:24) = dNdX(:,2);	% e_23 = u_2,3 + u_3,2
        B(5,1:8) = dNdX(:,3);	B(5,17:24) = dNdX(:,1);	% e_13 = u_1,3 + u_3,1
        B(6,1:8) = dNdX(:,2);	B(6,9:16) = dNdX(:,1);	% e_12 = u_1,2 + u_2,1
        
        Ke = Ke + B' * C * B * intpoints.wgts(ip) * det(J);
        
        F(17:24) = F(17:24) + N{ip} * intpoints.wgts(ip) * det(J);

    end	% for each ip

    Ke = 0.5*(Ke + Ke');
    
    F = -F * rho * g;
    
    if sum(sum(abs(Ke-Ke'))) > 1e-12
        error('Element Stiffness Matrix Not Symmetric');
    end
    
    Ke = Ke(:); % Convert to a column vector
    
    
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


% function Ke = elementStiffness(coords,elem,dN,nip,C,ang,intpoints)
% 
%     Ke = zeros(24);
% 
% 	for ip = 1 : nip
% 
% 	  	J	=	coords'*dN{ip};
% 	    dNdX	=	dN{ip}*inv(J); %#ok<MINV>
% 
% 	    % For Composite Elements
% 
% 		C = Rotation_Matrices(C,[0,0,ang]);
% 
% 		% Local Coordinates to Global Coordinates
% 
% 		th = zeros(3,1);
% 
% 		C = Rotation_Matrices(C,th);
%         
%         % Calculate Strain Matrix
%         
%         B = zeros(6,24);
% 
% 		B(1,1:8) = dNdX(:,1); % e_11 = u_1,1
% 		B(2,9:16) = dNdX(:,2); % e_22 = u_2,2 
% 		B(3,17:24) = dNdX(:,3); % e_33 = u_3,3
% 
% 		B(4,9:16) = dNdX(:,3);	B(4,17:24) = dNdX(:,2);	% e_23 = u_2,3 + u_3,2
% 		B(5,1:8) = dNdX(:,3);	B(5,17:24) = dNdX(:,1);	% e_13 = u_1,3 + u_3,1
% 		B(6,1:8) = dNdX(:,2);	B(6,9:16) = dNdX(:,1);	% e_12 = u_1,2 + u_2,1
%         
% 	    Ke = Ke + B'*C*B*intpoints.wgts(ip)*det(J);
% 
% 	end	% for each ip
% 
%     Ke = 0.5*(Ke + Ke');
%     
%     if sum(sum(abs(Ke-Ke'))) > 1e-12
%         error('Element Stiffness Matrix Not Symmetric');
%     end
%     
%     Ke = Ke(:); % Convert to a column vector
% 
% end	


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

function C = Rotation_Matrices(C,th)

R = cell(3,1);

% Rotation about x_1

c = cos(th(1)); s = sin(th(1));

R{1} = zeros(6);

R{1}(1,1) = 1;
R{1}(2,2) = c * c; 
R{1}(2,3) = s * s;	
R{1}(2,4) = 2 * c * s;
R{1}(3,2) = s * s;
R{1}(3,3) = c * c;
R{1}(3,4) = -2 * c * s;
R{1}(4,2) = -c * s;
R{1}(4,3) = c * s;
R{1}(4,4) = (c * c) - (s * s);
R{1}(5,5) = c;
R{1}(5,6) = -s;
R{1}(6,5) = s;
R{1}(6,6) = c;

% Rotation about x_2

c = cos(th(2)); s = sin(th(2));

R{2} = zeros(6);

R{2}(1,1) = c * c;
R{2}(1,3) = s * s;
R{2}(1,5) = 2 * c * s;
R{2}(2,2) = 1;
R{2}(3,1) = s * s;
R{2}(3,3) = c * c;
R{2}(3,5) = -2 * c * s;
R{2}(4,4) = c;
R{2}(4,6) = -s;
R{2}(5,1) = -c * s;
R{2}(5,3) = c * s;
R{2}(5,5) = (c * c) - (s * s);
R{2}(6,4) = s;
R{2}(6,6) = c;

% Rotation about x_3

c = cos(th(3)); s = sin(th(3));

R{3} = zeros(6);

R{3}(1,1) = c * c;
R{3}(1,2) = s * s;
R{3}(1,6) = 2 * c * s;
R{3}(2,1) = s * s;
R{3}(2,2) = c * c;
R{3}(2,6) = -2 * c * s;
R{3}(3,3) = 1;
R{3}(4,4) = c;
R{3}(4,5) = s;
R{3}(5,4) = -s;
R{3}(5,5) = c;
R{3}(6,1) = -c * s;
R{3}(6,2) = c * s;
R{3}(6,6) = (c * c) - (s * s);

T = R{3}*R{2}*R{1};

C = T*C*T';


end