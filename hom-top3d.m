clear; close all; clc; workspace; tic;
%% SECTION 1 - USER-DEFINED DESIGN PARAMETERS
scale_law = [1598.0      657.68      29.65;...     % Scale-law constants for Gyroid TPMS structure                                
             661.93      406.03      15.69;...
             515.83      167.02      4.478];
base_material  = [2100; 0.42];                      % Young's Modulus in MPa, and Poisson's Ratio                                                        
% Preprocess parameters
cell_size = 10;                                     % mm
cell_counts = [5, 1, 1];                            % Number of unit cell is x, y, and z directions, respectively                                                           
plate_thickness = 1;                                % Upper and lower plate thickness in mm                                                                     
total_load = 2000;                                  % Total load in N
% Solution parameters
minumum_volfrac  = 0.1;                             % Minimum volume fraction constraint
maximum_volfrac  = 0.5;                             % Maximum volume fraction constraint
target_volfrac = 0.3;                               % Total volume fraction constraint
r_min = 1.5;
max_loop = 20;                                      % Maximum number of iterations
tol_x = 0.001;                                      % Terminarion criterion
displayflag_loop = 1;                               % Display density distribution at each iteration (this will increase the termination time)
% Lattice morphology 
f = @(x,y,z,k) (cos(k(1)*x).*sin(k(2)*y) + cos(k(2)*y).*sin(k(3)*z)...    % Gyroid TPMS function
    + cos(k(3)*z).*sin(k(1)*x));
ftf= @(p) 1.52*p;                                   % t to rho* mapping
% Reconcsturction parameters
mesh_scaling = 5;                                   % Defines mesh density of reconstructed geometry
fem_reconstruction = 1;                             % Switch for fem reconstruction on: 1 and off: 0
stl_reconstruction = 1;                             % Switch for stl reconstruction on: 1 and off: 0
%% SECTION 2 - START OF FEA AND OPTIMIZATION
nelx = cell_counts(1)*cell_size;
nely = cell_counts(2)*cell_size + (2*plate_thickness);
nelz = cell_counts(3)*cell_size;
VF = ones(nely, nelx, nelz);                                                                  % create a VF array
VF(plate_thickness + 1: end - plate_thickness, :, :) = (maximum_volfrac + minumum_volfrac)/2; % assign initial volume fraction to design domain
type = zeros(nely, nelx, nelz);                                                               % create a type array
type(plate_thickness + 1: end - plate_thickness, :, :) = 1;                                   % assign type to design domain 0 means non-design domain, 1 means design domain
% User-defined load DOFs
[il,jl,kl] = meshgrid(nelx,nely,0:nelz);					% Coordinates			
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); 	% Node IDs
loaddof = 3*loadnid(:) - 1;    								% DOFs			                         
% User-defined support DOFs
[iifL,jfL,kfL] = meshgrid(0,0,0:nelz);       					% Coordinates Roller BC                 
fixednidL = kfL*(nelx+1)*(nely+1)+iifL*(nely+1)+(nely+1-jfL);	% Node IDs 	   
[iifR,jfR,kfR] = meshgrid(nelx, 0:nely, 0:nelz);         		% Coordinates Symmetry Plane	          
fixednidR = kfR*(nelx+1)*(nely+1)+iifR*(nely+1)+(nely+1-jfR); 	% Node IDs   
fixeddof = [3*fixednidL(:)-1; 3*fixednidR(:)-2]; 				% DOFs		
% Prepare for finite element analysis
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1, -total_load/(nelz+1) ,ndof,1);            
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
% KE = lk_H8(nu);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
% Prepare filter
iH = ones(nele*(2*(ceil(r_min)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
            for k2 = max(k1-(ceil(r_min)-1),1):min(k1+(ceil(r_min)-1),nelz)
                for i2 = max(i1-(ceil(r_min)-1),1):min(i1+(ceil(r_min)-1),nelx)
                    for j2 = max(j1-(ceil(r_min)-1),1):min(j1+(ceil(r_min)-1),nely)
                        e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,r_min-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
x = VF;
passive = x;
xPhys = x;
loop = 0; 
change = 1;
ce=zeros(nele,1);
ces=zeros(nele,1);
while change > tol_x && loop < max_loop
    loop = loop+1;
    % Finite element analysis
    for i=1:nele
        KE=hex(xPhys(i),1,1,1,type(i),base_material,scale_law);
        sK((((i-1)*24*24)+1):i*24*24,1)=reshape(KE,24*24,1);
    end  
    K = sparse(iK,jK,sK); K = (K+K')/2;
    tolit = 1e-8;
    maxit = 8000;
    M = diag(diag(K(freedofs,freedofs)));
    U(freedofs,:) = pcg(K(freedofs,freedofs),F(freedofs,:),tolit,maxit,M);
    % Objective function and sensitivity analysis
    for m=1:nele
        ce(m,1)= U(edofMat(m,:))'*hex(xPhys(m),1,1,1,type(m),base_material,...
            scale_law)*U(edofMat(m,:));
    end
    c = sum(ce)/2;
    for n=1:nele
         ces(n,1)= U(edofMat(n,:))'*(hex_sens(xPhys(n),1,1,1,type(n),...
             base_material,scale_law))*U(edofMat(n,:));
    end
    dc=-reshape(ces,[nely,nelx,nelz]);
    dv = ones(nely,nelx,nelz);
    % Filtering and modification of sensitivities
    dc(:) = H*(dc(:)./Hs);  
    dv(:) = H*(dv(:)./Hs);
    % Optimality criteria update
    l1 = 0; l2 = 1e9; move = 0.2;
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        xnew = max(minumum_volfrac,max(x-move,min(maximum_volfrac,min(x+move,x.*sqrt(-dc./dv/lmid)))));
        xPhys(:) = (H*xnew(:))./Hs;
        xPhys(passive==1) = 1;
        if (sum(xPhys(:))-(2*plate_thickness)*nelx*nelz) > target_volfrac*(nelx*(nely-2*plate_thickness)*nelz), l1 = lmid; else l2 = lmid; end
    end
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,(sum(xPhys(:))-(2*plate_thickness)*nelx*nelz)/(nelx*(nely-2*plate_thickness)*nelz),change);
    if displayflag_loop, clf; display_3D(xPhys,minumum_volfrac,maximum_volfrac); end %#ok<UNRCH>
end
clf; display_3D(xPhys,minumum_volfrac,maximum_volfrac); % Display density results
toc
%% SECTION 3 - CALL FOR RECONSTRUCTION
if fem_reconstruction || stl_reconstruction
elem_connectivity_contour = element_setup(edofMat,xPhys,nele);
averaged_densities = averaging_densities(elem_connectivity_contour);
reconstruct("test", averaged_densities, cell_size, cell_counts, mesh_scaling, f, ftf, fem_reconstruction,stl_reconstruction);
end
%% SECTION 4 - DISPLAY 3D TOPOLOGY (ISO-VIEW)
function display_3D(rho,minimum_volfrac,maximum_volfrac)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;            % User-defined unit element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) >= minimum_volfrac - 1E3  && rho(j,i,k) <= maximum_volfrac + 1E3)  % User-defined display to exclude plates
                vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                patch('Faces', face, 'Vertices', vert, 'FaceColor', max(min([0.2 + 0.8 * (1 - (rho(j, i, k) - minimum_volfrac) / (maximum_volfrac - minimum_volfrac)), ...
                                                              0.2 + 0.8 * (1 - (rho(j, i, k) - minimum_volfrac) / (maximum_volfrac - minimum_volfrac)), ...
                                                              0.2 + 0.8 * (1 - (rho(j, i, k) - minimum_volfrac) / (maximum_volfrac - minimum_volfrac))], 1), 0));
                hold on;
            end
        end
    end
end
axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end
%% SECTION 5 - GENERATE ELEMENT STIFFNESS MATRIX 
function Kelem=hex(xPhys,length_x,length_y,length_z,type,base_material,scale_law)
if type == 1   % Lattice Material Model
    C11=polyval(scale_law(1, :), xPhys);
    C12=polyval(scale_law(2, :), xPhys);
    C44=polyval(scale_law(3, :), xPhys);
else            % Full Material Model
    C12 = (base_material(1)*base_material(2))/((1+base_material(2))*(1-2*base_material(2)));
    C44 = base_material(1)/(2*(1+base_material(2)));
    C11 = C12 + 2*C44;
end
C =[C11 C12 C12 0   0   0;
    C12 C11 C12 0   0   0;
    C12 C12 C11 0   0   0; 
    0   0   0   C44 0   0; 
    0   0   0   0   C44 0;
    0   0   0   0   0   C44];
GaussPoint = [-1/sqrt(3), 1/sqrt(3)];
coordinates = zeros(8,3);
coordinates(1,:) = [-length_x/2 -length_y/2 -length_z/2];
coordinates(2,:) = [length_x/2 -length_y/2 -length_z/2];
coordinates(3,:) = [length_x/2 length_y/2 -length_z/2];
coordinates(4,:) = [-length_x/2 length_y/2 -length_z/2];
coordinates(5,:) = [-length_x/2 -length_y/2 length_z/2];
coordinates(6,:) = [length_x/2 -length_y/2 length_z/2];
coordinates(7,:) = [length_x/2 length_y/2 length_z/2];
coordinates(8,:) = [-length_x/2 length_y/2 length_z/2];
Kelem = zeros (24,24);
for xi1=GaussPoint
    for xi2=GaussPoint
        for xi3=GaussPoint
            dShape = (1/8)*[-(1-xi2)*(1-xi3),(1-xi2)*(1-xi3),...
                (1+xi2)*(1-xi3),-(1+xi2)*(1-xi3),-(1-xi2)*(1+xi3),...
                (1-xi2)*(1+xi3),(1+xi2)*(1+xi3),-(1+xi2)*(1+xi3);...
                -(1-xi1)*(1-xi3),-(1+xi1)*(1-xi3),(1+xi1)*(1-xi3),...
                (1-xi1)*(1-xi3),-(1-xi1)*(1+xi3),-(1+xi1)*(1+xi3),...
                (1+xi1)*(1+xi3),(1-xi1)*(1+xi3);-(1-xi1)*(1-xi2),...
                -(1+xi1)*(1-xi2),-(1+xi1)*(1+xi2),-(1-xi1)*(1+xi2),...
                (1-xi1)*(1-xi2),(1+xi1)*(1-xi2),(1+xi1)*(1+xi2),...
                (1-xi1)*(1+xi2)];
            JacobianMatrix = dShape*coordinates;
            auxiliar = inv(JacobianMatrix)*dShape;
            B = zeros(6,24);
            for i=1:3
                for j=0:7
                    B(i,3*j+1+(i-1)) = auxiliar(i,j+1);
                end
            end
            for j=0:7
                B(4,3*j+1) = auxiliar(2,j+1);
            end
            for j=0:7
                B(4,3*j+2) = auxiliar(1,j+1);
            end
            for j=0:7
                B(5,3*j+3) = auxiliar(2,j+1);
            end
            for j=0:7
                B(5,3*j+2) = auxiliar(3,j+1);
            end
            for j=0:7
                B(6,3*j+1) = auxiliar(3,j+1);
            end
            for j=0:7
                B(6,3*j+3) = auxiliar(1,j+1);
            end     
            Kelem = Kelem + B'*C*B*det(JacobianMatrix);
        end
    end
end
end
%% SECTION 6 - GENERATE ELEMENT SENSITIVITY MATRIX
function Ksens=hex_sens(xPhys,length_x,length_y,length_z,type,base_material,scale_law)
if type == 1   % Lattice Material Model Sensitivity 
    C11 = polyval(polyder(scale_law(1,:)), xPhys);
    C12 = polyval(polyder(scale_law(2,:)), xPhys);
    C44 = polyval(polyder(scale_law(3,:)), xPhys);
else            % Full Material Model Sensitivity
    [C11, C12, C44] = deal(0, 0, 0);
end
C =[C11 C12 C12 0   0   0; 
    C12 C11 C12 0   0   0;
    C12 C12 C11 0   0   0; 
    0   0   0   C44 0   0;
    0   0   0   0   C44 0;
    0   0   0   0   0   C44];
GaussPoint = [-1/sqrt(3), 1/sqrt(3)];
coordinates = zeros(8,3);
coordinates(1,:) = [-length_x/2 -length_y/2 -length_z/2];
coordinates(2,:) = [length_x/2 -length_y/2 -length_z/2];
coordinates(3,:) = [length_x/2 length_y/2 -length_z/2];
coordinates(4,:) = [-length_x/2 length_y/2 -length_z/2];
coordinates(5,:) = [-length_x/2 -length_y/2 length_z/2];
coordinates(6,:) = [length_x/2 -length_y/2 length_z/2];
coordinates(7,:) = [length_x/2 length_y/2 length_z/2];
coordinates(8,:) = [-length_x/2 length_y/2 length_z/2];
Ksens = zeros (24,24);
for xi1=GaussPoint
    for xi2=GaussPoint
        for xi3=GaussPoint
            dShape = (1/8)*[-(1-xi2)*(1-xi3),(1-xi2)*(1-xi3),...
                (1+xi2)*(1-xi3),-(1+xi2)*(1-xi3),-(1-xi2)*(1+xi3),...
                (1-xi2)*(1+xi3),(1+xi2)*(1+xi3),-(1+xi2)*(1+xi3);...
                -(1-xi1)*(1-xi3),-(1+xi1)*(1-xi3),(1+xi1)*(1-xi3),...
                (1-xi1)*(1-xi3),-(1-xi1)*(1+xi3),-(1+xi1)*(1+xi3),...
                (1+xi1)*(1+xi3),(1-xi1)*(1+xi3);-(1-xi1)*(1-xi2),...
                -(1+xi1)*(1-xi2),-(1+xi1)*(1+xi2),-(1-xi1)*(1+xi2),...
                (1-xi1)*(1-xi2),(1+xi1)*(1-xi2),(1+xi1)*(1+xi2),...
                (1-xi1)*(1+xi2)];
            JacobianMatrix = dShape*coordinates;
            auxiliar = inv(JacobianMatrix)*dShape;
            B = zeros(6,24);
            for i=1:3
                for j=0:7
                    B(i,3*j+1+(i-1)) = auxiliar(i,j+1);
                end
            end
            for j=0:7
                B(4,3*j+1) = auxiliar(2,j+1);
            end
            for j=0:7
                B(4,3*j+2) = auxiliar(1,j+1);
            end
            for j=0:7
                B(5,3*j+3) = auxiliar(2,j+1);
            end
            for j=0:7
                B(5,3*j+2) = auxiliar(3,j+1);
            end
            for j=0:7
                B(6,3*j+1) = auxiliar(3,j+1);
            end
            for j=0:7
                B(6,3*j+3) = auxiliar(1,j+1);
            end 
            Ksens = Ksens + B'*C*B*det(JacobianMatrix);
        end
    end
end
end
%% SECTION 7 - DATA PREPARATION FOR RECONSTRUCTION
function elem_connectivity_contour = element_setup(edofMat,xPhys,nele)
xPhysClmn = reshape(xPhys, [], 1);
x_edofMat = zeros(nele, 8);
x_edofMat(:, 1:8) = edofMat(:, 3:3:24) / 3;
elem_connectivity_contour = [(1:nele)', x_edofMat, xPhysClmn];
elem_connectivity_contour(elem_connectivity_contour(:, 10) == 1, :) = [];
end
%% SECTION 8 - RELATIVE DENSITY AVERAGING
function averaged_densities = averaging_densities(elem_connectivity_contour)
data = elem_connectivity_contour(:, 2:end-1);
nodeNumbers = unique(data(:));
numNodes = numel(nodeNumbers);
averaged_densities = zeros(numNodes, 2);
for idx = 1:numNodes
    currentNode = nodeNumbers(idx); 
    rows = any(data == currentNode, 2);
    values = elem_connectivity_contour(rows, 10);
    if ~isempty(values)
        averaged_densities(idx, 1) = currentNode;
        averaged_densities(idx, 2) = mean(values);
    end
end
averaged_densities(all(averaged_densities == 0, 2), :) = [];
end
%% SECTION 9 - STL AND FEM RECONSTRUCTION
function reconstruct(file_name, averaged_densities, cell_size, cell_counts, mesh_scaling, f, ftf, fem_reconstruction, stl_reconstruction) 
le = cell_counts * cell_size;
node_inc = 1;
[X1,Y1,Z1] = ndgrid(0:node_inc:le(1), 0:node_inc:le(2),0:node_inc:le(3));
fine_size = mesh_scaling.*le; 
rhomat = reshape(averaged_densities(:,2), le(2) + 1, le(1) + 1, le(3) + 1);
rhomat = permute(rhomat, [2 1 3]);
F=griddedInterpolant(X1,Y1,Z1,rhomat);
[X_finer,Y_finer,Z_finer] = meshgrid(linspace(0,le(1),fine_size(1)),...
                                     linspace(0,le(2),fine_size(2)),...
                                     linspace(0,le(3),fine_size(3)));
rho_finer=F(X_finer,Y_finer,Z_finer);
k = 2*pi*cell_counts./le; % Coefficents for the TPMS function
[x,y,z] = meshgrid(linspace(0,le(1),fine_size(1)),linspace(0,le(2),fine_size(2)),linspace(0,le(3),fine_size(3)));
t = ftf(rho_finer); 
as = f(x,y,z,k);
as_stl = (as - t).*(as + t);
as_stl = flip(as_stl, 1);
if stl_reconstruction % STL write
[F,V]=isosurface(x, y, z, as_stl, 0);
[FC,VC]=isocaps(x, y, z, as_stl,0,'below');
F= [F; FC+length(V(:,1))];
V= [V; VC];
clear FC VC;
stlwrite(triangulation(F,V), strcat("./", file_name, ".stl"), "text");
end
if fem_reconstruction % FEM write
as = (as - t).*(as + t);
[x,y,z] = meshgrid(linspace(0,le(1),fine_size(1) + 1),linspace(0,le(2),fine_size(2) + 1),...
          linspace(0,le(3),fine_size(3) + 1));   
temp_nodes = [(1:(fine_size(1)+1)*(fine_size(2)+1)*(fine_size(3)+1))', reshape(x,[],1), reshape(y,[],1), reshape(z,[],1)];
max_elem_id = numel(find(as<0));
elements = zeros(max_elem_id, 9);
total = numel(find(as<0));
as = rot90(as,3);
K = find(as<0);
for iter=1:total
    i = K(iter);
    X = mod(i-1, fine_size(1)) + 1;
    Y = mod(floor((i-1)/fine_size(1)), fine_size(2)) + 1;
    Z = mod(floor((i-1)/fine_size(2)/fine_size(1)), fine_size(3)) + 1;
    index = Y + (X-1)*fine_size(2) + (Z-1)*fine_size(2)*fine_size(1) + X - 1 + (Z-1)*(fine_size(2)+fine_size(1)) +Z - 1;
    connectivity = [index, index + 1, index + fine_size(2) + 2,...
                   index + fine_size(2) + 1, index + fine_size(1) + fine_size(2) + fine_size(1)*fine_size(2) + 1,...
                   index + fine_size(1) + fine_size(2) + fine_size(1)*fine_size(2) + 2, index + fine_size(1) + 2*fine_size(2) + fine_size(1)*fine_size(2) + 3,...
                   index + fine_size(1) + 2*fine_size(2) + fine_size(1)*fine_size(2) + 2];
    elements(iter,:) = [0, connectivity];
end
check = elements(:,2:end);
node_bool = ismember((1:size(temp_nodes,1))', check);
elements = [(1:max_elem_id)',elements];
fileID = fopen(strcat(file_name, '.fem'), 'w');
fprintf(fileID,'%s\n','BEGIN BULK');
fprintf(fileID,'GRID  %10d        %-8.4f%-8.4f%-8.4f\n', temp_nodes(find(node_bool)',:)');
fprintf(fileID,'CHEXA %10d%8d%8d%8d%8d%8d%8d%8d\n+     %10d%8d\n', elements');
fprintf(fileID,'%s\n','ENDDATA');
fclose(fileID);
end
end
%% Scale Law Arrays (Change lines 3-5)
% scale_law = [2758.20     156.45      54.48;...     % Scale-law constants for Primitive TPMS structure                             
%              1116.90     277.47      34.01;...
%              457.73      259.96      6.43];

% scale_law = [1598.0      657.68      29.65;...     % Scale-law constants for Gyroid TPMS structure                                
%              661.93      406.03      15.69;...
%              515.83      167.02      4.478];

% scale_law = [1556.30     738.08      37.14;...     % Scale-law constants for Diamond TPMS structure                             
%              840.07      302.84      21.61;...
%              604.60      118.89      4.54];
%% CASE 1 -- Messerschmitt–Bölkow–Blohm beam (MBB) (Default)
    % % User-defined load DOFs
    % [il,jl,kl] = meshgrid(nelx,nely,0:nelz);					            % Coordinates			
    % loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); 	            % Node IDs
    % loaddof = 3*loadnid(:) - 1;    								        % DOFs			                         
    % % User-defined support DOFs
    % [iifL,jfL,kfL] = meshgrid(19,0,0:nelz);       					    % Coordinates Roller BC                 
    % fixednidL = kfL*(nelx+1)*(nely+1)+iifL*(nely+1)+(nely+1-jfL);	        % Node IDs 	   
    % [iifR,jfR,kfR] = meshgrid(nelx, 0:nely, 0:nelz);         		        % Coordinates Symmetry Plane	          
    % fixednidR = kfR*(nelx+1)*(nely+1)+iifR*(nely+1)+(nely+1-jfR); 	    % Node IDs   
    % fixeddof = [3*fixednidL(:)-1; 3*fixednidR(:)-2]; 				        % DOFs		
%% CASE 2 -- Cantilever beam    (Change lines 36-45 from MBB Beam Case which is default case)
    % % User-defined load DOFs
    % [il,jl,kl] = meshgrid(119, nely, 0:nelz);                        		% Coordinates
    % loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); 				% Node IDs
    % loaddof = 3*loadnid(:) - 1;                             				% DOFs
    % % User-defined support DOFs
    % [iifL,jfL,kfL] = meshgrid(0:cell_size, 0, 0:nelz);                    % Coordinates Bottom Face
    % fixednidL = kfL*(nelx+1)*(nely+1)+iifL*(nely+1)+(nely+1-jfL);   		% Node IDs
    % [iifR,jfR,kfR] = meshgrid(0:cell_size, nely, 0:nelz);                 % Coordinates Top Face
    % fixednidR = kfR*(nelx+1)*(nely+1)+iifR*(nely+1)+(nely+1-jfR);   		% Node IDs
    % fixeddof = [3*fixednidR(:); 3*fixednidR(:)-1; 3*fixednidR(:)-2;...	% DOFs
    %             3*fixednidL(:); 3*fixednidL(:)-1; 3*fixednidL(:)-2]; 
%% CASE 3 -- Flatwise compression (Change lines 36-45 from MBB Beam Case which is default case, and line 49 for distribuded load)
    % % User-defined load DOFs
    % [il,jl,kl] = meshgrid(0:nelx, nely, 0:nelz);						    % Coordinates
    % loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl);				% Node IDs
    % loaddof = 3*loadnid(:) - 1;                             			    % DOFs
    % % User-defined support DOFs
    % [iifR,jfR,kfR] = meshgrid(0:nelx, 0, 0:nelz);                  		% Coordinates
    % fixednidR = kfR*(nelx+1)*(nely+1)+iifR*(nely+1)+(nely+1-jfR);  		% Node IDs
    % fixeddof = [3*fixednidR(:); 3*fixednidR(:)-1; 3*fixednidR(:)-2]; 	    % DOFs

    % F = sparse(loaddof,1,-total_load/((nelx+1)*(nelz+1)),ndof,1);	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This MATLAB code was written by:                                         %
% M. Ozdemir - Karlsruhe Institute of Technology                           %
% U. Simsek - University of Michigan                                       %
% C. Gayir - Koc University                                                %
% Please sent your comments to:                                            % 
% ugursims@umich.edu,                                                      %
% ugur.simsek.16339@ozu.edu.tr,                                            %
% cgayir21@ku.edu.tr                                                       %
% mirhan.oezdemir@kit.edu                                                  %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% along with the usage procedures are discussed in the paper:              %
% "A Matlab code: anisotropic homogenization-based topology optimization   %
% for generating functionally graded surface-based lattices"               %
% M. Ozdemir, U. Simsek, C. Gayir, K. Gunaydin, O. Gulcan  [PAPER PENDING] %
%                                                                          %
% This code is available at https://github.com/cemal-efe-g/homo-top3d with % 
% any post-script additions or fixes.                                      %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserve all rights to this software. The software is         %
% provided  "as is", without warranty of any kind and the authors          %
% will not be liable for any claim that arises from its use.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
