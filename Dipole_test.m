%First 2D drawing of the Dipole coil surface

clc, clear, close all;

%% Input Parametrs for the GaTo MAIN Line

%Geometric parameters
dL=5.0E-3;       %[m]

Ncoils = 1;

aper = 0.07;        %Magnet aperture (diameter) [m]
R0 = 6.6/4;         %Reference curvature radius [m]
R_ext = R0+aper/2;  %External profile radius [m]
R_int = R0-aper/2;  %Internal profile radius [m]

R_h = aper/2;       %Magnet's head radius [m]

t=10.0E-2;          %Coil Thickness [m]
h=2.0E-2;          %Coil Height    [m]

defl = 45*pi/180;  %Angle of beam deflection [rad]

z_h = aper/2;      %Height of the heads [m]

th = 30*pi/180;      %Angle of the cos-theta section [rad]

CS = th/2*((aper/2+t)^2-(aper/2)^2);    % Conductor crosse-section, assuming it is a perfect circle sector [m**2]

%% Setting up the geometry

% External profile - section 1
Arc=R_ext*defl;
NR(1)=ceil(Arc/dL);
theta = linspace(0,defl,NR(1));
rho = linspace(R_ext,R_ext,NR(1));
[x1,y1] = pol2cart(theta,rho);
x1 = x1 - x1(1) + aper/2;
% y1 = y1 - y1(1) -;
z1 = zeros(size(x1));

% top magnet head - section 2
Arc=R_h*pi;
NR(2)=ceil(Arc/dL);
if mod(NR(2),2)>0; NR(2)=NR(2)+1;end
theta = linspace(defl,defl+pi,NR(2));
rho = linspace(R_h,R_h,NR(2));
[x2,y2] = pol2cart(theta,rho);
x2 = x2 - x2(1) + x1(end);
y2 = y2 - y2(1) + y1(end);
z2 = z_h*sin(linspace(0,pi,NR(2)));

% Internal profile - section 3
Arc=R_int*defl;
NR(3)=ceil(Arc/dL);
theta = linspace(defl,0,NR(3));
rho = linspace(R_int,R_int,NR(3));
[x3,y3] = pol2cart(theta,rho);
x3 = x3 - x3(1) + x2(end);
y3 = y3 - y3(1) + y2(end);
z3 = zeros(size(x3));

% bottom magnet head - section 4
Arc=R_h*pi;
NR(4)=ceil(Arc/dL);
if mod(NR(4),2)>0; NR(4)=NR(4)+1;end
theta = linspace(-pi,0,NR(4));
rho = linspace(R_h,R_h,NR(4));
[x4,y4] = pol2cart(theta,rho);
x4 = x4 - x4(1) + x3(end);
y4 = y4 - y4(1) + y3(end);
z4 = z_h*sin(linspace(0,pi,NR(4)));

%combine all sections
x = {x1, x2, x3, x4};
y = {y1, y2, y3, y4};
z = {z1, z2, z3, z4};
Nsections = length(x);

% combine coordinates into 3 vectors
X1 = cell(1,Nsections);
for i=1:Nsections
    X1{i}=[x{i};y{i};z{i}];
end


% Definition of rotation matrix
M_rot = @(v,rot)[cos(rot)+v(1)*v(1)*(1-cos(rot)), ...
    v(1)*v(2)*(1-cos(rot))-v(3)*sin(rot), ...
    v(1)*v(3)*(1-cos(rot))+v(2)*sin(rot); ...
    v(2)*v(1)*(1-cos(rot))+v(3)*sin(rot), ...
    cos(rot)+v(2)*v(2)*(1-cos(rot)), ...
    v(2)*v(3)*(1-cos(rot))-v(1)*sin(rot); ...
    v(3)*v(1)*(1-cos(rot))-v(2)*sin(rot), ...
    v(3)*v(2)*(1-cos(rot))+v(1)*sin(rot), ...
    cos(rot)+v(3)*v(3)*(1-cos(rot))];

%% plot for checking
figure; axes; hold on;
C = lines(Nsections);
for i=1:Nsections
    X1{i} = M_rot([0,0,1],45*pi/180)*X1{i};
    plot3(X1{i}(1,:),X1{i}(2,:),X1{i}(3,:),'s-','Color',C(i,:));
    X1{i} = M_rot([0,0,1],-45*pi/180)*X1{i};
end
axis equal;
axis tight;


%% Vector math for finding other edges
% find direction vector
Dt = cell(1,Nsections);
for i=1:Nsections
    Dt{i} = diff(X1{i},[],2);
end

L_wire = 0;
for i=1:Nsections
      L_wire = L_wire + sum(sqrt(sum(Dt{i}.^2,1)));
end


% take average direction vector
D = cell(1,Nsections);
for i=1:Nsections
    if i==1; idxp = Nsections; else idxp = i-1; end
    if i==Nsections; idxn = 1; else idxn = i+1; end

    Dtemp = [Dt{idxp}(:,end),Dt{i},Dt{idxn}(:,1)];
    D{i} = (Dtemp(:,1:end-1)+Dtemp(:,2:end))/2;
end

% normalize
for i=1:Nsections
   D{i} = bsxfun(@rdivide,D{i},sqrt(sum(D{i}.^2,1)));
end

% vertical vector
R = cell(1,Nsections);
for i=1:Nsections
    R{i} = zeros(size(D{i}));
    R{i}(3,:) = 1;
end



dTh = cell(1,1);
for i=1:Nsections
    dTh{i} = asin(Dt{i}(3,:)./sqrt(sum(Dt{i}.^2,1)));
    dTh{i} = [0,dTh{i}];
end


for i=2
    for j=1:size(X1{i},2)
        R{i}(:,j) = M_rot([cos(defl),sin(defl),0],(dTh{i}(2)))*R{i}(:,j);
    end
end

for i=4
    for j=1:size(X1{i},2)
        R{i}(:,j) = M_rot([1,0,0],-(dTh{i}(2)))*R{i}(:,j);
    end
end

% normal vector
N = cell(1,Nsections);
for i=1:Nsections
    N{i} = cross(D{i},R{i});
end

%% find other cable edges
X2 = cell(1,Nsections);
X3 = cell(1,Nsections);
X4 = cell(1,Nsections);
TH = cell(1,Nsections);

for i=1:Nsections
    X2{i} = X1{i} + h*R{i};
    X3{i} = X1{i} + h*R{i} + t*N{i};
    X4{i} = X1{i} + t*N{i};
end

%% vectors' check
figure; hold on;
for i=1:Nsections
    quiver3(X1{i}(1,:),X1{i}(2,:),X1{i}(3,:),N{i}(1,:),N{i}(2,:),N{i}(3,:),'k');
    quiver3(X1{i}(1,:),X1{i}(2,:),X1{i}(3,:),R{i}(1,:),R{i}(2,:),R{i}(3,:),'r');

%     quiver3(X1{i}(1,:),X1{i}(2,:),X1{i}(3,:),R2{i}(1,:),R2{i}(2,:),R2{i}(3,:),'k');
%     quiver3(X1{i}(1,:),X1{i}(2,:),X1{i}(3,:),D{i}(1,:),D{i}(2,:),D{i}(3,:),'b');
end
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
axis equal;


%% to create non-square cross-section
for i=1:4
    X2{i}(1,:) = X1{i}(1,:) - aper/2*(1-cos(th));
    X2{i}(3,:) = X1{i}(3,:) + aper/2*sin(th);

    X3{i}(1,:) = X4{i}(1,:) - (aper/2+t)*(1-cos(th));
    X3{i}(3,:) = X4{i}(3,:) + (aper/2+t)*sin(th);
end

% combine into surface matrix per section
xmat = cell(Nsections,1); ymat = cell(Nsections,1); zmat = cell(Nsections,1);
for i=1:Nsections
    % turn into surface matrices
    xmat{i} = [X1{i}(1,:)',X2{i}(1,:)',X3{i}(1,:)',X4{i}(1,:)'];
    ymat{i} = [X1{i}(2,:)',X2{i}(2,:)',X3{i}(2,:)',X4{i}(2,:)'];
    zmat{i} = [X1{i}(3,:)',X2{i}(3,:)',X3{i}(3,:)',X4{i}(3,:)'];
end

%% Plot edges
% Make figure
f = figure('Color','w','Position',[100,100,700,400]);
ax = axes('Parent',f);
hold(ax,'on'); grid(ax,'on');

for i=1:Nsections
    plot3(X1{i}(1,:),X1{i}(2,:),X1{i}(3,:),'k-');
    plot3(X2{i}(1,:),X2{i}(2,:),X2{i}(3,:),'r-');
    plot3(X3{i}(1,:),X3{i}(2,:),X3{i}(3,:),'b-');
    plot3(X4{i}(1,:),X4{i}(2,:),X4{i}(3,:),'g-');
end

axis(ax,'equal'); axis(ax,'tight');
view(ax,3);
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');

%% Plot as surface
% Make figure
f = figure('Color','w','Position',[100,100,700,400]);
ax = axes('Parent',f);
hold(ax,'on'); grid(ax,'on');

% note to self; need to remove last element of each section here

xsurf = cell2mat(xmat);
ysurf = cell2mat(ymat);
zsurf = cell2mat(zmat);

surface(xsurf(:,[1:4,1]),ysurf(:,[1:4,1]),zsurf(:,[1:4,1]),'FaceColor',[0.3,0.3,0.3],'FaceAlpha',0.5,'EdgeColor','none');
plot3(xsurf,ysurf,zsurf,'w-');

% surface(xsurf(:,[1:4,1]),ysurf(:,[1:4,1]),-zsurf(:,[1:4,1]),'FaceColor',[0.3,0.3,0.3],'FaceAlpha',0.5,'EdgeColor','none');
% plot3(xsurf,ysurf,-zsurf,'w-');

axis(ax,'equal'); axis(ax,'tight');
view(ax,3);
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');

autoArrangeFigures(2,3,1)


%% export coil to field
% create field model
mymodel = create_model();
mymodel.verbose = -1;
mymodel.device = 'CPU';
mymodel.waitbar = false;

% create coils as MP file class
mycoil = create_MP;
mycoil.scale = 1; % dimensions are in m
mycoil.name = 'Coil up';

% shape of coil
mycoil.X = cell2mat(xmat);
mycoil.Y = cell2mat(ymat);
mycoil.Z = cell2mat(zmat);
mycoil.nsec = size(mycoil.X,1);
mycoil.Across = CS;
mycoil.gamma = pi/4;
mycoil.theta = 0;

% number of elements in cross section
mycoil.nh = 10;
mycoil.nd = 20;

% current settings
mycoil.st = 'nJ'; % set number of windings and current density
mycoil.nwindings = 1;
mycoil.J = 88E+6; % [A/mm^2]
mycoil.calcnIJ();

% add coil to model
mymodel.add(mycoil);

%%
mycoil = create_MP;
mycoil.scale = 1; % dimensions are in m
mycoil.name = 'Coil down';

% shape of coil
mycoil.X = cell2mat(xmat);
mycoil.Y = cell2mat(ymat);
mycoil.Z = -cell2mat(zmat);
mycoil.nsec = size(mycoil.X,1);
mycoil.Across = CS ;
mycoil.gamma = pi/4;
mycoil.theta = 0;

% number of elements in cross section
mycoil.nh = 10;
mycoil.nd = 20;

% current settings
mycoil.st = 'nJ'; % set number of windings and current density
mycoil.nwindings = 1;
mycoil.J = 88E+6; % [A/mm^2]
mycoil.calcnIJ();

% add coil to model
mymodel.add(mycoil);

%% Start Field GUI
geo = tool_geometry;
geo.add(mymodel);
geo.updateModel;


%% ---- Peak Field Calculation ----
% mysurf = calculate_surface;
% mysurf.msinsource = mymodel;
% mysurf.msintarget = mymodel;
% mysurf.calculate;
%
% fprintf(['Peak field on the conductors is: ', num2str(max(mysurf.peakfield.Bmax)), ' [T] \n \n']);
% Bx_surf = reshape(mysurf.field{1}.B.x, 1, []);
% By_surf = reshape(mysurf.field{1}.B.y, 1, []);
% Bz_surf = reshape(mysurf.field{1}.B.z, 1, []);

%% ---- Field Plane Calculation ----
myplane = calculate_plane;
myplane.add(mymodel);
set(myplane,'planetype','xy','L',1.5,'H',1,'nstepL',500,'nstepH',500,'offsety',0.2,'offsetx',-0.6);
myplane.calculate();

%% plot plane calculation
figure;
Bmag = sqrt(myplane.B.x.^2 + myplane.B.y.^2 + myplane.B.z.^2);
surface(myplane.B.plotxlim,myplane.B.plotylim,myplane.B.z,'FaceColor','interp','EdgeColor','none');
colormap(jet); colorbar('EastOutside');
axis('equal'); axis('tight');

% interpolation of the field at coordinate
Bmagi=interp2(myplane.B.plotxlim,myplane.B.plotylim,Bmag,0,1);
