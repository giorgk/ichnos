%% Processor polys
c2vsim_outline = shaperead('F:\UCDAVIS\C2VSIM_FG_OR\C2Vsim_FG_v2\wrkspc\C2VsimMesh_Outline_3310.shp'); 
[Xs, Ys] = polysplit(c2vsim_outline.X, c2vsim_outline.Y);
simplify_threshold = 1000;
[xx,yy] = reducem(Xs{1,1}', Ys{1,1}', simplify_threshold);
%%
Dy = (max(yy) - min(yy))/4;
ymin = min(yy);
for ii = 1:4
    proc_poly(ii,1).x = [min(xx) max(xx) max(xx) min(xx)]';
    proc_poly(ii,1).y = [ymin ymin ymin+Dy ymin+Dy]';
    ymin = ymin + Dy;
end
%% Extended processor polygons
buffer = 5000;
for ii = 1:4
    ext_poly(ii,1).x = proc_poly(ii,1).x + [-buffer; buffer; buffer ; -buffer ];
    ext_poly(ii,1).y = proc_poly(ii,1).y + [-buffer; -buffer; buffer ; buffer ];
end
%%
clf
hold on
plot(BCX_3310, BCY_3310,'.')
plot(xx,yy,'linewidth',2)
for ii = 1:4
    plot(proc_poly(ii,1).x([1 2 3 4 1]), proc_poly(ii,1).y([1 2 3 4 1]),'--k','linewidth',1) 
    plot(ext_poly(ii,1).x([1 2 3 4 1]), ext_poly(ii,1).y([1 2 3 4 1]),':r','linewidth',1) 
end
axis equal
axis off
%% Write processor polygons
fid = fopen('c2vsim_4proc_polys.ich', 'w');
fprintf(fid,'%d\n', length(proc_poly));
for ii = 1:length(proc_poly)
   fprintf(fid, '%d %d\n', [ii-1 length(proc_poly(ii,1).x)]); 
   fprintf(fid, '%.2f %.2f\n', [proc_poly(ii,1).x proc_poly(ii,1).y]');
end
fclose(fid);
%% Write extended processor polygons
fid = fopen('c2vsim_4ext_polys.ich', 'w');
fprintf(fid,'%d\n', length(ext_poly));
for ii = 1:length(ext_poly)
   fprintf(fid, '%d %d\n', [ii-1 length(ext_poly(ii,1).x)]); 
   fprintf(fid, '%.2f %.2f\n', [ext_poly(ii,1).x ext_poly(ii,1).y]');
end
fclose(fid);
%% assign processor ids
for ii = 1:length(proc_poly)
    in = inpolygon(VEL_DATA(:,1), VEL_DATA(:,2), proc_poly(ii,1).x, proc_poly(ii,1).y);
    VEL_DATA(in,4) = ii - 1;
end
%% write data as ASCII files
for ii = 1:length(ext_poly)
    in = inpolygon(VEL_DATA(:,1), VEL_DATA(:,2), ext_poly(ii,1).x, ext_poly(ii,1).y);
    fid = fopen(['c2vsim_SS_05_15_4proc_' num2str(ii-1,'%04d') '.ich'],'w');
    fprintf(fid, '%.3f %.3f %.3f %d %.2f %.2f %.5f %.5f %.5f\n', VEL_DATA(in,:)');
    fclose(fid);
end
%% write data as HDF5
for ii = 1:length(ext_poly)
    in = inpolygon(VEL_DATA(:,1), VEL_DATA(:,2), ext_poly(ii,1).x, ext_poly(ii,1).y);
    h5create(['c2vsim_SS_05_15_4proc_' num2str(ii-1,'%04d') '.h5'],'/XYZDR',[sum(in) 5], 'Datatype','single');
    h5write(['c2vsim_SS_05_15_4proc_' num2str(ii-1,'%04d') '.h5'], '/XYZDR', VEL_DATA(in,[1 2 3 5 6]));
    h5create(['c2vsim_SS_05_15_4proc_' num2str(ii-1,'%04d') '.h5'],'/PROC',[sum(in) 1], 'Datatype','uint32');
    h5write(['c2vsim_SS_05_15_4proc_' num2str(ii-1,'%04d') '.h5'], '/PROC', VEL_DATA(in,4));
    h5create(['c2vsim_SS_05_15_4proc_' num2str(ii-1,'%04d') '.h5'],'/VXYZ',[sum(in) 3], 'Datatype','single');
    h5write(['c2vsim_SS_05_15_4proc_' num2str(ii-1,'%04d') '.h5'], '/VXYZ', VEL_DATA(in,[7 8 9]));
end
%% Generate particels
nParticles = 1000;
min_dist = 5000;
cnt = 0;
xmin = min(xx);
xmax = max(xx);
ymin = min(yy);
ymax = max(yy);
cv_shape = polyshape(c2vsim_outline.X, c2vsim_outline.Y);
particles = [];
while cnt < nParticles
   xr = xmin + (xmax - xmin)*rand; 
   yr = ymin + (ymax - ymin)*rand; 
   in = cv_shape.isinterior(xr, yr);
   if in
       if cnt == 0
           particles = [particles; xr yr];
           cnt = cnt + 1;
       else
           dst = pdist2(particles, [xr yr]);
           if min(dst) > min_dist
               particles = [particles; xr yr];
               cnt = cnt + 1;
           end
       end
   end
end
%% Calculate elevation
Ftop = scatteredInterpolant(X_3310, Y_3310, TopElev, 'linear');
pTop = Ftop(particles(:,1),particles(:,2));
particles(:,3) = pTop - 10;
%% Write scattered particles
fid = fopen('CV_scattered_particles.ich','w');
fprintf(fid, '# Scatter particles 10m below gse\n');
fprintf(fid, '#\n');
fprintf(fid, '%d %d %.3f %.3f %.3f\n', [ones(size(particles,1),1) [1:size(particles,1)]' particles]');
fclose(fid);
%% Read the printed velocity field
fid = fopen('c2vsim_SS_05_15_4proc_0000.ich','r');
tmp = textscan(fid,'%f %f %f %f %f %f %f %f %f\n');
fclose(fid);
plot(tmp{1,1},tmp{1,2},'.')
%% find velocity points inside the search box
ll = [-136535,138035,-12.8326];
uu = [-131097,143473,63.9158];
idx = find(VEL_DATA(:,1) > ll(1) & VEL_DATA(:,1) < uu(1) & ...
     VEL_DATA(:,2) > ll(2) & VEL_DATA(:,2) < uu(2) & ...
     VEL_DATA(:,3) > ll(3) & VEL_DATA(:,3) < uu(3));