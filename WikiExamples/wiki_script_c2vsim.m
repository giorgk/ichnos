%% Domain files
c2vsim_outline = shaperead('F:\UCDAVIS\C2VSIM_FG_OR\C2Vsim_FG_v2\wrkspc\C2VsimMesh_Outline_3310.shp'); 
[Xs, Ys] = polysplit(c2vsim_outline.X, c2vsim_outline.Y);
%% Write the domain outline
simplify_threshold = 1000;
fid = fopen('c2vsimFG_outline.ich' ,'w');
for ii = 1:length(Xs)
    [xx,yy] = reducem(Xs{ii,1}', Ys{ii,1}', simplify_threshold);
    fprintf(fid, '%d %d\n',[length(xx) ispolycw(xx,yy)]);
    fprintf(fid, '%.3f %.3f\n', [xx yy]');
end
fclose(fid);
%% Write Top and bottom
%addpath(fullfile('..','..','GWToolsRpack','gwtools','Matlab'))
c2vsim_path = fullfile('..','..','C2VsimV1','c2vsim-working');
ND = readIWFM_Nodes(fullfile(c2vsim_path,'Preprocessor','C2VSimFG_Nodes.dat'));
% convert Nodes from EPSG 26910 to 3310
[lat,lon] = projinv(projcrs(26910),ND(:,1), ND(:,2));
[X_3310, Y_3310] = projfwd(projcrs(3310),lat, lon);

Strat = readIWFM_Stratigraphy(fullfile(c2vsim_path,'Preprocessor','C2VSimFG_Stratigraphy.dat'), size(ND,1), 105);
TopElev = Strat(:,2) * 0.3048; 
Bottom_elev = (Strat(:,2) - sum(Strat(:,3:end),2)) * 0.3048;

fid = fopen('c2vsim_TopBottom.ich','w');
fprintf(fid, '%.3f %.3f %.4f %.4f\n', [X_3310 Y_3310 TopElev Bottom_elev]');
fclose(fid);
%% Processor polygons for a single processor simulation.
offset = 1000;
proc_poly = [min(Xs{1,1}) - offset min(Ys{1,1}) - offset; ...
             min(Xs{1,1}) - offset max(Ys{1,1}) + offset; ...
             max(Xs{1,1}) + offset max(Ys{1,1}) + offset; ...
             max(Xs{1,1}) + offset min(Ys{1,1}) - offset];
fid = fopen('c2vsim_Single_proc_poly.ich', 'w');
fprintf(fid, '%d\n', 1);
fprintf(fid, '%d %d\n', [0 size(proc_poly,1)]);
fprintf(fid, '%.2f %.2f\n', proc_poly');
fclose(fid);
%% Read velocity
VelOut = readIWFM_Velocity(fullfile(c2vsim_path,'Results','C2VSimFG_GW_VELOUTFL.out'), 32537, 4);
%% Average velocity
for ii = 1:size(VelOut.VX,3)
    Vx_av(:,ii) = mean(VelOut.VX(:,385:504,ii),2);
    Vy_av(:,ii) = mean(VelOut.VY(:,385:504,ii),2);
    Vz_av(:,ii) = mean(VelOut.VZ(:,385:504,ii),2);
end
%% Convert the velocity units from AC-FT/month to m^3/day
days_per_month = mean(eomday(VelOut.YMD(385:504,1), VelOut.YMD(385:504,2)));
Vx_av = 1233.48 .* Vx_av ./ days_per_month;
Vy_av = 1233.48 .* Vy_av ./ days_per_month;
Vz_av = 1233.48 .* Vz_av ./ days_per_month;

%% Convert the velocity coordinates to EPSG.
[lat,lon] = projinv(projcrs(26910),VelOut.ND(:,1)*0.3048, VelOut.ND(:,2)*0.3048);
[BCX_3310, BCY_3310] = projfwd(projcrs(3310),lat, lon);

%% Calculate node elevations.
ND_ELEV(:,1) = Strat(:,2)*0.3048;
ND_ELEV(:,2) = ND_ELEV(:,1) - sum(Strat(:,3:4),2)*0.3048;
ND_ELEV(:,3) = ND_ELEV(:,2) - sum(Strat(:,5:6),2)*0.3048;
ND_ELEV(:,4) = ND_ELEV(:,3) - sum(Strat(:,7:8),2)*0.3048;
ND_ELEV(:,5) = ND_ELEV(:,4) - sum(Strat(:,9:10),2)*0.3048;
%%
msh = readIWFM_Elements(fullfile(c2vsim_path,'Preprocessor','C2VSimFG_Elements.dat'), 32537, 142);
%%
for ii = 1:size(msh,1)
    nv = 4;
    if msh(ii,5) == 0
        nv = 3;
    end
    for j = 1:4
        t = mean(ND_ELEV(msh(ii,2:nv+1),j));
        b = mean(ND_ELEV(msh(ii,2:nv+1),j+1));
        MSH_ELEV(ii,j) = (t + b)/2;
        av_thick(ii,j) =  mean(ND_ELEV(msh(ii,2:nv+1),j) - ND_ELEV(msh(ii,2:nv+1),j+1));
    end
    DIAM(ii,1) = max(pdist([X_3310(msh(ii,2:nv+1)) Y_3310(msh(ii,2:nv+1))]));
    
    %RATIO(ii,1) = DIAM(ii,1)/av_thick(ii,1); 
    
    av_dist(ii,1) = mean(sqrt((X_3310(msh(ii,2:nv+1))-BCX_3310(ii,1)).^2 + ...
        (Y_3310(msh(ii,2:nv+1))-BCY_3310(ii,1)).^2));
    elem_area(ii,1) = polyarea(X_3310(msh(ii,2:nv+1)), Y_3310(msh(ii,2:nv+1)));
end
%% Convert the velocity to m/day
Vx_av = Vx_av ./ (2.*av_dist.*av_thick);
Vy_av = Vy_av ./ (2.*av_dist.*av_thick);
Vz_av = Vz_av ./ elem_area;
%% Reshare tha data
nvel_p = size(MSH_ELEV,1)*size(MSH_ELEV,2);
mult = 1000000;
VEL_DATA = [...
    repmat(BCX_3310,4,1) ... % X
    repmat(BCY_3310,4,1) ... % X
    reshape(MSH_ELEV, nvel_p, 1) ... % Z
    zeros(nvel_p,1) ... % Proc ID
    reshape(Vx_av, nvel_p, 1)*mult ... % VX
    reshape(Vy_av, nvel_p, 1)*mult ... % VY
    reshape(Vz_av, nvel_p, 1)*mult ... % VZ
    repmat(DIAM,4,1) ... % DIAM 
    reshape(bsxfun(@rdivide, DIAM, av_thick), nvel_p, 1) ... % RATIO
];
%% Write velocity file
fid = fopen('c2vsim_SS_05_15_vel_000.ich','w');
fprintf(fid, '%.3f %.3f %.3f %d %.5f %.5f %.5f 0 %.2f %.2f\n', VEL_DATA');
fclose(fid);
%% Generate Individual particles
bbx = c2vsim_outline.BoundingBox(:,1)';
bby = c2vsim_outline.BoundingBox(:,2)';
particles_locations = [];
for ii = 1:20
    while 1
        xr = bbx(1) + (bbx(2) - bbx(1))*rand;
        yr = bby(1) + (bby(2) - bby(1))*rand;
        if inpolygon(xr, yr, Xs{1,1}, Ys{1,1})
            isin = false;
            for k = 2:4
               if inpolygon(xr, yr, Xs{k,1}, Ys{k,1})
                   isin = true;
                  break 
               end
            end
            if ~isin
                particles_locations = [particles_locations; xr yr];
               break 
            end
        end
    end
end
%%
plot(c2vsim_outline.X, c2vsim_outline.Y)
hold on
plot(particles_locations(:,1),particles_locations(:,2),'.')
%% Manually set the particle locations
plot(c2vsim_outline.X, c2vsim_outline.Y)
particles_locations = ginput(20);
%% Create an interpolation function for the top elevation
Ftop = scatteredInterpolant(X_3310, Y_3310, TopElev, 'linear');
pTop = Ftop(particles_locations(:,1),particles_locations(:,2));
%% 
particleData = [ ...
    10*ones(length(pTop),1) (1:length(pTop))' particles_locations pTop - 10];
particleData = [particleData; ...
    20*ones(length(pTop),1) (1:length(pTop))' particles_locations pTop - 20];
particleData = [particleData; ...
    30*ones(length(pTop),1) (1:length(pTop))' particles_locations pTop - 30];
particleData = [particleData; ...
    40*ones(length(pTop),1) (1:length(pTop))' particles_locations pTop - 40];
%% Write particles
fid = fopen('CV_particles.ich','w');
fprintf(fid, '# Test particle file\n');
fprintf(fid, '#\n');
fprintf(fid, '%d %d %.3f %.3f %.3f\n', particleData');
fclose(fid);
%% Read the particles
S = readICHNOStraj('c2vsim_out01__ireal_0000_iter_0000_proc_0000.traj');
%% plot streamlines
clf
hold on
for ii = 1:length(S)
    plot3(S(ii,1).p(:,1), S(ii,1).p(:,2), S(ii,1).p(:,3),'-')
end
plot(c2vsim_outline.X, c2vsim_outline.Y, 'linewidth',2)
axis equal
axis off
%% Prepare input files for Transient State velocity
VelOut = readIWFM_Velocity(fullfile(c2vsim_path,'Results','C2VSimFG_GW_VELOUTFL.out'), 32537, 4);
%% Convert the velocity units from AC-FT/month to m^3/day
days_per_month = eomday(VelOut.YMD(:,1), VelOut.YMD(:,2));
Vx_m3_day = VelOut.VX;
Vy_m3_day = VelOut.VY;
Vz_m3_day = VelOut.VZ;
for ii = 1:length(days_per_month)
    ii
    Vx_m3_day(:,ii,:) = 1233.48 .* Vx_m3_day(:,ii,:) ./ days_per_month(ii);
    Vy_m3_day(:,ii,:) = 1233.48 .* Vy_m3_day(:,ii,:) ./ days_per_month(ii);
    Vz_m3_day(:,ii,:) = 1233.48 .* Vz_m3_day(:,ii,:) ./ days_per_month(ii);
end
%% Convert the velocity to m/day
Vx_m_day = Vx_m3_day;
Vy_m_day = Vy_m3_day;
Vz_m_day = Vz_m3_day;
for ii = 1:size(Vx_m_day,2)
    for jj = 1:size(Vx_m_day,3)
        Vx_m_day(:,ii,jj) = Vx_m_day(:,ii,jj) ./ (2.*av_dist.*av_thick(:,jj));
        Vy_m_day(:,ii,jj) = Vx_m_day(:,ii,jj) ./ (2.*av_dist.*av_thick(:,jj));
        Vz_m_day(:,ii,jj) = Vx_m_day(:,ii,jj) ./ elem_area;
    end
end
%% Split velocity data
% The velocity files will be printed into multiple files for easier
% manipulation and debug
%% XYZ files
nvel_p = size(MSH_ELEV,1)*size(MSH_ELEV,2);
XYZ_DATA = [...
    repmat(BCX_3310,4,1) ... % X
    repmat(BCY_3310,4,1) ... % X
    reshape(MSH_ELEV, nvel_p, 1) ... % Z
    zeros(nvel_p,1) ... % Proc ID
    repmat(DIAM,4,1) ... % DIAM 
    reshape(bsxfun(@rdivide, DIAM, av_thick), nvel_p, 1) ... % RATIO
];
%% write the XYZ file
prefix = 'c2vsim_TR_DBG_';
fid = fopen([prefix 'XYZ_000.ich'],'w');
fprintf(fid, '%.3f %.3f %.3f %d %.2f %.2f\n', XYZ_DATA');
fclose(fid);
%% prepare the velocity files
VX_DATA = [];
VY_DATA = [];
VZ_DATA = [];
for ii = 1:size(Vx_m_day,3)
    VX_DATA = [VX_DATA; Vx_m_day(:,:,ii)];
    VY_DATA = [VY_DATA; Vy_m_day(:,:,ii)];
    VZ_DATA = [VZ_DATA; Vz_m_day(:,:,ii)];
end
%% extract a small subset of the data for debugging
plot(XYZ_DATA(:,1),XYZ_DATA(:,2),'.')
dbg_domain = ginput;
hold on
plot(dbg_domain(:,1), dbg_domain(:,2))
%% create domain outline for debug area
simplify_threshold = 1000;
[xx,yy] = reducem(Xs{1,1}', Ys{1,1}', simplify_threshold);
poly1 = polyshape(xx, yy);
poly2 = polyshape(dbg_domain(:,1), dbg_domain(:,2));
polyout = intersect(poly1,poly2);
%% Print dbg domain outline
fid = fopen('c2vsimFG_outlineDBG.ich' ,'w');
fprintf(fid, '%d %d\n',[size(polyout.Vertices,1) 1]);
fprintf(fid, '%.3f %.3f\n', polyout.Vertices');
fclose(fid);
%% extract the indices of the point inside the dbg domain
idx = find(inpolygon(XYZ_DATA(:,1),XYZ_DATA(:,2),dbg_domain(:,1), dbg_domain(:,2)));
plot(XYZ_DATA(idx,1),XYZ_DATA(idx,2),'o')
XYZ_DATA = XYZ_DATA(idx,:);
VX_DATA = VX_DATA(idx,:);
VY_DATA = VY_DATA(idx,:);
VZ_DATA = VZ_DATA(idx,:);

%% write the VEL files
mult = 1000000;
frmt = ['%.5f' repmat(' %.5f',1,size(VX_DATA,2)-1) '\n'];

% VX
fid = fopen([prefix 'VX_000.ich'],'w');
fprintf(fid, frmt, mult*VX_DATA');
fclose(fid);
% VY
fid = fopen([prefix 'VY_000.ich'],'w');
fprintf(fid, frmt, mult*VY_DATA');
fclose(fid);
% VZ
fid = fopen([prefix 'VZ_000.ich'],'w');
fprintf(fid, frmt, mult*VZ_DATA');
fclose(fid);
%% Time stamp
fid = fopen('c2vsim_STEP_monthly.ich','w');
fprintf(fid, '%.0f\n', [0; cumsum(days_per_month(1:length(days_per_month)-1))]);
fclose(fid);
%% create a few particles
plot(polyout)
part_TR = ginput;
part_TR(:,3) = Ftop(part_TR(:,1),part_TR(:,2)) - 50;
%% Write particles
fid = fopen('CV_particlesTR.ich','w');
fprintf(fid, '# Test particle file in transient velocity field\n');
fprintf(fid, '# Transient tracing requires the location and the release time in the field\n');
fprintf(fid, '# Eid Sid X Y Z Time\n');
fprintf(fid, '#\n');
fprintf(fid, '%d %d %.3f %.3f %.3f %.2f\n', [ones(size(part_TR,1),1) (1:size(part_TR,1))' part_TR 8035*ones(size(part_TR,1),1)]');
fclose(fid);
