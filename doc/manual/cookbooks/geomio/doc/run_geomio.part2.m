% define the bounding box for the output mesh
% (this should be the X extent and Y extent in your ASPECT model)
xmin = 0; xmax = opt.svg.width;
ymin = 0; ymax = opt.svg.height;

% set the resolution in the output file: 
% [Xp,Yp] = ndgrid(xmin:your_steplength_x:xmax,ymin:your_steplength_y:ymax);
[Xp,Yp]      = ndgrid(xmin:15:xmax,ymin:15:ymax);
Phase        = zeros(size(Xp));

% assign a phase to each grid point according to your drawing
Phase = assignPhase2Markers(PathCoord, opt, Xp, Yp, Phase);

% plot your output
figure(2)
scatter(Xp(:),Yp(:),10,Phase(:),'filled');
axis equal
axis([xmin xmax ymin ymax])
