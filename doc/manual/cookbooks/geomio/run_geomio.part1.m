% set options for geomIO
opt                 = geomIO_Options();
opt.inputFileName   = ['/path/to/aspect/doc/manual/cookbooks/geomio/jellyfish.svg'];
opt.DrawCoordRes    = 21; % optionally change resolution with opt.DrawCoordRes = your value;

% run geomIO
[PathCoord]         = run_geomIO(opt,'2D');

