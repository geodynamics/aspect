% the headers ASPECT needs for the ascii data plugin
header1 = 'x';
header2 = 'y';
header3 = 'phase';

% create an array in the correct format for the ascii data plugin
Vx = Xp(:);
Vy = Yp(:);
VPhase = Phase(:);
[m,n] = size(Phase);

% write the data into the output file
fid=fopen('jelly.txt','w');
fprintf(fid, '# POINTS: %d %d \n',[m n]);
fprintf(fid, ['# Columns: ' header1 ' ' header2 ' ' header3 '\n']);
fprintf(fid, '%f %f %f \n', [Vx Vy VPhase]');
fclose(fid);
