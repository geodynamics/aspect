format long;
[radius]=textread('layers-radius.txt','%f',28,'headerlines',0);
knots=zeros(28,1);
for i=1:28
    knots(i)=-1+((radius(i)-3480)/(6346-3480))*2;
end
fid=fopen('Spline_knots.txt','w');
fprintf(fid,'%f\n',knots);
fclose(fid);
