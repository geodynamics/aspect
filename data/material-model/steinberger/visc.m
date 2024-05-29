mat=[0.186E+21 0.5448
0.354E+22 0.5668
0.332E+23 0.5888
0.885E+23 0.6107
0.112E+24 0.6327
0.111E+24 0.6547
0.101E+24 0.6767
0.875E+23 0.6986
0.731E+23 0.7206
0.586E+23 0.7426
0.451E+23 0.7646
0.332E+23 0.7865
0.235E+23 0.8085
0.159E+23 0.8305
0.103E+23 0.8525
0.624E+22 0.8744
0.651E+21 0.8964
0.431E+21 0.9184
0.432E+21 0.9356
0.290E+21 0.9466
0.241E+21 0.9655
0.264E+23 0.9843
];

% values are for the middle of the interval, so shift the x values

len = size(mat,1);
for i=1:size(mat,1)-1
    mat(i,2)=mat(i,2) + 0.5*(mat(i+1,2)-mat(i,2));
    
end
% last point is special
mat(len,2)=mat(len,2) + 0.5*(1.0-mat(len,2));
    
%interpolate
range=mat(1,2):0.0005:1.0;
yi = interp1(mat(:,2), mat(:,1), range,'linear','extrap');

plot(mat(:,2), mat(:,1),'x')
hold on
plot(range, yi, 'r-')



range = (1.0-range)*6371;

format longeng
A=[yi(length(yi):-1:1)' range(length(range):-1:1)']

save radial-visc.txt A -ASCII
