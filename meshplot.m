function meshplot(T, nx, ny)
figure('Name','Temperature Mesh Plot','NumberTitle','off');
%[M,N] = size(T);
[x,y] = meshgrid(1:ny,1:nx);
surf(x,y,T);
 
end