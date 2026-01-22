%%Requires dS3(t,s) already read into workspace
x=linspace(-2,2,78);
y=linspace(-2,2,78);
[X,Y]=meshgrid(x,y);

Z=zeros(size(X));
for i=1:size(Z,1)
    for j=1:size(Z,2)
        t=X(i,j);
        s=Y(i,j);
        Z(i,j)=subs(dS3);
    end
end
%%
figure()
pcolor(Z)
figure()
image('XData',x,'YData',y,'CData',Z/max(Z)*250)
