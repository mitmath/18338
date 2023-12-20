t=sym("t","real");
s=sym("s","real");
g=[s^3+t^2*s;-2*s*t;s];
dg=[[t^2-s^2;-2*t;1],[-2*s*t;2*s;0]];
dftilde=(g'*g)*[[1,0,0];[0,1,0];[0,0,1]]- g*g'; %unnormalised
df=dftilde ./ (g'*g)^1.5;
dX3=df*dg; % 3by3 * 3by2 = 3by2
dS3=norm(cross(dX3(:,1),dX3(:,2))); %symbolic area element

dS3_handle=@(t,s) double(subs(dS3)); %reads dS3 into a double
%% 
area=integral2(dS3_handle,-Inf,Inf,-Inf,Inf);%,RelTol=1e-9); 

%%
c=[1;t;t^2];
dc=[0;1;2*t];
dhtilde=(c'*c)*eye(3)-c*c';
dh=dhtilde ./ (c'*c)^1.5;
dl=norm(dh*dc);
dl_handle=@(t) double(subs(dl));

length=integral(dl_handle,-Inf,Inf);

%l correponds to pi*E(#real roots)
%a corresponds to 2*pi*E(#complex roots)
length/pi + area/2/pi