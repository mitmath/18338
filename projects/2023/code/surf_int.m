
f=@(t,s) trivialweight(t,s);
area=quad2d(f,-1000,1000,-1000,1000)

g=@(t,s) areaweight_2dcomplex(t,s);
area=quad2d(g,1,1000,1,1000)