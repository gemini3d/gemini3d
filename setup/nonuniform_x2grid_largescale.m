minx2=-1000e3;
maxx2=1000e3;
x2ref1=-300e3;
x2ref2=300e3;
sigx2=25e3;


%FIRST COUNT NUMBER OF ELEMENTS NEEDED
ix2=1;
x2tmp=0;
while(x2tmp<maxx2)
  dx2=20e3-16.667e3/2e0*(tanh((x2tmp-x2ref1)/sigx2)-tanh((x2tmp-x2ref2)/sigx2));
  x2tmp=x2tmp+dx2;
  ix2=ix2+1;
end
lx2=(ix2-1)*2+1;

x2=zeros(1,lx2);

%START AT ZERO AND WORK OUT, INSURING THAT THE GRID IS SYMMETRIC/CENTERED PROPERLY
%FILL POSITIVE HALF FIRST
ix2=floor(lx2/2)+1;
x2(ix2)=0;
while(x2(ix2)<maxx2)
  dx2=20e3-16.667e3/2e0*(tanh((x2(ix2)-x2ref1)/sigx2)-tanh((x2(ix2)-x2ref2)/sigx2));
  x2(ix2+1)=x2(ix2)+dx2;
  ix2=ix2+1;
end

%NOW JUST MIRROR GRID TO GET NEGATIVE HALF
for (ix2=1:floor(lx2/2))
  x2(ix2)=-1*x2(lx2+1-ix2);
end


minx3=-1000e3;
maxx3=1000e3;
x3ref1=-200e3;
x3ref2=200e3;
sigx3=25e3;


%FIRST COUNT NUMBER OF ELEMENTS NEEDED
ix3=1;
x3tmp=0;
while(x3tmp<maxx3)
  dx3=20e3-18e3/2e0*(tanh((x3tmp-x3ref1)/sigx3)-tanh((x3tmp-x3ref2)/sigx3));
  x3tmp=x3tmp+dx3;
  ix3=ix3+1;
end
lx3=(ix3-1)*2+1;

x3=zeros(1,lx3);

%START AT ZERO AND WORK OUT, INSURING THAT THE GRID IS SYMMETRIC/CENTERED PROPERLY
%FILL POSITIVE HALF FIRST
ix3=floor(lx3/2)+1;
x3(ix3)=0;
while(x3(ix3)<maxx3)
  dx3=20e3-18e3/2e0*(tanh((x3(ix3)-x3ref1)/sigx3)-tanh((x3(ix3)-x3ref2)/sigx3));
  x3(ix3+1)=x3(ix3)+dx3;
  ix3=ix3+1;
end

%NOW JUST MIRROR GRID TO GET NEGATIVE HALF
for (ix3=1:floor(lx3/2))
  x3(ix3)=-1*x3(lx3+1-ix3);
end



% x1min=80e3;
% x1max=1000e3;
% x1(1)=x1min;
% ix1=1;
% while(x1(ix1)<x1max)
%   ix1=ix1+1;
%   dx1tmp=30d3+25d3*tanh((x1(ix1-1)-400d3)/150d3);
%   x1(ix1)=x1(ix1-1)+dx1tmp;
% end
% lx1=numel(x1);


%DO SOME PLOTS
%figure, plot(diff(x1)/1e3,x1(1:end-1)/1e3,'o'), xlabel('dx1 (km)'), ylabel('x1 (km)')
figure, plot(x2(1:end-1)/1e3,diff(x2)/1e3,'o'), xlabel('dx2 (km)'), ylabel('x2 (km)')
figure, plot(x3(1:end-1)/1e3,diff(x3)/1e3,'o'), xlabel('dx3 (km)'), ylabel('x3 (km)')