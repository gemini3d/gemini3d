function newdata=interpolate(data,grid,newgrid,exmode,imode)

%function newdata=interpolate(data,grid,newgrid,exmode,imode)
%
%(ex,i)mode controls the extrapolation behavior:
%     exmode='hold' --> hold values for newgrid < grid(1) at data(1)
%                hold values for newgrid > grid(length(grid)) at data(length(data))
%     exmode='login' --> log linear extrapolation at high and low ends
%     exmode='lin'  --> linear extrapolation at both ends
%
%     imode='loglin'
%     imode='lin'

       switch exmode
         case 'hold',	exflag=0;
         case 'loglin',	exflag=1;
         case 'lin',	exflag=2;
         otherwise,
           fprintf('\n INTERPOLATE --> exmode not valid!')
       end
       switch imode
         case 'loglin',	inflag=1;      
         case 'lin',	inflag=2;
         otherwise,
           fprintf('\n INTERPOLATE --> imode not valid!')  
       end


  lgr=length(grid);
  newdata=zeros(size(newgrid));
  l=1;
  for k=1:length(newgrid)  
    %Compute the bin # that the new point falls in
    while  l<=length(grid) && newgrid(k)>grid(l)
       l=l+1; 
    end
    
    %Low end extrap.
    if l==1
       switch exflag
         case 0,
           newdata(k)=data(1);
         case 1,
           x0=log10(grid(1)); x1=log10(grid(2));
           y0=log10(data(1)); y1=log10(data(2));
           slope=(y1-y0)/(x1-x0);
           newdata(k)=10^(y0+slope*(log10(newgrid(k))-x0));
         case 2,
           x0=grid(1); x1=grid(2);
           y0=data(1); y1=data(2);
           slope=(y1-y0)/(x1-x0);
           newdata(k)=y0+slope*(newgrid(k)-x0);
         otherwise,
           fprintf('\n INTERPOLATE --> exmode not valid!')
       end
    %High end extrap.
    elseif l>lgr
       switch exflag
         case 0,
           newdata(k)=data(length(data));
         case 1,
           x0=log10(grid(lgr-1)); x1=log10(grid(lgr));
           y0=log10(data(lgr-1)); y1=log10(data(lgr));
           slope=(y1-y0)/(x1-x0);
           newdata(k)=10^(y0+slope*(log10(newgrid(k))-x0));
         case 2,
           x0=grid(lgr-1); x1=grid(lgr);
           y0=data(lgr-1); y1=data(lgr);
           slope=(y1-y0)/(x1-x0);
           newdata(k)=y0+slope*(newgrid(k)-x0);
         otherwise,
           fprintf('\n INTERPOLATE --> exmode not valid!')
       end
    %Interp.
    else
       switch inflag
         case 1,
           x0=log10(grid(l-1)); x1=log10(grid(l));
           y0=log10(data(l-1)); y1=log10(data(l));
           slope=(y1-y0)/(x1-x0);
           newdata(k)=10^(y0+slope*(log10(newgrid(k))-x0));          
         case 2,
           x0=grid(l-1); x1=grid(l);
           y0=data(l-1); y1=data(l);
           slope=(y1-y0)/(x1-x0);
           newdata(k)=y0+slope*(newgrid(k)-x0);
         otherwise,
           fprintf('\n INTERPOLATE --> imode not valid!')  
       end
    end %conditional
    
    l=1;
  end %loop over data points
  
%  newdata=reshape(newdata,size(newgrid));
end %function
