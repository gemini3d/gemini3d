function f=qp2robj(r,q,p)

Re=6370e3;
f=abs(q^2*(r/Re)^4+1/p*(r/Re)-1);

end