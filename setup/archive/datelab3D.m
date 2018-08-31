function strlab=datelab3D(ymd,UTsec)


if UTsec==86400
    UTstr='00000';
    [tmp,dmynew]=dateinc(UT,dmy,1);
else
    UTstr=num2str(round(UTsec));
    lstr=numel(UTstr);
    while lstr<5
        UTstr=['0',UTstr];
        lstr=lstr+1;
    end
    ymdnew=ymd;
end

yearstr=num2str(ymdnew(1));
lstr=numel(yearstr);
while lstr<4
    yearstr=['0',yearstr];
    lstr=lstr+1;
end

monthstr=num2str(ymdnew(2));
lstr=numel(monthstr);
while lstr<2
    monthstr=['0',monthstr];
    lstr=lstr+1;
end

daystr=num2str(ymdnew(3));
lstr=numel(daystr);
while lstr<2
    daystr=['0',daystr];
    lstr=lstr+1;
end

strlab=[yearstr,monthstr,daystr,'_',UTstr,'.000000.dat'];

end
