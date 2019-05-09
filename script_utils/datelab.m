function strlab=datelab(ymd,UTsec)

narginchk(2,2)
validateattr(ymd, {'numeric'}, {'vector', 'positive', 'numel', 3}, mfilename, 'year month day', 1)
validateattr(UTsec, {'numeric'}, {'scalar', 'positive'}, mfilename, 'UTC second', 2)


%SECONDS
if UTsec==86400
    UTstr='00000';
    [ymdnew,tmp]=dateinc(1,ymd,UT);
else
    UTstr=sprintf('%.6f',UTsec);    %print the UT seconds with 6 decimal points to a string variable

    UTstrtmp=num2str(round(UTsec));
    lstr=numel(UTstrtmp);
    while lstr<5                    %do some zero padding
        UTstr=['0',UTstr];
        lstr=lstr+1;
    end
    ymdnew=ymd;
end


%YEAR
yearstr=num2str(ymdnew(1));
lstr=numel(yearstr);
while lstr<4                        %zero padding
    yearstr=['0',yearstr];
    lstr=lstr+1;
end


%MONTH
monthstr=num2str(ymdnew(2));
lstr=numel(monthstr);
while lstr<2                        %zero padding
    monthstr=['0',monthstr];
    lstr=lstr+1;
end


%DAY
daystr=num2str(ymdnew(3));
lstr=numel(daystr);
while lstr<2                        %padding
    daystr=['0',daystr];
    lstr=lstr+1;
end

strlab=[yearstr,monthstr,daystr,'_',UTstr];

end
