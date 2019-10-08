function txt = time2str(ymd, UTsec)

narginchk(2,2)
t = datenum(ymd(1), ymd(2), ymd(3), 0, 0, UTsec);
txt = {datestr(t,1), [datestr(t,13),' UT']};

end