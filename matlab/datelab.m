function strlab = datelab(ymd,UTsec)
% convert gemini time format to string
narginchk(2,2)
validateattr(ymd, {'numeric'}, {'vector', 'positive', 'numel', 3}, mfilename, 'year month day', 1)
validateattr(UTsec, {'numeric'}, {'scalar', 'nonnegative'}, mfilename, 'UTC second', 2)

%% SECONDS
if UTsec >= 86400
  error('0 <= UTsec < 86400')
end

% microsecond resolution
UTstr = num2str(UTsec, '%012.6f');
ymdnew=ymd;

strlab=[num2str(ymdnew(1), '%04d'), ...
        num2str(ymdnew(2), '%02d'), ...
        num2str(ymdnew(3), '%02d'), '_', UTstr];

end
