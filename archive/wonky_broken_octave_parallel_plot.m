  if isoctave
    for i = 1:Nt
      cmd = ['octave --eval "plotframe(''',direc,''',[',int2str(ymd(i,:)),'],',num2str(UTsec(i))];
      
      if ~isempty(saveplots)
        cmd = [cmd,",'",saveplots,"')"""]; %#ok<AGROW>
      else
        cmd = [cmd,')"']; %#ok<AGROW>
      end
      disp(cmd)

      % set to "sync" for debugging
      system(cmd, false, "async");
      % don't overload system RAM
      pause(2)
      ramfree = memfree();
      while ramfree < 1e9
        disp(['waiting for enough RAM to plot, you have MB free: ',num2str(ramfree/1e6,'%7.1f')])
        pause(10)
      end

    end