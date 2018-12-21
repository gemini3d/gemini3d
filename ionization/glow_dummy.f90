submodule (ionization) glow_dummy

contains

module procedure glow_run
  error stop 'There appears to be a misconfiguration in your build'
end procedure glow_run

end submodule glow_dummy
