submodule (ionization) glow_dummy

implicit none (type, external)

contains

module procedure glow_run
  ! these parameters are never used due to error stop.
  ! They just stop lint output for stringent compiler warning flags
  ionrate = 0.
  eheating = 0.
  iver = 0.
  error stop 'There appears to be a misconfiguration in your build'
end procedure glow_run

end submodule glow_dummy
