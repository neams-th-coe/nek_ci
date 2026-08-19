[Mesh]
  type = NekRSMesh
  volume = true
[]

[Problem]
  type = NekRSProblem
  casename = 'ethier'
[]

[Executioner]
  type = Transient

  [TimeStepper]
    type = NekTimeStepper
  []
[]

[Outputs]
  csv = true
  show = 'pass'
  execute_on = final
  console = false
  file_base = 'nek_out'
  #[console]
  #  type = Console
  #  time_step_interval = 100
  #[]
[]

P_U0 = 0.5
P_V0 = 0.1
P_W0 = 0.2
P_A0 = 0.025
P_D0 = 0.5
a = ${fparse pi * P_A0}
nu = ${fparse 1/100}
d = ${fparse pi * P_D0}

TOL_V  = 4.78E-08
TOL_P  = 1.08E-07
TOL_T  = 9.48E-10
TOL_S  = 1.07E-09
TOL    = 1.00E-11


[Functions]
  [unitFunction]
    type = ParsedFunction
    expression = '1.0'
  []
  [uexact]
    type = ParsedFunction
    expression = '-${a} * (ex * syz + ez * cxy) * e2t + ${P_U0}'
    symbol_names = 'ex syz ez cxy e2t'
    symbol_values = 'ex syz ez cxy e2t'
  []
  [pexact]
    type = ParsedFunction
    expression = '-0.5 * ${a} * ${a} * e2t * e2t * (ex * ex + 2.0 * sxy * czx * eyz
                                              +ey * ey + 2.0 * syz * cxy * ezx
                                              +ez * ez + 2.0 * szx * cyz * exy)'
    symbol_names = 'e2t ex ey ez sxy syz szx cxy cyz czx exy eyz ezx'
    symbol_values = 'e2t ex ey ez sxy syz szx cxy cyz czx exy eyz ezx'
  []
  [pexactScaled]
    type = ParsedFunction
    expression = 'pexact + pscale'
    symbol_names = 'pexact pscale'
    symbol_values = 'pexact pscale'
  []

  # we create some intermediate functions to make the formulas above a bit more
  # human-readable (and to match the intermediate variables used in the original
  # NekRS test)
  [xx]
    type = ParsedFunction
    expression = 'x - ${P_U0} * t'
  []
  [yy]
    type = ParsedFunction
    expression = 'y - ${P_V0} * t'
  []
  [zz]
    type = ParsedFunction
    expression = 'z - ${P_W0} * t'
  []
  [ex]
    type = ParsedFunction
    expression = 'exp(${a} * xx)'
    symbol_names = 'xx'
    symbol_values = 'xx'
  []
  [ey]
    type = ParsedFunction
    expression = 'exp(${a} * yy)'
    symbol_names = 'yy'
    symbol_values = 'yy'
  []
  [ez]
    type = ParsedFunction
    expression = 'exp(${a} * zz)'
    symbol_names = 'zz'
    symbol_values = 'zz'
  []
  [exy]
    type = ParsedFunction
    expression = 'exp(${a} * (xx + yy))'
    symbol_names = 'xx yy'
    symbol_values = 'xx yy'
  []
  [eyz]
    type = ParsedFunction
    expression = 'exp(${a} * (yy + zz))'
    symbol_names = 'yy zz'
    symbol_values = 'yy zz'
  []
  [ezx]
    type = ParsedFunction
    expression = 'exp(${a} * (zz + xx))'
    symbol_names = 'zz xx'
    symbol_values = 'zz xx'
  []
  [e2t]
    type = ParsedFunction
    expression = 'exp(-${nu} * ${d} * ${d} * t)'
  []
  [sxy]
    type = ParsedFunction
    expression = 'sin(${a} * xx + ${d} * yy )'
    symbol_names = 'xx yy'
    symbol_values = 'xx yy'
  []
  [syz]
    type = ParsedFunction
    expression = 'sin(${a} * yy + ${d} * zz )'
    symbol_names = 'zz yy'
    symbol_values = 'zz yy'
  []
  [szx]
    type = ParsedFunction
    expression = 'sin(${a} * zz + ${d} * xx )'
    symbol_names = 'xx zz'
    symbol_values = 'xx zz'
  []
  [cxy]
    type = ParsedFunction
    expression = 'cos(${a} * xx + ${d} * yy )'
    symbol_names = 'xx yy'
    symbol_values = 'xx yy'
  []
  [cyz]
    type = ParsedFunction
    expression = 'cos(${a} * yy + ${d} * zz )'
    symbol_names = 'zz yy'
    symbol_values = 'zz yy'
  []
  [czx]
    type = ParsedFunction
    expression = 'cos(${a} * zz + ${d} * xx )'
    symbol_names = 'xx zz'
    symbol_values = 'xx zz'
  []
[]

[Postprocessors]
  # Calculate parameters for exact solution
  [volume]
    type = NekVolumeIntegral
    field = unity
    function = unitFunction
    execute_on = final
  []
  [pbar]
    type = NekVolumeIntegral
    field = pressure
    execute_on = final
  []
  [pbre]
    type = NekVolumeIntegral
    field = unity
    function = pexact
    execute_on = final
  []
  [pscale]
    type = ParsedPostprocessor
    expression = '(pbar - pbre)/volume'
    pp_names = 'pbre pbar volume'
    execute_on = final
  []
  
  # Calculate L2 errors
  [uxerrl2]
    type = NekVolumeNorm
    field = velocity_x
    function = uexact
    execute_on = final
  []
  [terrl2]
    type = NekVolumeNorm
    field = scalar01
    function = uexact
    execute_on = final
  []
  [serrl2]
    type = NekVolumeNorm
    field = scalar02
    function = uexact
    execute_on = final
  []
  [perrl2]
    type = NekVolumeNorm
    field = pressure
    function = pexactScaled
    execute_on = final
  []
  
  # Check if all tests passed
  [pass]
    type = ParsedPostprocessor
    expression = 'if((uxerrl2 < ${TOL_V} | uxerrl2 < ${TOL}) &
                     (perrl2 < ${TOL_P} | perrl2 < ${TOL}) &
                     (terrl2 < ${TOL_T} | terrl2 < ${TOL}) &
                     (serrl2 < ${TOL_S} | serrl2 < ${TOL}), 1, 0)'
    pp_names = 'uxerrl2 perrl2 terrl2 serrl2'
    execute_on = final
  []
[]
