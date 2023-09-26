sources = files(
  'Precision.F90',
  'Constants.f90',
  'Histogram.F90',
  'KernelMultiGaussian.f90',
  'GridCell.f90', 
  'GridProjectedKDE.F90',
  'compilerversion.fpp', 
  'openspec.f90', 
  'utl7.f',
  'utl8.f90', 
)

core   = static_library('core',
                      sources,  
                      dependencies: dependencies)
theexe = executable('gpkde', 
                    'GPKDE.F90', 
                    link_with: [core],  
                    dependencies: dependencies, 
                    install: true)
