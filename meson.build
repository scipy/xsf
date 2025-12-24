project(
  'xsf',
  'cpp',
  version : '0.1.3',
  license : 'BSD-3-Clause AND MIT AND BSD-3-Clause-LBNL AND Apache-2.0 WITH LLVM-exception',
  license_files : ['LICENSE', 'LICENSES_bundled.txt'],
  meson_version : '>=1.5.0'
)


# Find dependencies
lapack = dependency('lapack', required : true)
blas = dependency('blas', required : true)

_args_blas_lapack = []  # Your compile args

lapack_lp64_dep = declare_dependency(
  dependencies : [lapack, blas],
  compile_args : _args_blas_lapack
)

# For header-only library
xsf_inc = include_directories('include')

xsf_dep = declare_dependency(
  include_directories : xsf_inc,
  dependencies : lapack_lp64_dep
)