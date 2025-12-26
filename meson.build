project(
  'xsf',
  'cpp',
  version : '0.1.5',
  license : 'BSD-3-Clause AND MIT AND BSD-3-Clause-LBNL AND Apache-2.0 WITH LLVM-exception',
  license_files : ['LICENSE', 'LICENSES_bundled.txt'],
  meson_version : '>=1.5.0'
)



# Get compiler
cpp = meson.get_compiler('cpp')

# BLAS and LAPACK
# When used as a SciPy subproject, these are already set up by SciPy
# When built standalone, we need to find them

blas_name = get_option('blas')
lapack_name = get_option('lapack')

if blas_name == ''
  # Auto-detect
  blas = dependency('blas', required: false)
  if not blas.found()
    blas = cpp.find_library('blas', required: false)
  endif
  if not blas.found()
    blas = dependency('openblas', required: false)
  endif
  if not blas.found()
    blas = cpp.find_library('openblas', required: true)
  endif
else
  # Use specified BLAS
  blas = dependency(blas_name, required: true)
endif

if lapack_name == ''
  # Auto-detect
  lapack = dependency('lapack', required: false)
  if not lapack.found()
    lapack = cpp.find_library('lapack', required: false)
  endif
  if not lapack.found()
    lapack = dependency('openblas', required: false)
  endif
  if not lapack.found()
    lapack = cpp.find_library('openblas', required: true)
  endif
else
  # Use specified LAPACK
  lapack = dependency(lapack_name, required: true)
endif

message('xsf: BLAS library: ' + blas.type_name())
message('xsf: LAPACK library: ' + lapack.type_name())

# Combined dependency for convenience
blas_lapack_dep = declare_dependency(
  dependencies: [blas, lapack]
)

configure_file(
  output: 'xsf_config.h',
  configuration: conf_data
)

xsf_inc = include_directories('include')

xsf_dep = declare_dependency(
  include_directories: xsf_inc,
  link_with: xsf_lib,
  dependencies: blas_lapack_dep,  # Important: propagate BLAS/LAPACK
)

# Make it findable
meson.override_dependency('xsf', xsf_dep)


# ========================================
# Tests (optional)
# ========================================

if get_option('build_tests')
  subdir('tests')
endif

