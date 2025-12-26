project(
  'xsf',
  'cpp',
  version : '0.1.5',
  license : 'BSD-3-Clause AND MIT AND BSD-3-Clause-LBNL AND Apache-2.0 WITH LLVM-exception',
  license_files : ['LICENSE', 'LICENSES_bundled.txt'],
  meson_version : '>=1.5.0',
  default_options : ['cpp_std=c++17']
)

# Get compiler
cpp = meson.get_compiler('cpp')

# ========================================
# BLAS/LAPACK Dependencies -- needed for Mathieu fcns.
# ========================================
blas_name = get_option('blas')
lapack_name = get_option('lapack')

if blas_name == ''
  # Auto-detect BLAS
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
  # Auto-detect LAPACK
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

# Combined dependency
blas_lapack_dep = declare_dependency(
  dependencies: [blas, lapack]
)

# ========================================
# Float128 Configuration
# ========================================

# Check for _Float128 support
float128_code = '''
int main() {
  _Float128 x = 1.0;
  return 0;
}
'''

has_float128 = cpp.compiles(float128_code, name: '_Float128 support')

# Configuration data
conf_data = configuration_data()
conf_data.set('HAVE_FLOAT128', has_float128)
configure_file(
  output: 'xsf_config.h',
  configuration: conf_data
)

# ========================================
# Includes
# ========================================
xsf_inc = include_directories('include')


# ========================================
# Declare Dependency
# ========================================

xsf_dep = declare_dependency(
  include_directories: xsf_inc,
  dependencies: blas_lapack_dep,
)

# Make it discoverable
meson.override_dependency('xsf', xsf_dep)

# ========================================
# Installation
# ========================================

if not meson.is_subproject()
  # Install headers
  install_headers(
    'include/xsf/mathieu/mathieu_fcns.h',
    subdir: 'xsf/mathieu'
  )
  
  # Install pkg-config file
  pkg = import('pkgconfig')
  pkg.generate(
    xsf_lib,
    name: 'xsf',
    description: 'XSF Special Functions Library',
    version: meson.project_version(),
  )
endif

# ========================================
# Tests (optional)
# ========================================

if get_option('build_tests')
  subdir('tests')
endif
