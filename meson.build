project('CTMF', 'cpp',
  default_options : ['buildtype=release', 'b_ndebug=if-release', 'cpp_std=c++14'],
  meson_version : '>=0.48.0',
  version : '4'
)

sources = [
  'CTMF/CTMF.cpp',
  'CTMF/CTMF.hpp',
  'CTMF/vectorclass/instrset.h',
  'CTMF/vectorclass/instrset_detect.cpp'
]

libs = []

vapoursynth_dep = dependency('vapoursynth').partial_dependency(compile_args : true, includes : true)

add_project_arguments('-ffast-math', language : 'cpp')

if host_machine.cpu_family().startswith('x86')
  sources += [
    'CTMF/CTMF_SSE2.cpp',
    'CTMF/vectorclass/vectorclass.h',
    'CTMF/vectorclass/vectorf128.h',
    'CTMF/vectorclass/vectorf256.h',
    'CTMF/vectorclass/vectorf256e.h',
    'CTMF/vectorclass/vectori128.h',
    'CTMF/vectorclass/vectori256.h',
    'CTMF/vectorclass/vectori256e.h'
  ]

  add_project_arguments('-DVS_TARGET_CPU_X86', '-mfpmath=sse', '-msse2', language : 'cpp')

  libs += static_library('avx2', 'CTMF/CTMF_AVX2.cpp',
    dependencies : vapoursynth_dep,
    cpp_args : ['-mavx2', '-mfma'],
    gnu_symbol_visibility : 'hidden'
  )
endif

shared_module('ctmf', sources,
  dependencies : vapoursynth_dep,
  link_with : libs,
  install : true,
  install_dir : join_paths(vapoursynth_dep.get_pkgconfig_variable('libdir'), 'vapoursynth'),
  gnu_symbol_visibility : 'hidden'
)
