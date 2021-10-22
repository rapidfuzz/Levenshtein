from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys

class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc', '-std=c++14', '/O2', '/W4'],
        'unix': ['-O3', '/std:c++14', '-Wextra', '-Wall', '-Wconversion', '-g0'],
    }
    l_opts = {
        'msvc': [],
        'unix': [],
    }

    if sys.platform == 'darwin':
        darwin_opts = ['-stdlib=libc++', '-mmacosx-version-min=10.9']
        c_opts['unix'] += darwin_opts
        l_opts['unix'] += darwin_opts

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        link_opts = self.l_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args += opts
            ext.extra_link_args += link_opts
        build_ext.build_extensions(self)

ext_modules = [
    Extension(
        name='Levenshtein._levenshtein',
        sources=[
            'src/Levenshtein-c/_levenshtein.cpp',
            'src/_levenshtein.cpp'
        ],
        include_dirs=[
            "src/Levenshtein-c/",
            "3rd-party/rapidfuzz-cpp/"
        ],
        language='c++'
    ),
    Extension(
        name='Levenshtein.c_levenshtein',
        sources=[
            'src/Levenshtein-c/_levenshtein.cpp',
            'src/c_levenshtein.cpp'
        ],
        include_dirs=[
            "src/Levenshtein-c/",
            "3rd-party/rapidfuzz-cpp/"
        ],
        language='c++'
    ),
]

if __name__ == "__main__":
    setup(
        cmdclass={'build_ext': BuildExt},
        ext_modules = ext_modules
    )
