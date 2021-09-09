from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys

class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc', '/O2', '/W4'],
        'unix': ['-O3', '-Wextra', '-Wall', '-Wconversion', '-g0'],
    }
    l_opts = {
        'msvc': [],
        'unix': [],
    }

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
            'src/Levenshtein-c/_levenshtein.c',
            'src/_levenshtein.c'
        ],
        include_dirs=[
            "src/Levenshtein-c/",
        ]
    ),
    Extension(
        name='Levenshtein.c_levenshtein',
        sources=[
            'src/Levenshtein-c/_levenshtein.c',
            'src/c_levenshtein.c'
        ],
        include_dirs=[
            "src/Levenshtein-c/",
        ]
    ),
]

if __name__ == "__main__":
    setup(
        cmdclass={'build_ext': BuildExt},
        ext_modules = ext_modules
    )
