import os
import inspect
import subprocess
from setuptools import setup, find_packages

import numpy as np


def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        git_revision = out.strip().decode('ascii')
    except OSError:
        git_revision = "Unknown"

    return git_revision


def get_version_info(version, is_released):
    fullversion = version
    if not is_released:
        git_revision = git_version()
        fullversion += '.dev0+' + git_revision[:7]
    return fullversion


def write_version_py(version, is_released, filename='pyaero/version.py'):
    fullversion = get_version_info(version, is_released)
    with open("./pyaero/version.py", "wb") as f:
        f.write(('__version__ = "%s"\n' % fullversion).encode())
    return fullversion


# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    setupdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    return open(os.path.join(setupdir, fname)).read()


#_____________________________________________________________________________


CLASSIFIERS = """\

Intended Audience :: Science/Research
Intended Audience :: Developers
Intended Audience :: Education
Intended Audience :: End Users/Desktop
Topic :: Scientific/Engineering
Topic :: Education
Topic :: Software Development
Topic :: Software Development :: Libraries :: Python Modules
Operating System :: Microsoft :: Windows
Operating System :: Unix
Operating System :: POSIX :: BSD
Programming Language :: Python :: 3.8
Programming Language :: Python :: 3.9
Programming Language :: Python :: 3.10
Programming Language :: Python :: 3.11
License :: OSI Approved :: BSD License

"""

is_released = True
version = '0.1.2'

fullversion = write_version_py(version, is_released)

package_data = {
        'pyaero': ['README.md', 'LICENSE', 'version.py', '*.pxd'],
        '': ['tests/*.*'],
        }

s = setup(
    name = "pyaero-dlm",
    version = fullversion,
    author = "Higor Luis Silva",
    author_email = "higor@ufu.br",
    description = ("Doublet lattice method"),
    long_description = read('README.md'),
    long_description_content_type = 'text/markdown',
    license = "BSD",
    keywords = "aeroelasticity",
    url = "https://github.com/saullocastro/pyaero",
    package_data = package_data,
    classifiers = [_f for _f in CLASSIFIERS.split('\n') if _f],
    packages = find_packages(),
)
