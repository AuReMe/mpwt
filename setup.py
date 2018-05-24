import os

from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='mpwt',
      description='Multiprocessing for Pathway-Tools',
      long_description=readme(),
      version='0.1.6.4a1',
      url='https://gitlab.inria.fr/abelcour/mpwt',
      author='A. Belcour',
      author_email='arnaud.belcour@gmail.com',
      classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Audience
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',

        # Environnement, OS, languages
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 2'
      ],
      packages=['mpwt'],
      install_requires=[
            'biopython',
      ],
      python_requires='>=2.7,>=3.0.*',
      entry_points={
          'console_scripts': [
              'mpwt = mpwt.__main__:main'
          ]
      },
)