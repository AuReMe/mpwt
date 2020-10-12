import os

from setuptools import setup

setup_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(setup_directory, 'README.rst'), encoding='utf-8') as readme_file:
    long_description = readme_file.read()

setup(name='mpwt',
      description='Multiprocessing for Pathway Tools',
      long_description=long_description,
      version='0.5.8',
      url='https://github.com/AuReMe/mpwt',
      author='Arnaud Belcour',
      python_requires='>=3.6',
      classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Audience
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',

        # Environnement, OS, languages
        'Programming Language :: Python :: 3'
      ],
      packages=['mpwt'],
      install_requires=[
            'biopython>=1.70',
            'chardet>=3.0.4',
            'docopt>=0.6.2',
            'gffutils>=0.9',
      ],
      entry_points={
          'console_scripts': [
              'mpwt = mpwt.__main__:run_mpwt'
          ]
      },
)
