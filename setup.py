import os

from setuptools import setup

with open(os.path.join(os.path.dirname(__file__), 'README.rst'), encoding='utf-8') as readme_file:
      readme = readme_file.read()

setup(name='mpwt',
      description='Multiprocessing for Pathway-Tools',
      long_description=readme,
      version='0.1.4.3a1',
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
        'Programming Language :: Python :: 3'
      ],
      packages=['mpwt'],
      install_requires=[
            'biopython',
      ],
      entry_points={
          'console_scripts': [
              'mpwt = mpwt.__main__:main'
          ]
      },
)