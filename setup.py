from setuptools import setup

setup(name='eccentricdark',
      version='0.0.0',
      description='Simulating eccentricity distributions of massive binary formation channels.',
      keywords=['cosmology black hole binary eccentric eccentricity gravitational wave LIGO VIRGO DECIGO LISA'], 
      url='https://github.com/ndeporzio/eccentricdark',
      author='Nicholas DePorzio',
      author_email='nicholasdeporzio@g.harvard.edu',
      license='MIT',
      packages=['eccentricdark'],
      include_package_data=True,
      zip_safe=False,
      test_suite='nose.collector',
      tests_require=['nose', 'coverage'],
      # dependency_links=[''],
      scripts=[
      #  '' 
      ], 
      install_requires=[
          'numpy',
          'matplotlib',
          'pandas',
          'seaborn', 
          'scipy',
          'dill'
      ])
