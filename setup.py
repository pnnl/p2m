from setuptools import setup, find_packages

setup(
      name='p2m',
      version='0.1',
      description='Identify metabolites associated with a list of protein identifiers',
      author='Bryan J. Killinger',
      author_email='brykpnl@gmail.com',
      url='https://github.com/brykpnl/p2m',
      packages=find_packages(),
      entry_points={
          'console_scripts': ['p2m=p2m.main:main']
          },
      )
      
      
