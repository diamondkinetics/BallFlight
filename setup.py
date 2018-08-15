from setuptools import setup
# also saw some to import distutils.core?

setup(name='ballflight',
      version='0.1.0',
      author='Devon Goetz',
      author_email='dkgoetz@mit.edu',
      packages=['ballflight', ],
      url='http://pypi.python.org/pypi/ballflight/',
      license='MIT',
      description='A package that allows you to model the flight of a ball from either a pitch or a swing.',
      long_description=open('README.md').read(),
      install_requires=[]
      )
