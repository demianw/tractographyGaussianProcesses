from numpy.distutils.misc_util import Configuration
from numpy.distutils.system_info import get_info
import os, sys

        
config = Configuration('tractographyGP',parent_package=None,top_path=None)

config.set_options(delegate_options_to_subpackages=True)

if len(sys.argv)==1:
  sys.argv.append('-h')

config.add_extension(name='bounded_thinPlate',
  sources=['tractographyGP/_bounded_thinPlate.f'],
  f2py_options=['no-lower'],
)

config_dict = config.todict()
try:
    config_dict.pop('packages')
except:
    pass



if __name__ == '__main__':
    from numpy.distutils.core import setup,Extension
    setup(  version="0.0",
            description="GP",
            author="Demian Wassermann",
            author_email="demian.wassermann@sophia.inria.fr",
            url="None",
            #download_url="",
            license="None yet",
            classifiers=[
                'Programming Language :: Python',
                'Programming Language :: Fortran',
                'Topic :: Scientific/Engineering',
                 ],
            requires=['NumPy (>=1.2)','pymc (==2.1)'],
            long_description="""
""",
            packages=["tractographyGP"],
            scripts=[
              'scripts/generateAggloClustering',
              'scripts/tractQuerier',
              'scripts/getPrototypeTract',
              ],
            **(config_dict))


