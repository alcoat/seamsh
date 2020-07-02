import setuptools 
import sys
import os

#with open("README.md", "r") as fh:
#    long_description = fh.read()

long_description =''

version = "0.0.1"
commit_tag = os.environ.get("CI_COMMIT_TAG")
if commit_tag and (commit_tag.startswith("v-") or commit_tag.startswith("w-")):
    version = commit_tag[2:]




setuptools.setup(
    name="msea",
    version=version,
    author="Jonathan Lambrechts",
    author_email="jonathan.lambrechts@uclouvain.be",
    description="Ocean mesh generation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    url="https://www.migflow.be",
    packages=["msea"],
    package_dir={"msea":"msea"},
    #ext_modules=[CMakeExtension("msealib")],
    #ext_modules=[setuptools.Extension("seameshlib",["seamesh.c"])],
    package_data={"msea":["*.so","*.dll","*.dll.a","*.dylib"]},
    classifiers=[
        "Environment :: Console",
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS :: MacOS X",
        "Programming Language :: C",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering"
        ],
    install_requires=["scipy","numpy","gdal","gmsh-dev"],
    python_requires='>=3.6'
)




