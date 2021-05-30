import os
from setuptools import setup
from setuptools import find_packages

here = os.path.abspath(os.path.dirname(__file__))

version = {}
with open(os.path.join(here, "src", "__version__.py")) as f:
    exec(f.read(), version)


with open("README.md") as readme_file:
    readme = readme_file.read()

setup(
    name="satay",
    version=version["__version__"],
    description="A libray for processing sequencing data for SAturated Transposon Analysis in Yeast (SATAY)",
    long_description=readme,
    url="https://github.com/leilaicruz/LaanLab-SATAY-DataAnalysis",
    author="Liedewij Laan Lab",
    author_email="L.M.InigoDeLaCruz@tudelft.nl",
    license="Apache Software License 2.0",
    packages=find_packages(exclude=["*tests*"]),
    key_words=["transposon-mapping", "Saccharomyces Cerevisiae",],
    classifiers=[
        "Development Status :: 1 - Beta",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Bash",
    ],
    test_suite="tests",
    package_data={"satay": ["data_files/*"]},
)
