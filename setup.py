import setuptools
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='dspyields',
    version='0.0.1',
    packages=setuptools.find_packages(),
    license='MIT',
    author='Jia Xue',
    author_email='jiaxue@tju.edu.cn',
    description='Calculation of Maximum Theoretical Yield',
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=['path', 'requests', 'rdkit', 'cobra'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",    ],
    python_requires='>=3.7',
    include_package_data=True,
)