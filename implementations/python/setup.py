"""Setup project."""

from setuptools import setup, find_packages

setup(
    name='mzlib',
    packages=find_packages(),
    version='0.1.0-alpha',
    description='HUPO-PSI Spectral library format',
    entry_points={
        'console_scripts': [
            "mzspeclib = mzlib.tools.cli:main"
        ],
    },
    classifiers=[
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 3 - Alpha"
    ],
    install_requires=[
        "sqlalchemy",
        "click",
        "psims >= 0.1.41",
        "pyteomics >= 4.5.3",
    ],
    test_requires=[
        "jsonschema"
    ]
)
