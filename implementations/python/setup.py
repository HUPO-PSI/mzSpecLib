"""Setup project."""

from setuptools import setup

setup(
    name='mzlib',
    version='0.1.0-alpha',
    description='HUPO-PSI Spectral library format',
    packages=['mzlib'],
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
    ],
    test_requires=[
        "jsonschema"
    ]
)
