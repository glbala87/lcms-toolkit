#!/usr/bin/env python
"""
Setup script for PyLCMS - Python LC-MS Data Analysis Toolkit.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
readme_path = Path(__file__).parent.parent / "README.md"
long_description = ""
if readme_path.exists():
    long_description = readme_path.read_text(encoding="utf-8")

setup(
    name="pylcms",
    version="1.0.0",
    author="LCMS Toolkit Contributors",
    author_email="",
    description="Python LC-MS Data Analysis Toolkit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lcms-toolkit/lcms-toolkit",
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    keywords=[
        "mass-spectrometry",
        "proteomics",
        "metabolomics",
        "lcms",
        "mzml",
        "mzxml",
        "bioinformatics",
    ],
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
    ],
    extras_require={
        "viz": [
            "matplotlib>=3.4.0",
        ],
        "dataframe": [
            "pandas>=1.3.0",
        ],
        "full": [
            "matplotlib>=3.4.0",
            "pandas>=1.3.0",
        ],
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=3.0.0",
            "black>=22.0.0",
            "isort>=5.10.0",
            "mypy>=0.950",
            "flake8>=4.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "pylcms=pylcms.cli:main",
        ],
    },
    project_urls={
        "Bug Reports": "https://github.com/lcms-toolkit/lcms-toolkit/issues",
        "Source": "https://github.com/lcms-toolkit/lcms-toolkit",
        "Documentation": "https://lcms-toolkit.readthedocs.io/",
    },
)
