# coding utf8
import setuptools
from hugep2g.versions import get_versions

with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

setuptools.setup(
    name="HugeP2G",
    version=get_versions(),
    author="Yuxing Xu",
    author_email="xuyuxing@mail.kib.ac.cn",
    description="Aligning a large number of protein sequences to a genome (use genblasta and genewise)",
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url="https://github.com/SouthernCD/HugeP2G",
    include_package_data = True,

    entry_points={
        "console_scripts": ["HugeP2G = hugep2g.cli:main"]
    },    

    # package_data={
    #     "hugep2g": ['dep/*'],
    # },

    packages=setuptools.find_packages(),

    install_requires=[
        "toolbiox>=0.0.12",
        "bcbio-gff>=0.6.6",
        "biopython>=1.76",
        "interlap>=0.2.6",
        # "pandas>=1.0.1",
        "pyfaidx>=0.5.5.2",
    ],

    python_requires='>=3.5',
)
