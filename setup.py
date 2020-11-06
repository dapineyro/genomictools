import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="genomictools",
    version="0.0.1",
    author="David Piñeyro",
    author_email="dapineyro.dev@gmail.com",
    license="GPL-3",
    description="Genomic Tools to handle genes and sequences",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dapineyro/genomictools",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires='>=3.8',
)
