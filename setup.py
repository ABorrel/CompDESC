import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="CompDesc", # Replace with your own username
    version="0.8.4",
    author="Alexandre Borrel",
    author_email="a.borrel@gmail.com",
    description="Compute molecular descriptors 1D, 2D and 3D from a SMILES string ligand",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ABorrel/CompDESC",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    test_suite='nose.collector',
    tests_require=['nose'],
)
