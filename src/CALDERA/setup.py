import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="caldera",
    version="0.0.2",
    author="Hector Roux de Bezieux",
    author_email="hector.rouxdebezieux@berkeley.edu",
    description="CALDERA",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/HectoRDB/CALDERA",
    packages=setuptools.find_packages(),
    scripts=['bin/caldera-script',
             'bin/Pre-Process/toMajor.py',
             'bin/Pre-Process/cutNewick.R'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    zip_safe=False
)