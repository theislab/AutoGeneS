import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="autogenes",
    version="1.0",
    author="Hananeh Aliee, Maxim Schmidt",
    author_email="author@helmholtz-muenchen.de",
    description="Automatic Gene Selection",
    long_description=long_description,
    url="https://github.com/theislab/AutoGeneS",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
