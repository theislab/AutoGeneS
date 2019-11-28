import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="autogenes",
    version="0.9.1",
    author="Hananeh Aliee, Maxim Schmidt",
    author_email="author@example.com",
    description="Automatic Gene Selection",
    long_description=long_description,
    url="https://github.com/lila167/AutoGeneS",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
