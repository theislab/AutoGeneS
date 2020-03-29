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
    long_description_content_type="text/markdown",
    url="https://github.com/theislab/AutoGeneS",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
      'pandas>=0.25.1',
      'anndata>=0.6.22.post1',
      'numpy>=1.17.2'
      'dill>=0.3.1.1'
      'deap>=1.3.0'
      'scipy>=1.3.1'
      'cachetools>=3.1.1'
      'scikit-learn>=0.21.3'
      'matplotlib>=3.0.3'
    ]
)
