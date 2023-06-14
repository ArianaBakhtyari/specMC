import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt") as f:
    requirements = [x.strip() for x in f.read().splitlines() if x.strip()]

setuptools.setup(
    name="", 
    version="1.1",
    author="Ariana Bakhtyari",
    author_email="ariana.bakhtyari@queensu.ca",
    description="Bayesian Spectral Model",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ArianaBakhtyari/specMC",
    packages=setuptools.find_packages(),
    classifiers=[],
    python_requires='>=3.6',
)