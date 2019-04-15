import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="aspect_postprocess_utils",
    version="0.0.1",
    author="The authors of the ASPECT code",
    author_email="author@example.com",
    description="Utilites to help with some types of postprocessing for ASPECT",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/geodynamics/aspect",
    packages=setuptools.find_packages(),
    classifiers=(
        "Classifier :: Invalid"
    ),
    install_requires=[
        "numpy",
    ],
    entry_points = {
        'console_scripts' : [
            "aspect-postprocess=aspect_postprocess_utils.cli:base_cli",
        ],
    }
)
