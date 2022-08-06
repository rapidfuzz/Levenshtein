from skbuild import setup

with open('README.md', 'rt', encoding="utf8") as f:
    readme = f.read()

setup(
    name="Levenshtein",
    version="0.20.2",
    url="https://github.com/maxbachmann/Levenshtein",
    author="Max Bachmann",
    install_requires=["rapidfuzz >= 2.3.0, < 3.0.0"],
    author_email="contact@maxbachmann.de",
    description="Python extension for computing string edit distances and similarities.",
    long_description=readme,
    long_description_content_type="text/markdown",

    license="GPL",
    license_file = "COPYING",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)"
    ],

    packages=["Levenshtein"],
    package_dir={'':'src'},
    python_requires=">=3.6"
)
