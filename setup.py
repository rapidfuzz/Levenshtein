from __future__ import annotations

from skbuild import setup

with open("README.md", encoding="utf8") as f:
    readme = f.read()

setup(
    name="Levenshtein",
    version="0.25.1",
    url="https://github.com/rapidfuzz/Levenshtein",
    author="Max Bachmann",
    install_requires=["rapidfuzz >= 3.8.0, < 4.0.0"],
    author_email="pypi@maxbachmann.de",
    description="Python extension for computing string edit distances and similarities.",
    long_description=readme,
    long_description_content_type="text/markdown",
    license="GPL",
    license_file="COPYING",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
    ],
    packages=["Levenshtein"],
    package_dir={"": "src"},
    package_data={"Levenshtein": ["*.pyi", "py.typed"]},
    python_requires=">=3.8",
)
