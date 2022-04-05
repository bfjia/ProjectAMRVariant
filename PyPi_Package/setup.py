import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
    
# specify requirements of your package here
REQUIREMENTS = ['requests==2.27.1','biopython==1.72','urllib3==1.26.9','numpy==1.21.5','pandas==1.3.5']


setuptools.setup(
    name="STAPAMRTime",
    version="0.0.1",
    author="Baofeng Jia",
    author_email="Contact@bfjia.net",
    description="Metagenomic AMR variant detection tool, part of ARMTime.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/imasianxd/ProjectSTAPAMRTime",
    project_urls={
        "Bug Tracker": "https://github.com/imasianxd/ProjectSTAPAMRTime/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    install_requires=REQUIREMENTS,
    python_requires=">=3.7",
    entry_points = {
        'console_scripts': [
            'STAPAMRTime = STAPAMRTime.__main__:main'
        ]
    },

)
