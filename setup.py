from distutils.core import setup

setup(
    name="vbind",
    packages=["vbind"],
    version="1.0.1",
    license="BSD Clause 3",
    description="sRNA profiler",
    long_description="vbind is a userfriedly tool for RNA sequencing. In particular, it can be used to compute and visualize the bindings between a pool of sRNA nuleotides and a genome sequence. Please visit https://paviudes.github.io/vbind/ for more information.",
    author="Pavithran Iyer and Charith Adkar",
    author_email="pavithran.iyer@uwaterloo.ca,charith.adkar@usherbrooke.ca",
    url="https://github.com/paviudes/vbind",
    download_url="https://github.com/paviudes/vbind/archive/refs/tags/v1.1.tar.gz",
    keywords=["viroids", "sequencing", "sRNA profiling"],
    install_requires=[
        "scipy",
        "numpy",
        "multiprocessing",
        "matplotlib",
        "tqdm",
        "datetime",
    ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Build Tools",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
    ],
)
