from distutils.core import setup

setup(
    name="vbind",
    packages=["vbind"],
    version="1.0",
    license="BSD Clause 3",
    description="sRNA profiler",
    author="Pavithran Iyer and Charith Adkar",
    author_email="pavithran.iyer@uwaterloo.ca,charith.adkar@usherbrooke.ca",
    url="https://github.com/paviudes/vbind",
    download_url="https://github.com/paviudes/vbind/archive/refs/tags/v1.0.tar.gz",
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
        "Intended Audience :: Researchers and graduate students",
        "Topic :: Software Development :: Build Tools",
        "License :: OSI Approved :: BSD Clause 3",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
    ],
)
