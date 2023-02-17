from setuptools import setup, find_packages

with open("requirements.txt") as f:
    required = [x for x in f.read().splitlines() if not x.startswith("#")]

setup(
    name="p2m",
    version="1.0.0",
    description="Identify metabolites associated with a list of protein identifiers",
    author="Pacific Northwest National Laboratory",
    author_email="brykpnl@gmail.com",
    url="https://github.com/pnnl/p2m",
    packages=find_packages(exclude=("test")),
    install_requires=required,
    entry_points={"console_scripts": ["p2m=p2m.main:main"]},
)
