from setuptools import setup, find_packages

setup(
    name="pyemfield",
    version="0.1.5",
    author="Lin, Ming Chih",
    author_email="mingchih.lin8@gmail.com",
    description="A Python library for HFSS antenna radiation pattern analysis and gain optimization",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/linmingchih/pyemfield",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "numpy>=1.18.0",
        "scipy>=1.4.0",
        "matplotlib>=3.1.0",
        "pyaedt>=0.10.1",  # 確保這個版本符合你的需求
    ],
)
