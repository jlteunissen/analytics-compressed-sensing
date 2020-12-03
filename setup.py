import json
from setuptools import setup, find_packages

with open('metainfo.json') as file:
    metainfo = json.load(file)

# # the tutorial was tested for the following packages
# pandas==0.20.3
# matplotlib==2.1.1
# scikit-learn==0.19.1
# bokeh==1.0.4
# ase==3.17.0
# numpy==1.13.3

setup(
    name='compressed_sensing',
    version='1.0',
    author=', '.join(metainfo['authors']),
    author_email=metainfo['email'],
    url=metainfo['url'],
    description=metainfo['title'],
    long_description=metainfo['description'],
    packages=find_packages(),
    install_requires=['pandas', 'matplotlib', 'scikit-learn', 'bokeh', 'ase', 'numpy', 'colorcet', 'jupyter_jsmol'],
)


