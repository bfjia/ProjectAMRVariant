#!/bin/bash -e 
rm -rf dist
python -m build
python -m twine upload --repository pypi dist/* --verbose
source test/bin/activate
#pip install --no-cache-dir STAPAMRTime==$1
