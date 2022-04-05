# eoas-pyutils
A set of utilities for EOAS including IO, Visualization, Preprocessing, Metrics, etc. 

# Usage

The idea of this repo is to use it as a submodule. 
To include it in your project you need to:

1. Add it as a submodule together with *hycom-utils*
```shell
git submodule add git@github.com:HYCOM/HYCOM-utilities.git hycom_utils
git submodule add git@github.com:olmozavala/eoas-pyutils.git eoas_pyutils
```
2. Recursively download all the files
```shell
git submodule update --init --recursive
```
3. Include the proper paths in your python files
```python
import sys
sys.path.append("eoas_pyutils/")
sys.path.append("eoas_pyutils/hycom_utils/python")
```
4. Depending on the IDE you are using you also need to include
the `eoas_pyutils` and the `eoas_pyutils/hycom_utils/python` folder as **source** folders.