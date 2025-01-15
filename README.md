# aot-cg-simulations
Scripts to perform simulations using coarse-grained models of AOT.

To use this repository, install [cookiecutter](https://www.cookiecutter.io/).
Cookiecutter can be used to copy the templates for the bilayer and self-assembly
simulations (e.g. `cookiecutter self-assembly-template`). It will prompt you for
the topology name, which is used to format the topology files, and a 'friendly
name', which is used in the job scripts.

To analyse the properties of bilayers, use the included `bilayer-analysis.py`
script. For properties of the self-assembly simulations (dilute isotropic
phase), see <https://github.com/a-ws-m/aot-analysis>.
