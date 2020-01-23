venv:
source ~/.virtualenvs/rocketPy/bin/activate

Test development:
install a local version of the package using
``pip install -e .``

To save the output of a example file to a text file:
``python examples/simple/simple.py  &> examples/simple/output.txt``




cookiecutter:
https://cookiecutter-pypackage.readthedocs.io/en/latest/tutorial.html


pytest and coverage:
1. Activate the virtual env
2. from root directory run
``python -m pytest -v``
and check that all the tests run
3. Check coverage by running
``coverage run -m --source=. pytest -v``
4. Report it by:
``coverage html``
or
``coverage report``
5. if using html, the report is created and available at
``htmlcov/index.html``
which is easily accessed by calling ``open htmlcov/index.html``
