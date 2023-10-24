build:
	python3 setup.py build_ext --inplace
	python3 setup_wheel.py bdist_wheel

install:
	python3 -m pip install --ignore-installed dist/*.whl

sudo-install:
	sudo python3 -m pip install --ignore-installed dist/*.whl

clean:
	rm -rf build
	rm -rf dist
	rm -f libpointing/*.so
