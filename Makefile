# HEC-RAS 2D on Linux — common tasks.
#   make test            preprocessing tests, then every engine on every example
#   make test-fast       preprocessing only (no engines, ~1 min)
#   make test-a          Workflow A only (the validated path)
#   make docker          build hecras-linux:6.6, :7.0.1 and :7.0.2
#   make clean
PY      ?= python3
VERSIONS = 6.6 7.0.1 7.0.2

.PHONY: test test-fast test-a docker clean

test:
	$(PY) -m pytest tests -rfEX

test-fast:
	$(PY) -m pytest tests -m "not engine"

test-a:
	$(PY) -m pytest tests -k "TestWorkflowA or workflow_a"

docker:
	for v in $(VERSIONS); do \
	  docker build --platform linux/amd64 --build-arg HECRAS_VERSION=$$v -t hecras-linux:$$v . || exit 1; \
	done

clean:
	rm -rf .pytest_cache tests/__pycache__ tools/__pycache__ __pycache__
