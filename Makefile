.PHONY: clean format build docs test lint

clean:
	@rm -rf build dist .eggs *.egg-info
	@rm -rf .benchmarks .coverage coverage.xml htmlcov report.xml
	@rm -rf .mypy_cache
	@rm -rf .ipynb_checkpoints
	@rm -rf docs/build docs/generated
	@find . -type d -name '__pycache__' -exec rm -rf {} +
	@find . -type d -name '*pytest_cache*' -exec rm -rf {} +
	@find . -type f -name "*.py[co]" -exec rm -rf {} +

format: clean
	@black wind_stats tests
	@isort wind_stats tests

build:
	@python -m build

install:
	@python -m pip install -e .[test]

publish:
	@poetry publish

docs:
	@rm -rf docs/build docs/generated
	@sphinx-build docs docs/build

test:
	@pytest --cov=wind_stats --cov-report=term-missing --cov-report=xml tests

lint:
	@mypy wind_stats
	@ruff check wind_stats tests
	@black wind_stats tests --check
	@isort wind_stats tests --check-only