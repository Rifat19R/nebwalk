# PyPI Release Checklist

`nebwalk` is not published on PyPI yet. The repository is prepared for
trusted publishing through GitHub Actions.

## Local checks already expected before release

```bash
pytest tests/ -v
ruff check nebwalk/
python -m build
twine check dist/*
```

## Configure PyPI trusted publishing

Create a pending publisher on PyPI:

- PyPI project name: `nebwalk`
- GitHub owner: `Rifat19R`
- GitHub repository: `nebwalk`
- Workflow filename: `publish.yml`
- Environment name: `pypi`

PyPI supports pending publishers for projects that do not exist yet. The first
successful trusted-publishing run creates the project and turns the pending
publisher into a normal publisher.

## Publish v0.5.0

1. Confirm `pyproject.toml` and `nebwalk/__init__.py` both say `0.5.0`.
2. Push `main`.
3. Create GitHub release `v0.5.0`.
4. The `Publish` workflow builds the sdist/wheel, runs `twine check`, then
   publishes to PyPI using OIDC. No PyPI token should be committed or stored.

After publication:

```bash
python -m pip install nebwalk
python -c "import nebwalk; print(nebwalk.__version__)"
```
