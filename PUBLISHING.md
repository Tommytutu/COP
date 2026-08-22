# Publishing to PyPI

The repository builds the distributions
`ahpcop-1.0.0-py3-none-any.whl` and
`ahpcop-1.0.0.tar.gz` from `pyproject.toml`.

## Recommended: PyPI trusted publishing

1. Sign in to PyPI and create a pending trusted publisher for the project
   `AHPCOP`.
2. Use owner `Tommytutu`, repository `COP`, workflow
   `.github/workflows/publish.yml`, and environment `pypi`.
3. In GitHub, create the `pypi` environment for this repository.
4. Create and publish a GitHub release such as `v1.0.0`.

The release workflow builds both distributions and publishes them through
OpenID Connect, so no long-lived PyPI API token is stored in GitHub.

## Local validation

```powershell
python -m pip install build twine
python -m build
python -m twine check dist\*
```

To install the built wheel before publishing:

```powershell
python -m pip install dist\ahpcop-1.0.0-py3-none-any.whl
```

After publishing, users install the release with:

```powershell
python -m pip install AHPCOP
```
