# Release Process (GitHub Releases + Optional PyPI)

This project publishes install artifacts (`.whl`, `.tar.gz`, `SHA256SUMS`) to GitHub Releases.

## 1. Bump version

Update `version` in `pyproject.toml`.

## 2. Tag and push

```bash
git add pyproject.toml
git commit -m "release: vX.Y.Z"
git tag vX.Y.Z
git push origin main --tags
```

## 3. Create GitHub Release

On GitHub:
1. Open `Releases` -> `Draft a new release`
2. Choose tag `vX.Y.Z`
3. Publish release

When published, `.github/workflows/python-publish.yml` will:
- Build wheel/sdist with `python -m build`
- Generate `SHA256SUMS`
- Upload all files from `dist/` to the release assets
- Publish to PyPI only if `PYPI_API_TOKEN` is configured

## 4. User install from release (no clone)

```bash
cd /content
curl -fsSL https://raw.githubusercontent.com/paulshamrat/ColabMDA/main/scripts/install_colabmda_release.sh -o install_colabmda_release.sh
bash install_colabmda_release.sh latest /content/colabmda
```
