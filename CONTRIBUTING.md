# Contributing to xsf

Install [Pixi](https://pixi.sh/latest/#installation) first. Pixi provisions the compiler toolchain and other dependencies defined in `pixi.toml`.

> [!NOTE]
> The documented commands below have only been tested on macOS and Linux.

## Clone the repository

Fork the repository, then set `origin` to your fork and add `upstream`:

```bash
git clone https://github.com/your-github-id/xsf.git
cd xsf
git remote add upstream https://github.com/scipy/xsf.git
```

Check remotes with:

```bash
git remote -v
```

You should see:

```bash
origin  https://github.com/your-github-id/xsf.git (fetch)
origin  https://github.com/your-github-id/xsf.git (push)
upstream    https://github.com/scipy/xsf.git (fetch)
upstream    https://github.com/scipy/xsf.git (push)
```

To sync your branch with upstream `main`:

```bash
git fetch upstream
git checkout main
git merge upstream/main
```

## Run tests

Run the default C++ test workflow:

```shell
pixi run tests
```

For incremental local work, you can run the steps separately:

```shell
pixi run clone-xsref
pixi run configure-tests
pixi run build-only -j8
pixi run --skip-deps tests -j2
```

Useful variants:

- Rebuild after changing code: `pixi run build-tests`
- Debug test build with extra assertions: `pixi run tests-debug`

## Formatting

Check formatting without modifying files:

```shell
pixi run format-dry-error
```

Apply formatting:

```shell
pixi run format
```

## Coverage

Coverage support is wired through the `tests-ci` environment and CMake coverage
targets.

Generate a local HTML coverage report with:

```shell
pixi run --environment=tests-ci clone-xsref
pixi run --environment=tests-ci configure-coverage
pixi run --skip-deps --environment=tests-ci build-tests-ci
pixi run --skip-deps --environment=tests-ci tests-ci
pixi run --skip-deps --environment=tests-ci coverage
```

The generated report is written to `build/coverage_report/index.html`.

For subsequent runs in the same build tree, it is usually enough to rerun:

```shell
pixi run --skip-deps --environment=tests-ci tests-ci
pixi run --skip-deps --environment=tests-ci coverage
```

## Pull requests

Run the relevant tests and formatting checks before opening a pull request.
