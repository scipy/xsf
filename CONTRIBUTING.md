# Contributing to xsf

Development of xsf is made easy with [Pixi](https://pixi.prefix.dev/).
As a first step, [install Pixi](https://pixi.prefix.dev/latest/installation).

## Clone the repository

Clone the xsf repository, following
<https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository>.
In summary, fork the repository on GitHub,
then set `origin` to your fork and add `upstream`:

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

## Development tasks

All development tasks are then available via `pixi run`:

```bash
pixi run tests     # run the tests
pixi run coverage  # generate test coverage report
pixi run format    # format source files
```

Run `pixi task list` for a full list of available tasks.

Multiple environments are available for the `tests` task:

```bash
pixi run --environment=tests-release tests  # test with release build type
pixi run --environment=tests-debug tests    # test with debug build type
pixi run --environment=coverage tests       # test with coverage build type
```

Run specific tests by passing `ctest` flags to `pixi run tests`,
e.g. `pixi run tests -R xsf_tests_test_boxcox`,
following <https://cmake.org/cmake/help/latest/manual/ctest.1.html#cmdoption-ctest-R>.

Run `pixi info` for a full list of environments and their tasks.

## Pull requests

Run the relevant tests and formatting checks before opening a pull request.
