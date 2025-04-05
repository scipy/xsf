# xsf
Special function implementations

## Tests

To run the tests:
- [clone this repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)
- `cd xsf`
- [install Pixi](https://pixi.sh/latest/#installation)
- `pixi run tests`

For subsequent test runs, you can skip re-cloning [`xsref`](https://github.com/scipy/xsref) with:

```shell
pixi run --skip-deps tests
```

You can trigger a rebuild inbetween test runs with:

```shell
pixi run build-tests
```

> [!NOTE]  
> This has currently only been tested on Linux.
