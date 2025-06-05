"""Test xsf functions work on GPU through CuPy.

Beyond cupy and pytest, also requires xsref package to be installed.
Should be run with pytest-isolate or pytest-forked to isolate tests
in separate processes since memory corruption on GPU can cause
failures to occur in unrelated tests.

TODO:
Set this up to run through pixi, and to run in CI.
"""

import os
import numpy as np
import polars as pl
import pytest
import shutil
import tempfile

from glob import glob
from packaging.version import Version
from pathlib import Path

from xsref.float_tools import extended_relative_error
from xsref.tables import get_input_rows, get_output_rows, get_in_out_types


#------------------------------------------------------------------------------
# Check if a module is present to be used in tests
#
# Copied from
# https://github.com/scipy/scipy/blob/1cbfa1c894557041d9825d6754a2c48fc3bec484/scipy/special/_testutils.py
#------------------------------------------------------------------------------

class MissingModule:
    def __init__(self, name):
        self.name = name


def check_version(module, min_ver):
    if type(module) is MissingModule:
        return pytest.mark.skip(reason=f"{module.name} is not installed")
    return pytest.mark.skipif(
        Version(module.__version__) < Version(min_ver),
        reason=f"{module.__name__} version >= {min_ver} required"
    )

#------------------------------------------------------------------------------

try:
    import cupy  # type: ignore
except (ImportError, AttributeError):
    cupy = MissingModule('cupy')


@pytest.fixture(scope="function", autouse=True)
def manage_cupy_cache():
    # Temporarily change cupy kernel cache location so kernel cache will not be polluted
    # by these tests. Remove temporary cache in teardown.
    temp_cache_dir = tempfile.mkdtemp()
    original_cache_dir = os.environ.get('CUPY_CACHE_DIR', None)
    os.environ['CUPY_CACHE_DIR'] = temp_cache_dir

    yield

    if original_cache_dir is not None:
        os.environ['CUPY_CACHE_DIR'] = original_cache_dir
    else:
        del os.environ['CUPY_CACHE_DIR']
    shutil.rmtree(temp_cache_dir)


def _get_cols_helper(table_path, xp):
    # This is more complicated than need be due to an oversight in xsref
    # functions for getting tables. The table metadata for input, output,
    # and err tables has the input and output types, but does not have
    # the column types of the current table, and does not contain the info
    # of whether current table is an input, output, or err table. It also
    # provides separate functions to get input rows from an input table,
    # output rows from an output table, and provides no function to get
    # rows from an err table (though this last task is easy to do because
    # there are no complex types involved in an err table).
    # TODO: Update xsref tables to contain the table's column types
    # in the metadata, and provide a single function for getting rows
    # from tables as lists of tuples that uses this metadata.
    table_name = table_path.name.lower()
    if table_name.startswith("in_"):
        rows = get_input_rows(table_path)
    elif table_name.startswith("out_"):
        rows = get_output_rows(table_path)
    elif table_name.startswith("err_"):
        table = pl.read_parquet(table_path)
        rows = table.to_numpy()
    result = tuple(xp.asarray(col) for col in zip(*rows))
    if len(result) == 1:
        result = result[0]
    return result


def get_cols_as_cupy(table_path):
    return _get_cols_helper(table_path, cupy)


def get_cols_as_numpy(table_path):
    return _get_cols_helper(table_path, np)


HERE = Path(__file__)


def get_tables_for_func(func_name):
    tables_path = (
        HERE.parent.parent.resolve() / "xsref" / "tables" / "scipy_special_tests"
    )
    tables_path /= func_name
    input_tables = list(tables_path.glob("In_*.parquet"))
    output_tables = [
        path.parent / path.name.replace("In_", "Out_") for path in input_tables
    ]
    err_tables = []
    for path in input_tables:
        types = path.name.removesuffix(".parquet").replace("In_", "")
        name = f"Err_{types}_other.parquet"
        err_tables.append(path.parent / name)
    return list(zip(input_tables, output_tables, err_tables))


def get_preamble(header):
    header_path = (HERE.parent.parent / "include" / Path(header)).resolve()
    return f'#include "{header_path}"'


@pytest.mark.usefixtures("manage_cupy_cache")
@check_version(cupy, "13.0.0")
class TestCuPy:
    def _adjust_tol(self, tol, *, wiggle=16):
        return wiggle * np.maximum(tol, np.finfo(tol.dtype).eps)

    @pytest.mark.parametrize(
        "tables_paths", get_tables_for_func("beta")
    )
    def test_beta(self, tables_paths):
        input_path, output_path, tol_path = tables_paths
        beta = cupy._core.create_ufunc(
            "cupyx_scipy_beta",
            ("ll->d", "LL->d", "ee->d", "ff->f", "dd->d"),
            "out0 = out0_type(xsf::beta(in0, in1));",
            preamble=get_preamble("xsf/beta.h"),
        )

        a, b = get_cols_as_cupy(input_path)
        out = cupy.asnumpy(beta(a, b))

        desired = get_cols_as_numpy(output_path)
        rtol = get_cols_as_numpy(tol_path)
        assert np.all(
            extended_relative_error(out, desired) <= self._adjust_tol(rtol)
        )

    @pytest.mark.parametrize(
        "tables_paths", get_tables_for_func("binom")
    )
    def test_binom(self, tables_paths):
        input_path, output_path, tol_path = tables_paths
        binom = cupy._core.create_ufunc(
            "cupyx_scipy_binom",
            ("ff->f", "dd->d"),
            "out0 = out0_type(xsf::binom(in0, in1));",
            preamble=get_preamble("xsf/binom.h"),
        )

        n, k = get_cols_as_cupy(input_path)
        out = cupy.asnumpy(binom(n, k))

        desired = get_cols_as_numpy(output_path)
        rtol = get_cols_as_numpy(tol_path)
        assert np.all(
            extended_relative_error(out, desired) <= self._adjust_tol(rtol)
        )

    @pytest.mark.parametrize(
        "tables_paths", get_tables_for_func("digamma")
    )
    def test_digamma(self, tables_paths):
        input_path, output_path, tol_path = tables_paths
        digamma = cupy._core.create_ufunc(
            'cupyx_scipy_special_digamma',
            (
                ('l->d', 'out0 = xsf::digamma(double(in0))'),
                ('e->d', 'out0 = xsf::digamma(double(in0))'),
                'f->f',
                'd->d',
                'F->F',
                'D->D',
            ),
            'out0 = xsf::digamma(in0)',
            preamble=get_preamble("xsf/digamma.h")
        )

        x = get_cols_as_cupy(input_path)
        out = cupy.asnumpy(digamma(x))

        desired = get_cols_as_numpy(output_path)
        rtol = get_cols_as_numpy(tol_path)
        assert np.all(
            extended_relative_error(out, desired) <= self._adjust_tol(rtol)
        )


    @pytest.mark.parametrize(
        "tables_paths", get_tables_for_func("ellipkinc")
    )
    def test_ellipkinc(self, tables_paths):
        input_path, output_path, tol_path = tables_paths
        ellipkinc = cupy._core.create_ufunc(
            'cupyx_scipy_special_ellipkinc', ('ff->f', 'dd->d'),
            'out0 = xsf::cephes::ellik(in0, in1)',
            preamble=get_preamble("xsf/cephes/ellik.h"),
        )

        phi, m = get_cols_as_cupy(input_path)
        out = cupy.asnumpy(ellipkinc(phi, m))

        desired = get_cols_as_numpy(output_path)
        rtol = get_cols_as_numpy(tol_path)
        assert np.all(
            extended_relative_error(out, desired) <= self._adjust_tol(rtol)
        )

    @pytest.mark.parametrize(
        "tables_paths", get_tables_for_func("ellipeinc")
    )
    def test_ellipeinc(self, tables_paths):
        input_path, output_path, tol_path = tables_paths
        ellipeinc = cupy._core.create_ufunc(
            'cupyx_scipy_special_ellipeinc', ('ff->f', 'dd->d'),
            'out0 = xsf::cephes::ellie(in0, in1)',
            preamble=get_preamble("xsf/cephes/ellie.h"),
        )

        phi, m = get_cols_as_cupy(input_path)
        out = cupy.asnumpy(ellipeinc(phi, m))

        desired = get_cols_as_numpy(output_path)
        rtol = get_cols_as_numpy(tol_path)
        assert np.all(
            extended_relative_error(out, desired) <= self._adjust_tol(rtol)
        )

    @pytest.mark.parametrize(
        "tables_paths", get_tables_for_func("expn")
    )
    def test_expn(self, tables_paths):
        input_path, output_path, tol_path = tables_paths
        expn = cupy._core.create_ufunc(
            'cupyx_scipy_special_expn',
            ('ff->f', 'dd->d'),
            'out0 = xsf::cephes::expn(in0, in1)',
            preamble=get_preamble("xsf/cephes/expn.h"),
        )

        x, n = get_cols_as_cupy(input_path)
        out = cupy.asnumpy(expn(x, n))

        desired = get_cols_as_numpy(output_path)
        rtol = get_cols_as_numpy(tol_path)
        assert np.all(
            extended_relative_error(out, desired) <= self._adjust_tol(rtol)
        )

    @pytest.mark.parametrize(
        "tables_paths", get_tables_for_func("gdtrib")
    )
    @pytest.mark.xfail(reason="Requires cpp_std >= 17")
    def test_gdtrib(self, tables_paths):
        input_path, output_path, tol_path = tables_paths
        gdtrib = cupy._core.create_ufunc(
            'cupyx_scipy_special_gdtrib',
            ('fff->f', 'ddd->d'),
            'out0 = xsf::gdtrib(in0, in1, in2)',
            preamble=get_preamble("xsf/cdflib.h"),
        )

        a, p, x = get_cols_as_cupy(input_path)
        out = cupy.asnumpy(gdtrib(a, p, x))

        desired = get_cols_as_numpy(output_path)
        rtol = get_cols_as_numpy(tol_path)
        assert np.all(
            extended_relative_error(out, desired) <= self._adjust_tol(rtol)
        )

    @pytest.mark.parametrize(
        "tables_paths", get_tables_for_func("hyp2f1")
    )
    def test_hyp2f1(self, tables_paths):
        input_path, output_path, tol_path = tables_paths
        hyp2f1 = cupy._core.create_ufunc(
            'cupyx_scipy_special_hyp2f1',
            ('ffff->f', 'dddd->d', 'fffF->F', 'dddD->D'),
            'out0 = xsf::hyp2f1(in0, in1, in2, in3)',
            preamble=get_preamble("xsf/hyp2f1.h"),
        )

        a, b, c, z = get_cols_as_cupy(input_path)
        if not cupy.iscomplexobj(z):
            pytest.xfail(
                "Real valued hyp2f1 currently broken on GPU due to use of"
                " recursion."
            )
        out = cupy.asnumpy(hyp2f1(a, b, c, z))

        desired = get_cols_as_numpy(output_path)
        rtol = get_cols_as_numpy(tol_path)
        assert np.all(
            extended_relative_error(out, desired) <= self._adjust_tol(rtol)
        )

    @pytest.mark.parametrize(
        "tables_paths", get_tables_for_func("lambertw")
    )
    def test_lambertw(self, tables_paths):
        input_path, output_path, tol_path = tables_paths
        _lambertw_scalar = cupy._core.create_ufunc(
            "cupyx_scipy_lambertw_scalar",
            ("Dld->D", "Fif->f"),
            "out0 = xsf::lambertw(in0, in1, in2)",
            preamble=get_preamble("xsf/lambertw.h"),
        )

        # A parameter called tol, not to be confused with the rtol for assessing
        # accuracy.
        z, k, tol = get_cols_as_cupy(input_path)
        out = cupy.asnumpy(_lambertw_scalar(z, k, tol))

        desired = get_cols_as_numpy(output_path)
        rtol = get_cols_as_numpy(tol_path)
        assert np.all(
            extended_relative_error(out, desired) <= self._adjust_tol(rtol)
        )

    @pytest.mark.parametrize(
        "tables_paths", get_tables_for_func("sici")
    )
    def test_sici(self, tables_paths):
        input_path, output_path, tol_path = tables_paths
        sici = cupy._core.create_ufunc(
            'cupyx_scipy_special_sici',
            (
                (
                    'f->ff',
                    '''
                    float si, ci;
                    xsf::sici(in0, si, ci);
                    out0 = si; out1 = ci;
                    ''',
                ),
            (
                'd->dd',
                '''
                double si, ci;
                xsf::sici(in0, si, ci);
                out0 = si; out1 = ci;
                ''',
            ),
                (
                    'F->FF',
                    '''
                    complex<float> si, ci;
                    xsf::sici(in0, si, ci);
                out0 = si; out1 = ci;
                    ''',
                ),
                (
                    'D->DD',
                    '''
                    complex<double> si, ci;
                    xsf::sici(in0, si, ci);
                    out0 = si; out1 = ci;
                    ''',
                ),
            ),
            preamble=get_preamble("xsf/sici.h"),
        )
            
        x = get_cols_as_cupy(input_path)
        if cupy.iscomplexobj(x):
            pytest.xfail("Known bug, returning nan instead of a complex infinity.")
        out0, out1 = map(cupy.asnumpy, sici(x))

        desired0, desired1 = get_cols_as_numpy(output_path)
        rtol0, rtol1 = get_cols_as_numpy(tol_path)
        error = extended_relative_error(out0, desired0)
        tol = self._adjust_tol(rtol0)
        assert np.all(
            extended_relative_error(out0, desired0) <= self._adjust_tol(rtol0)
        )
        assert np.all(
            extended_relative_error(out1, desired1) <= self._adjust_tol(rtol1)
        )

    @pytest.mark.parametrize(
        "tables_paths", get_tables_for_func("shichi")
    )
    def test_shichi(self, tables_paths):
        input_path, output_path, tol_path = tables_paths
        shichi = cupy._core.create_ufunc(
            'cupyx_scipy_special_shichi',
            (
                (
                    'f->ff',
                    '''
                    float shi, chi;
                    xsf::shichi(in0, shi, chi);
                    out0 = shi; out1 = chi;
                    ''',
                ),
                (
                    'd->dd',
                    '''
                    double shi, chi;
                    xsf::shichi(in0, shi, chi);
                    out0 = shi; out1 = chi;
                    ''',
                ),
                (
                    'F->FF',
                    '''
                    complex<float> shi, chi;
                    xsf::shichi(in0, shi, chi);
                    out0 = shi; out1 = chi;
                    ''',
                ),
                (
                    'D->DD',
                    '''
                    complex<double> shi, chi;
                    xsf::shichi(in0, shi, chi);
                    out0 = shi; out1 = chi;
                    ''',
                ),
            ),
            preamble=get_preamble("xsf/sici.h"),
        )
            
        x = get_cols_as_cupy(input_path)
        if cupy.iscomplexobj(x):
            pytest.xfail("Known bug, returning nan instead of a complex infinity.")
        out0, out1 = map(cupy.asnumpy, shichi(x))

        desired0, desired1 = get_cols_as_numpy(output_path)
        rtol0, rtol1 = get_cols_as_numpy(tol_path)
        assert np.all(
            extended_relative_error(out0, desired0) <= self._adjust_tol(rtol0)
        )
        assert np.all(
            extended_relative_error(out1, desired1) <= self._adjust_tol(rtol1)
        )
