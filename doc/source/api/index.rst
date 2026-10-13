C++ API reference
=================

The API reference lists the functions declared in the ``xsf`` namespace,
grouped by header. Pages are generated from the headers under
``include/xsf`` at build time, so new functions appear automatically.
Implementation details in internal namespaces such as ``xsf::detail`` and
``xsf::cephes`` are omitted.

Most functions are not documented yet; for those, only the signature is
shown. Documentation is added incrementally through Doxygen comments in the
headers, and the header layout is expected to change as xsf is restructured.

General kernels
---------------

Headers in ``include/xsf``.

.. toctree::
   :maxdepth: 1
   :glob:

   generated/*

CPU-only kernels
----------------

Headers in ``include/xsf/cpu``, which contain CPU-only code.

.. toctree::
   :maxdepth: 1
   :glob:

   generated/cpu/*
