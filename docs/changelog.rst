Changelog
=========

0.2.1
-----

* frcmod DIHE lines carry ``SCEE=``/``SCNB=`` only when a term's 1-4 scaling is not Amber's
  default (1.2 / 2.0). leap treats explicit values specially, so default terms now look like
  parmchk2's output.

0.2.0
-----

* New Python API: :func:`unleap.extract_ligand` and :func:`unleap.extract_ligand_from_structure`
  return a :class:`~unleap.LigandParams`. It holds every output as a string, the structured
  validation results (:class:`~unleap.CheckResult`), warnings, and JSON-able metadata
  (``to_dict``). See :doc:`user/python_api`.
* The library no longer prints or exits. It logs to the ``unleap`` logger and raises
  :class:`~unleap.UnleapError`.
* Validation scratch files go to a chosen ``workdir`` or a temporary directory. The CLI keeps
  using ``OUTDIR/check``.
* The ligand is stripped from a copy of the loaded topology instead of re-reading the parm7.
* ``import unleap`` no longer imports ParmEd or RDKit until the API is used.
* CLI output files are unchanged. In ``report.txt`` the ``wrote ...`` lines now come last.

0.1.0
-----

* First release: the ``unleap`` command line tool.
