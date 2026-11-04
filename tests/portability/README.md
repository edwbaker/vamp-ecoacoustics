# Portability checks

`.github/workflows/portability.yml` checks a release's binaries in real Vamp hosts on Linux, macOS and Windows:

- the native library (`.so`, `.dylib` or `.dll`) through [ReVAMP](https://cran.r-project.org/package=ReVAMP), installed
  from CRAN, with `check_release.R`;
- the WebAssembly library through [vamp-wasm](https://github.com/edwbaker/vamp-wasm) in Node.js, with `check_wasm.js`.

Both run every plugin at its defaults on the 91 pilot recordings of the Acoustic Index User's Guide
([Zenodo](https://doi.org/10.5281/zenodo.11004284)), checked against their md5, and compare each value with the R
function's value in the 0.3.0 benchmark ([NHM Data Portal](https://doi.org/10.5519/ijg223e8)), held in
`reference.csv`. A value passes within 1e-4 of the R value, relative to it, plus 1e-6 absolute, as soundecology rounds
some indices to six decimal places. The results also give each value's difference from the plugin's value in the
benchmark, where the plugins were built from source and run through ReVAMP on Linux.

The workflow runs after each release (from `build.yml`), on demand for any release (Actions, Release portability, Run
workflow) and when these files change. Each platform's results are uploaded as an artifact, and `summarise.py` combines
them into `portability-<tag>.csv`.

To run the checks locally, with `VAMP_PATH` naming the folder that holds the native library:

    Rscript tests/portability/check_release.R <recordings folder> results.csv
    node tests/portability/check_wasm.js <vamp-wasm folder> vamp-ecoacoustics.wasm <recordings folder> results.csv

`reference.csv` was made with `make_reference.py` from the dataset's `aiguide.csv` and the recordings.
