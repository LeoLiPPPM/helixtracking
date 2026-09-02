# TrackForge

TrackForge is a C++20 charged-particle track reconstruction pipeline. It
simulates a helical trajectory through a uniform solenoidal magnetic field,
adds detector noise and false hits, recovers the track with RANSAC and
least-squares refinement, and estimates transverse momentum from curvature.

## Build, test, and run

```bash
make -C trackforge test
make -C trackforge demo
python trackforge/analysis/plot_results.py --input trackforge/outputs
```

No third-party C++ library is required. The implementation uses small,
purpose-built linear algebra routines so the geometry is visible rather than
hidden behind a framework. For production detector work, those routines would
normally be replaced with a vetted package such as Eigen.

The demo uses a fixed random seed so its metrics and plots are reproducible.

