# Shell-and-Tube HX Designer — segmented single-phase + steam heater

## What this version adds (requested)
- Input fields for shell-side detailed geometry (tube layout, pitch, clearances, sealing strips)
- Optional **auto-estimation** of Bell–Delaware leakage/bypass correction factors **Jl and Jb** from those inputs
  - Defaults are sensible starting points
  - You can always override J-factors manually

## Notes
- This is a practical engineering tool. The auto Jl/Jb is an **approximation** until full Bell–Delaware area bookkeeping
  (window vs crossflow tube counts, exact leakage streams, sealing strip effectiveness, etc.) is implemented.
- For glycol: concentration is **mass %**, entered via number boxes (keyboard + +/- buttons).
