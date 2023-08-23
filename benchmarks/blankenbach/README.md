# The Blankenbach convection benchmark

This folder allows running the benchmark cases 1a, 1b, 1c, 2a, and 2b from {cite:t}`BBC89`:

Blankenbach, B., et al. "A benchmark comparison for mantle convection codes."
Geophysical Journal International 98.1 (1989): 23-38.

The reference values (given in [reference_statistics.txt](reference_statistics.txt)) are
from that paper (see Table 9). After building the plugin in the `plugin/` subfolder
and running the benchmarks by executing the script `run_all_models.sh`, you can
use the python script `get_statistics.py` to analyze and plot the results.
Please note that these computations take a long time to reach steady state,
especially on finer meshes, consider adjusting the number of parallel
processes and the choice of refinement levels in `run_all_models.sh`
and use ASPECT's optimized mode
(see {ref}`sec:run-aspect:run-faster:debug-vs-optimized`).

The reference output for this benchmark is stored in [statistics.txt](statistics.txt)
and shown in {numref}`fig:blankenbach`.

```{figure-md} fig:blankenbach
<img src="blankenbach.*" alt="Screenshot"  width="100%"/>

Convergence of the Blankenbach benchmark cases. Left panel: Nusselt number over
cell size h. Right panel: Relative error over cell size. Gray lines
indicate the reference values from {cite:t}`BBC89` (left panel) and
theoretical convergence predicted for different convergence orders
(right panel). Dashed lines show a simple
heat flux computed as the gradient of temperature, solid lines show the
heat flux computed using a consistent boundary flux method.
```

```{table} For reference the number of timesteps needed to compute these results:
:name: tab:blankenbach
| Case | Resolution | #Timesteps |
| :--- | ---------: | ---------: |
|   1a | 3          | 72         |
|   1a | 4          | 128        |
|   1a | 5          | 232        |
|   1a | 6          | 470        |
|   1a | 7          | 939        |
|   1b | 3          | 210        |
|   1b | 4          | 347        |
|   1b | 5          | 690        |
|   1b | 6          | 1368       |
|   1b | 7          | 2745       |
|   1c | 3          | 830        |
|   1c | 4          | 1916       |
|   1c | 5          | 3397       |
|   1c | 6          | 6728       |
|   1c | 7          | 13445      |
|   2a | 3          | 1857       |
|   2a | 4          | 3088       |
|   2a | 5          | 6638       |
|   2a | 6          | 13593      |
|   2a | 7          | 27244      |
|   2b | 3          | 4076       |
|   2b | 4          | 8404       |
|   2b | 5          | 17034      |
|   2b | 6          | 34164      |
|   2b | 7          | 68363      |
```
