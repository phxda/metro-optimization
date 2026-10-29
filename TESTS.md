# TESTS

Everything here was run before the sweep (unit tests, smoke test, visual inspection) or after it (audit of all 216 results).

## Unit tests with hand-computed answers (`tests/test_model.py`, 6 tests, all pass)
- **Cost.** A 3-station line with hops of 3 km and 5 km, with one water column on the second hop: cost must equal 2.5/km x 8 km + 1.5 x 3 stations + 4.0 x (water km sampled on the hop). A second line sharing a station must charge the shared station once and track once per line-hop.
- **Travel time.** One line, 60 km/h, 1 min dwell, fleet 14: headway hits the 2 min floor, so wait term = 2.0 x 2/2 = 2. Station-to-station generalized time 0->3 is 2 + 4 = 6, 0->8 is 2 + 4 + 6 = 12, reverse direction equal. With fleet 4 the headway is 28/4 = 7 min and the 0->3 time becomes 7 + 4 = 11.
- **Transfers.** Two crossing lines sharing a hub. Transfer counts from the shortest-path predecessor tree are 0 for same-line pairs and 1 for cross-line pairs. Times: same-line 5, cross-line 5 + 10 (penalty) + 5 = 20. A fleet too small for the line count (V_MIN = 2 per line) returns infeasible.
- **Fitness end to end.** A 2-station line on a 1x10 grid with all other zones outside the access radius. Fitness is compared with a closed-form value built from the drive time, access walk, wait, ride, logit share and gravity flow. Agreement to 1e-4.
- **Real-network encoding.** After snapping, no stylized Baltimore station is in water and no line repeats a station.
- **Penalty invariance.** A single straight line has no transfers, so its fitness does not depend on the transfer penalty.

## Optimizer invariants (`tests/test_evo.py`, 3 tests, all pass)
- Phase freezing: after a phased optimization step, every frozen station sequence is still a contiguous, unchanged slice of a line in the new network (extension at line ends is the only allowed change).
- 400 random mutation chains on the river city never produce a station in water, an invalid hop, a duplicate station, or a cost above the cap.
- A phased run at the low budget ends within budget and records three phase fitnesses.

## Smoke test and visual inspection
Planned runs (300-400 generations) on all six cities were rendered and inspected. This found three problems that were fixed before the sweep:
1. Transit never beat driving, so the optimizer's best network was empty. Fixed by adding congestion and parking pressure to the driving baseline, and by moving to the standard logit with a time-difference slope (see PAPER.md, Methods). The fitness definition was also changed so negative savings among logit riders are bounded.
2. Access walking was over-penalized by the 2 km zone size (zone offset too large), which made single serpentine lines the optimum. Fixed by reducing the zone offset and access weight.
3. The Triangle optimizer converged to a corner-of-the-map 2-station network because every network had fitness <= 0. Fixed by recalibrating access (feeder/park-and-ride speed) and driving speed, and adding a tiny cost tie-break.
A second inspection of the final grid (shape x budget) and of the Baltimore, Triangle and planned-vs-phased figures showed no lines through water except costed crossings, no single line covering the whole city except at the lowest budget, and no budget violations.

## Audit of all 216 results (`tests/check_runs.py`)
0 violations: every network is within its budget, has all stations on land, no duplicate stations, all hops within the min/max hop length, and a feasible fleet. Warnings (not errors): 105 runs are single-line networks (low budgets, and the Triangle at all budgets) and 8 have one line holding >80% of track. These are findings, not bugs: at the low budget the fixed fleet and the layover overhead favor one line.

Not tested: optimality (the optimizer is a heuristic; seed-to-seed coefficient of variation of planned fitness is 5% median and 24% max), and no calibration against observed ridership.
