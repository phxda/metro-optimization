# PREDICTIONS (written before any code; not to be edited)

Budgets are low / mid / high (roughly 1x, 2x, 4x of a base track-km allowance). "Transfer penalty" is low vs. high. Metrics are measured on end-state networks, mean over 3 seeds.

**Q1. Topology vs. budget and city shape**
- P1a: Line count rises with budget in every city, but sub-linearly: high-budget networks have between 1.5x and 3x the lines of low-budget networks (not 4x).
- P1b: Monocentric city produces a radial/star network with 0 closed rings at low and mid budgets; a ring (cycle in the station graph) appears in at least 2 of 3 seeds only at the high budget, if at all.
- P1c: Linear coastal city produces nearly collinear lines: the dominant line is within 20 degrees of the coast axis, and the network has at most 1 cycle at every budget.
- P1d: River-split city: crossings are concentrated; the number of distinct river crossings is at most 2 at low budget and at most 4 at high budget, and exceeds 1 line per crossing only at high budget.
- P1e: Polycentric city yields the highest transfer count (transfer stations per line) of the four synthetic shapes at mid and high budgets.
- P1f: Fraction of total demand captured by transit (mode share) rises with budget, with diminishing returns: the low-to-mid gain exceeds the mid-to-high gain.

**Q2. Transfer penalty**
- P2a: Mean line length (stations per line) is at least 15% greater under the high transfer penalty than the low one, pooled over cities and budgets.
- P2b: Line count is lower under the high penalty in at least 4 of 6 cities at the mid budget.
- P2c: Number of transfer stations is lower under the high penalty, with a typical reduction of 10 to 30%.
- P2d: Total rider time saved is lower under the high penalty (it is a harsher cost), but the optimizer-recovered fraction is within 10% of the low-penalty fitness after normalizing for the penalty.

**Q3. Phased vs. planned**
- P3a: Phased networks have lower fitness than planned in the large majority of conditions; the median gap is 3 to 10% of planned fitness, with a worst case below 25%.
- P3b: The gap is larger at the high budget than at the low budget, since early lock-in matters more with more tranches to place.
- P3c: Phased networks have more lines with fewer stations per line (more fragmented), and a lower share of the longest line, and fewer rings than planned networks.
- P3d: The phased gap is larger in the polycentric and river-split cities than in monocentric.

**Q4. Baltimore**
- P4a: At the real-world comparison budget, the optimizer's network beats the stylized real network (Metro SubwayLink + Light RailLink + Red Line) on total rider time saved by 20% or more, because the real lines serve legacy corridors and sparse rowhouse density, and the optimizer concentrates on the downtown-JHH-Homewood-Hopkins corridor.
- P4b: The optimizer builds a north-south spine through downtown, Mount Vernon/Charles Village, and Homewood (the central corridor) plus at least one east-west line serving East Baltimore (JHH) and West Baltimore. It does not reproduce the Owings Mills / Hunt Valley / BWI-heavy suburban tails of the real lines.
- P4c: The Red Line (east-west) accounts for a large share of the real-network benefit gain in the model versus Metro + Light Rail alone (at least 25% increase).
- P4d: Even the optimized network saves a modest absolute amount: transit mode share under 25% of all trips.

**Q5. Triangle**
- P5a: At a budget comparable to the proposed commuter rail (about 37 miles of track, a few stations), the optimizer builds a network that links Raleigh downtown, RTP, and Durham (the densest job pairs), and not primarily Wake Forest or Chapel Hill.
- P5b: The optimizer's network beats the stylized proposed commuter rail corridor on rider time saved by at least 30%, because commuter rail corridors follow the old rail alignment, not demand.
- P5c: Rail does not meaningfully beat driving across the region: transit mode share under 8% of trips even at the high budget, and the share of OD pairs where rail door to door beats driving is under 10%.
- P5d: Rider time saved per capita is at least 5x smaller in the Triangle than in Baltimore at matched budgets (low density, polycentric, long distances).
