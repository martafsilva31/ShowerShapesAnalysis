I# Cell Collision Study — Layer 2 (7×11 Grid)

Answers to supervisor's questions about cell collisions in the
CaloFillRectangularCluster 7×11 grid output.

## Context

When CaloFillRectangularCluster extracts cells in a 7η × 11φ window
(0.175 × 0.275), the result is reordered onto a regular grid with cell
size 0.025 × 0.025. At the barrel–endcap transition (~1.3 < |η| < 1.55),
the calorimeter geometry changes and more than 77 raw cells can fall
inside the window. `reorderToGrid()` handles this by **accumulating**
energies into the grid bin: `E_out[idx] += E[i]`.

Mohamed found that **topocluster cell links** can contain the same
physical cell shared between multiple clusters with fractional weights
(`CaloClusterCellLink::iterator::weight()`). The question is whether
the same mechanism applies here.

## Key Finding: It Does Not

| Topocluster cell links | CaloFillRectangularCluster |
|------------------------|----------------------------|
| Same physical cell shared between clusters | Fresh cell extraction from CaloCellContainer |
| Fractional weights (< 1.0) | All weights = 1.0 |
| Duplicates = same cell with split energy | "Collisions" = different cells at different positions mapping to same grid bin |

## Questions & Answers

### Q1: How often do collisions happen?

Run `analyse_collisions.py` for data-driven numbers. Results from
v9 ntuple (4145 photons, 30k EGAM3 events):

- **89/4145 events (2.1%)** have Lr2Size > 77 (more raw cells than grid slots)
- All have exactly **Lr2Size = 88** (11 extra cells = one full phi row)
- **98.9%** occur in the barrel–endcap crack (1.3 ≤ |η| ≤ 1.55)
- |η| distribution: 64% in [1.30, 1.37), 35% in [1.52, 1.60)

### Q2: Are cells weighted? Should we use weights?

**No.** `CaloFillRectangularCluster` builds a new cluster from the full
`CaloCellContainer` — every cell has weight 1.0. The `weight()` method
on `CaloClusterCellLink::iterator` exists but is never called in
NTupleMaker, and would return 1.0 anyway for these cells.

Mohamed's weighted-duplicate scenario is specific to **topocluster cell
links**, where the same physical cell can be shared between overlapping
topoclusters with fractional energy sharing. That does not happen with
CaloFillRectangularCluster, which extracts cells independently.

### Q3: CaloFillRectangularCluster or cluster cell links?

**Both, in two stages:**

1. **Stage 1 — Cell links**: Iterate the original cluster's
   `CaloClusterCellLink` to build a flat `vector<const CaloCell*>`.
   Used only to find the seed cell position.

2. **Stage 2 — CaloFillRectangularCluster**: Build a fresh rectangular
   cluster centred on the seed. Iterate its cell links to extract
   (E, η, φ) for each cell. This is the source of the 7×11 grid arrays.

### Q4: Are duplicates identical (E, η, φ) or just same grid bin?

They are **different physical cells** at different (η, φ) positions that
happen to fall into the same 0.025 × 0.025 grid bin. This occurs at the
barrel–endcap crack where the calorimeter cell geometry is irregular.

**0 exact duplicates** found (after excluding zero-padded empty grid slots).
All 89 collision events have distinct physical cells — confirmed by
`findDuplicatedCells()` in NTupleMaker and by this study.

Energy accumulation (`+=`) is the correct treatment: these are real cells
contributing real energy to the same region of the detector.

## How to Run

```bash
cd /project/atlas/users/mfernand/QT/ShowerShapes/ShowerShapesAnalysis
lsetup "root 6.32.08-x86_64-el9-gcc13-opt"
python scripts/cell_collision_study/analyse_collisions.py
```

## Code References

- `reorderToGrid()`: [NTupleMaker.cxx L1548–1596](../../NTupleMaker_workspace/source/NTupleMaker/src/NTupleMaker.cxx)
- `findDuplicatedCells()`: [NTupleMaker.cxx L1495–1511](../../NTupleMaker_workspace/source/NTupleMaker/src/NTupleMaker.cxx)
- `egammaCellDecorator.cxx` Stage 2: [L183–212](../../NTupleMaker_workspace/source/NTupleMaker/src/egammaCellDecorator.cxx)
