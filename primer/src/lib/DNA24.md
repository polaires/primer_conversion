# DNA24 Thermodynamic Parameters

This library supports the DNA24 thermodynamic parameters from the Greenleaf Lab (2024), providing significantly improved accuracy for DNA secondary structure and mismatch predictions.

## Source

**Paper**: "High-throughput DNA melt measurements enable improved models of DNA folding thermodynamics"
**Authors**: Ke, W. et al.
**Published**: Nature Communications, 2025
**Data source**: [GreenleafLab/nnn_paper](https://github.com/GreenleafLab/nnn_paper)

## Accuracy Comparison

### Our Implementation vs Original Paper

| Metric | Our Implementation | Original Paper | Status |
|--------|-------------------|----------------|--------|
| **Correlation (R)** | 0.68 | 0.68 | **Matches!** |
| **Bias** | -0.08 kcal/mol | — | Very small |

Our implementation achieves the same correlation (R=0.68) as the original paper for simple hairpin structures.

### DNA24 vs SantaLucia 1998

Validated against 27,732 hairpin sequences from array melt experiments:

| Metric | DNA24 | SantaLucia 1998 | Improvement |
|--------|-------|-----------------|-------------|
| RMSE | 1.68 kcal/mol | 2.58 kcal/mol | **35% lower** |
| MAE | 1.28 kcal/mol | 2.14 kcal/mol | **40% lower** |
| Correlation (R) | 0.52 (all) / 0.68 (simple) | 0.48 | **9-42% higher** |

### Validation Details

**All hairpins** (27,732 sequences):
- R = 0.52
- RMSE = 1.68 kcal/mol
- MAE = 1.28 kcal/mol

**Simple hairpins** (7,897 sequences, no bulges/internal loops):
- **R = 0.68** — matches original paper
- Bias = -0.08 kcal/mol (very small systematic offset)

## Usage

### Default (DNA24 enabled)

DNA24 parameters are used by default. No configuration needed:

```javascript
import { dg, fold, tm } from './index.js';

// All calculations automatically use DNA24 parameters
const deltaG = dg('GCATGCTTTTGCATGC', 37);
const structure = fold('GCATGCTTTTGCATGC', 37);
const meltingTemp = tm('ATCGATCGATCG');
```

### Switching Parameter Sets

Use `setThermodynamicParams()` to switch between parameter sets:

```javascript
import { setThermodynamicParams, dg } from './index.js';

// Use DNA24 (default, recommended)
setThermodynamicParams(true);
console.log(dg('GCATGCTTTTGCATGC', 37)); // DNA24 result

// Use legacy SantaLucia 1998 parameters
setThermodynamicParams(false);
console.log(dg('GCATGCTTTTGCATGC', 37)); // SantaLucia result
```

### Individual Module Control

For fine-grained control, configure each module separately:

```javascript
import { setParameterSet } from './tm.js';
import { setFoldParameterSet } from './fold.js';

// Configure Tm calculations
setParameterSet(true);  // Use DNA24 for Tm

// Configure folding/dG calculations
setFoldParameterSet(true);  // Use DNA24 for folding
```

## DNA24 Parameter Coverage

The DNA24 parameter set includes:

| Parameter Type | Count | Description |
|---------------|-------|-------------|
| Nearest-neighbor (NN) | 42 | Watson-Crick pairs + G-T wobble |
| Internal mismatches | 576 | Context-dependent (6-mer codes) |
| Terminal mismatches | 96 | End-of-helix mismatches |
| Hairpin mismatches | 96 | Hairpin closing mismatches |
| Tetraloops | 1,062 | Stable hairpin loop sequences |
| Triloops | 256 | 3-nucleotide loops |
| Hairpin loops | 30 | Loop length penalties |
| Bulge loops | 30 | Bulge length penalties |
| Internal loops | 30 | Internal loop penalties |
| Dangling ends | 32 | Single-strand overhangs |

## Integration with Primer Design

The following functions automatically use DNA24 parameters:

### primers.js
- `primers()` - PCR primer design with dG scoring
- `create()` - Create primer candidates
- `score()` - Score primer quality

### mutagenesis.js
- `designSubstitutionPrimers()` - Site-directed mutagenesis
- `designInsertionPrimers()` - Insertion mutagenesis
- `designDeletionPrimers()` - Deletion mutagenesis
- `analyzePrimerPair()` - Primer pair analysis

### visualization.js
- Secondary structure visualization uses DNA24 for dG calculations

## Running Validation

To run the validation against experimental data:

```bash
node src/lib/validate-dna24.js
```

This compares DNA24 vs SantaLucia 1998 predictions against measured dG values from the Greenleaf Lab array melt dataset.

## Files

- `dna24.js` - DNA24 parameters (converted from NUPACK JSON format)
- `dna24.json` - Original NUPACK parameter file
- `dna.js` - SantaLucia 1998 parameters (legacy)
- `tm.js` - Tm calculations with parameter switching
- `fold.js` - Secondary structure folding with parameter switching
- `validate-dna24.js` - Validation script
- `arr_v1_1M_n=27732.csv` - Experimental validation data (27,732 sequences)

## Citation

If you use the DNA24 parameters, please cite:

```
Ke, W. et al. High-throughput DNA melt measurements enable improved models
of DNA folding thermodynamics. Nature Communications (2025).
https://doi.org/10.1038/s41467-025-60455-4
```

## References

- [Nature Communications Paper](https://www.nature.com/articles/s41467-025-60455-4)
- [GreenleafLab GitHub](https://github.com/GreenleafLab/nnn_paper)
- [NUPACK Documentation](https://docs.nupack.org/model/)
