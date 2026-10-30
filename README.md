# cpDraw

**cpDraw v1.0.0** is a stand-alone Python command-line program for generating editable circular plastid genome maps from annotated GenBank files. It runs locally and produces SVG output without requiring a web server.

cpDraw was developed from the former `plastome_draw` codebase. SVG metadata and optional diagnostic CSV files preserve the lineage as `plastome_draw v1.1.8`.

## Main features

- Circular plastid genome maps from annotated GenBank files
- Gene blocks displayed according to transcriptional orientation
- LSC, SSC, IRA, and IRB labels when recognizable annotations are present
- Multipart and trans-spliced feature handling
- Suppression of redundant labels from overlapping GenBank annotations
- Asterisks for multipart or spliced features
- Optional pseudogene prefix (`psi`)
- Automatic gene-label placement and leader-line routing
- GC-content display on a 0–100% radial scale
- Solid 50% GC reference circle and dashed genome-wide mean-GC circle
- Optional genomic ruler with 1-kb ticks, 5-kb medium ticks, and 10-kb labels
- Editable SVG output with embedded provenance metadata
- Optional label-placement diagnostic CSV

## Requirements

- Python 3
- Biopython

Install Biopython with:

```bash
python -m pip install biopython
```

On Windows, the following command may also be used:

```bash
py -m pip install biopython
```

cpDraw creates SVG output only. CairoSVG, Cairo, Poppler, and `pdftops` are not required.

## Installation

Download `cpdraw.py` and place it in a working directory. Confirm that the program starts correctly:

```bash
python cpdraw.py --help
```

## Quick start

### Use the default input filename

If the current directory contains `sequence.gb`:

```bash
python cpdraw.py
```

Output:

```text
sequence.svg
```

### Specify a GenBank file

```bash
python cpdraw.py plastome.gb
```

Output:

```text
plastome.svg
```

### Specify an output filename

```bash
python cpdraw.py plastome.gb -o figure.svg
```

### Add the genomic coordinate ruler

```bash
python cpdraw.py plastome.gb --show-scale
```

The ruler contains 1-kb minor ticks, 5-kb medium ticks, and labeled 10-kb major ticks.

### Show pseudogene prefixes

```bash
python cpdraw.py plastome.gb --show-pseudogene
```

### Remove gene-block outlines

```bash
python cpdraw.py plastome.gb --block-outline-width 0
```

### Create a diagnostic CSV with the automatic filename

```bash
python cpdraw.py plastome.gb --label-report
```

Outputs:

```text
plastome.svg
plastome_label_report.csv
```

### Create a diagnostic CSV with a specified filename

```bash
python cpdraw.py plastome.gb --label-report diagnostics.csv
```

## Sample files

The public repository includes a sample GenBank file and its cpDraw output:

```text
examples/
├── sequence.gb
├── sequence.svg
└── sequence_label_report.csv.
```

Generate the sample map with:

```bash
python cpdraw.py examples/sequence.gb
```

Generate the sample map with the genomic coordinate ruler:

```bash
python cpdraw.py examples/sequence.gb --show-scale
```

Generate the optional diagnostic CSV using the automatic report name:

```bash
python cpdraw.py examples/sequence.gb --label-report
```

This creates:

```text
examples/sequence.svg
examples/sequence_label_report.csv
```

To avoid overwriting the distributed sample SVG while testing other settings, specify another output name:

```bash
python cpdraw.py examples/sequence.gb -o sequence_test.svg --show-scale
```

## Input

cpDraw accepts one annotated sequence in GenBank flat-file format. If no input filename is supplied, `sequence.gb` is used.

The program reads gene-related annotations from `gene`, `CDS`, `tRNA`, and `rRNA` features. Region annotations such as LSC, SSC, IRA, and IRB are displayed when recognizable annotations are available in the GenBank record.

## Output

### Standard output

- SVG genome map (`*.svg`)

### Optional output

- Label-placement diagnostic CSV, enabled with `--label-report`

SVG metadata records:

- program name and version
- legacy program name and version
- GenBank record ID
- organism
- genome length
- input filename

The diagnostic CSV records `cpDraw 1.0.0` and the `plastome_draw 1.1.8` lineage in every row.

## Public command-line options

| Option | Default | Description |
|---|---:|---|
| `genbank` | `sequence.gb` | Input GenBank file |
| `-o`, `--output` | input basename | Output SVG path or basename |
| `--size` | `1200` | SVG canvas size |
| `--font-size` | `10.0` | Gene-label font size |
| `--block-outline-width` | `0.25` | Gene-block outline width; `0` removes outlines |
| `--ir-region-width` | `3.0` | IRA/IRB width on the inner region ring; `0` disables the overlay |
| `--ir-backbone-width` | `6.0` | IRA/IRB width on the gene-block backbone; `0` disables thickening |
| `--show-pseudogene` | off | Prefix pseudogene labels with `psi` |
| `--show-scale` | off | Add the genomic coordinate ruler |
| `--leader-width` | `0.25` | Leader-line stroke width |
| `--leader-label-clearance` | `4.2` | Minimum leader-to-label clearance |
| `--leader-line-clearance` | `2.0` | Minimum clearance between leader lines |
| `--leader-promotion-clearance` | `6.0` | Promote a natural label when another leader enters this corridor |
| `--leader-elbow-clearance` | `5.0` | Minimum separation between leader elbows |
| `--label-report [FILE]` | off | Write diagnostic CSV; an omitted `FILE` uses `<input>_label_report.csv` |
| `--collision-angle-step` | `0.18` | Tangential search step in degrees |
| `--max-label-shift` | `7.0` | Maximum tangential shift before radial fallback |
| `--collision-padding` | `1.2` | Padding for rotated-label collision detection |
| `--gc-bars` | `3600` | Number of radial GC bars |
| `--gc-bar-width` | `1.0` | GC-bar stroke width at 1200 px |
| `--trans-gap` | `10000` | Gap threshold in bp for trans-spliced feature handling |
| `--center-text` | record-derived | Override center text |
| `-h`, `--help` | n/a | Show help and usage examples |

Run the following command for the authoritative option list in the installed version:

```bash
python cpdraw.py --help
```

## Multipart features and asterisks

Individual components of multipart features are represented in the map, while one label is normally assigned per biological feature and strand. An asterisk indicates a multipart or spliced feature. Independent same-name copies at different genomic positions are retained separately. Consequently, an intron-containing copy and an intronless copy of the same gene name may be displayed with and without an asterisk, respectively.

## Diagnostic CSV

The diagnostic CSV is intended for checking label placement and leader-line decisions. It is not generated during standard execution.

Use the automatic name:

```bash
python cpdraw.py plastome.gb --label-report
```

Or specify the report name:

```bash
python cpdraw.py plastome.gb --label-report diagnostics.csv
```

## Known limitations

- In locally dense annotation regions, automatic leader-line routing may occasionally produce visually suboptimal bends.
- GenBank annotation practices differ among records. Outputs should be checked against the source annotation, especially for partial, linear, legacy, or unusually structured records.
- SVG output can be refined in vector-graphics software such as Inkscape, Adobe Illustrator, or Affinity Designer when final manual adjustment is needed.
- PDF and EPS conversion are intentionally outside the scope of cpDraw v1.0.0.

## Citation

Citation metadata are provided in `CITATION.cff`.

DOI: https://doi.org/10.5281/zenodo.23771064

If you use cpDraw, please cite the archived software release:

Tsubota H. (2026) cpDraw: a stand-alone Python tool for drawing circular plastid genome maps. Version 1.0.0. Zenodo. https://doi.org/10.5281/zenodo.23771064

## Acknowledgements

The author thanks Dr. Phan Quynh Chi and Dr. Miho Nakahara-Tsubota for testing cpDraw with plastid genome data and for providing helpful comments on its graphical output during development. The author also thanks Dr. Yuya Inoue for providing the *Timmiella anomala* plastid genome data ([LC761301]https://www.ncbi.nlm.nih.gov/nuccore/LC761301.1/) used for testing cpDraw and as the example dataset distributed with the program.

## Supporting research and humanitarian activities

cpDraw is freely available under the MIT License. If you find this software useful, please consider supporting scientific research, natural history collections, or humanitarian activities through donations to universities, museums, botanical institutions, or organizations of the International Red Cross and Red Crescent Movement.

Such donations are entirely voluntary and are not required to use cpDraw.

## License

cpDraw is distributed under the MIT License.

Copyright (c) 2026 Hiromi Tsubota

See `LICENSE` for the complete license text.

## Release files

A release repository should contain at least:

```text
cpDraw/
├── cpdraw.py
├── README.md
├── LICENSE
├── CITATION.cff
├── requirements.txt
└── examples/
    ├── sequence.gb
    ├── sequence.svg
    └── sequence_label_report.csv
```

## Tested environment

cpDraw v1.0.0 was tested with:

- Windows 11
  - Python 3.12.1
  - Biopython 1.83

- Linux Mint 22.3
  - Python 3.9.18
  - Biopython 1.85

## Contact

Hiromi Tsubota  
Hiroshima University

Add a public contact email address before release if desired.
