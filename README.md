# Auto-Sleep-Screen

Code and data accompanying the manuscript:

**An automated genetic screen identifies modulators of stress-induced sleep in *Caenorhabditis elegans***
Zihao (John) Li *et al.*, Fang-Yen Laboratory, Department of Biomedical Engineering, The Ohio State University.

The manuscript will be deposited on bioRxiv.
Repository: <https://github.com/Zihao-celegans/Auto-Sleep-Screen>

---

## Repository layout

```
.
├── figure1/         Figure 1 plotting script + data
├── figure2/         Figure 2 plotting script + data
├── figure3/         Figure 3 plotting script + data
├── helpers/         Shared MATLAB helper functions used by the figure scripts
├── SKAT/            R/SKAT pipeline used for the association
│                    analysis (inputs, scripts, and result tables)
├── LICENSE          Apache License 2.0
└── README.md
```

Each `figureN/` folder is self-contained: open MATLAB, `cd` into the folder,
and run `Main_figureN.m`. The script will generate the figures shown in the paper.

## Software requirements

**MATLAB** (tested on R2023b or later)
- Statistics and Machine Learning Toolbox (`ranksum`, `vartestn`, `ttest2`)
- Bioinformatics Toolbox (`mafdr`, used for FDR control)

**R** (for the SKAT pipeline only; tested on R 4.3+)
- `SKAT`
- See [SKAT/README.md](SKAT/README.md) for the full list and usage.

## Reproducing the figures

1. Clone the repository.
2. Open MATLAB and set the working directory to one of the `figureN/` folders.
3. Run `Main_figureN.m`. Figures display on screen and are saved as SVG
   files into the same folder. SVG/PDF/PNG outputs are git-ignored.

## License

Apache License 2.0 — see [LICENSE](LICENSE).

## Citation

If you use this code or data, please cite the manuscript (citation will be
added once the bioRxiv preprint is live).

## Contact

Questions or issues: open a GitHub issue, or contact the Fang-Yen Lab at
The Ohio State University.