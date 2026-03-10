# mip-poolfreq
![License](https://img.shields.io/github/license/GE2POP/mip-poolfreq)
![CI](https://github.com/GE2POP/mip-poolfreq/actions/workflows/ci.yml/badge.svg)
[![Container](https://img.shields.io/badge/container-ghcr.io%2Fge2pop%2Fmip--poolfreq-blue)](https://github.com/GE2POP/mip-poolfreq/pkgs/container/mip-poolfreq)

Tools to estimate allele frequencies and genotype frequencies from pooled sequencing data generated with Molecular Inversion Probes (MIPs).

This repository contains two complementary command-line tools:

- ***Afreq*** – extract reference allele frequencies and read depths from VCF files

- ***Gfreq*** – estimate genotype frequencies in pooled samples using allele frequencies and sequencing depth information

The typical analysis workflow is illustrated below.

<pre>
VCF with AD field
      │
      ▼
<em><strong>Afreq estimate</strong></em>
      │
      ├── ref allele frequencies
      └── total depths
              │
              ▼
        <em><strong>Gfreq estimate</strong></em>
              │
              ▼
   genotype frequency estimates
</pre>

Both tools are distributed together in a container to simplify reproducible analyses.


## Container image
A pre-built OCI image is available on GitHub Container Registry:
```bash
ghcr.io/ge2pop/mip-poolfreq
```

It can be pulled with Apptainer to obtain a local .sif image:
```bash
apptainer pull mip-poolfreq.sif docker://ghcr.io/ge2pop/mip-poolfreq:latest
```

Show available commands:
```bash
apptainer exec mip-poolfreq.sif --help
```


## Tools
### Afreq

Compute reference allele frequencies and total read depths from the AD field of a VCF file.

Typical use case: transform raw variant calls into allele frequency matrices that can be used as input for downstream analyses.

Documentation: [Afreq/README.md](Afreq/README.md)

### Gfreq

Estimate genotype frequencies in pooled samples from allele frequencies and sequencing depth data produced by **Afreq**.

Includes several commands:

- estimate – estimate genotype frequencies

- eval – evaluate estimation accuracy in controlled mixtures with known composition

- maf – compute minor allele frequencies

- depth – summarize sequencing depth per marker

Documentation: [Gfreq/README.md](Gfreq/README.md)


## Citation

If you use **mip-poolfreq**, please cite [...]


## License

MIT License
