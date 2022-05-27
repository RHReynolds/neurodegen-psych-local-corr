
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
[![DOI](https://zenodo.org/badge/375954357.svg)](https://zenodo.org/badge/latestdoi/375954357) <!-- badges: end -->

# Background

This repository contains code to determine local and global genetic correlations between several neurodegenerative and neuropsychiatric disorders with LAVA and LDSC, respectively.

# Code contents

Within this repository you will find:

<table>
<colgroup>
<col width="11%" />
<col width="88%" />
</colgroup>
<thead>
<tr class="header">
<th>Directory</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><a href="docs" class="uri">docs</a></td>
<td>Contains all <code>.Rmd</code>s and their corresponding <code>.html</code>s describing analyses performed for this project. These can be view interactively at: <a href="https://rhreynolds.github.io/neurodegen-psych-local-corr/" class="uri">https://rhreynolds.github.io/neurodegen-psych-local-corr/</a></td>
</tr>
<tr class="even">
<td><a href="logs" class="uri">logs</a></td>
<td>For any scripts that were run outside of an <code>.Rmd</code> (e.g. scripts from the <a href="scripts" class="uri">scripts</a> directory), a log file was recorded and can be accessed here.</td>
</tr>
<tr class="odd">
<td><a href="man" class="uri">man</a></td>
<td>Figures used in <code>.Rmd</code>s.</td>
</tr>
<tr class="even">
<td><a href="R" class="uri">R</a></td>
<td>Various functions called in <a href="docs" class="uri">docs</a> and <a href="scripts" class="uri">scripts</a>.</td>
</tr>
<tr class="odd">
<td><a href="raw_data" class="uri">raw_data</a></td>
<td>External tables used in analyses.</td>
</tr>
<tr class="even">
<td><a href="renv" class="uri">renv</a></td>
<td><code>renv</code>-related scripts-</td>
</tr>
<tr class="odd">
<td><a href="results" class="uri">results</a></td>
<td>Results from all analyses.</td>
</tr>
<tr class="even">
<td><a href="scripts" class="uri">scripts</a></td>
<td>Contains analysis scripts. Each script contains a one-line description and is also referenced in its corresponding <code>.Rmd</code>.</td>
</tr>
</tbody>
</table>

# Reproducibility

This repository uses [`renv`](https://rstudio.github.io/renv/index.html) to create a reproducible environment for this R project.

1.  When you first launches this project, `renv` should automatically bootstrap itself, thereby downloading and installing the appropriate version of `renv` into the project library.
2.  After this has completed, you can use `renv::restore()` to restore the project library locally on your machine.

For more information on collaborating with `renv`, please refer to this [link](https://rstudio.github.io/renv/articles/collaborating.html).

# License

The code in this repository is released under an MIT license. This repository is distributed in the hope that it will be useful to the wider community, but without any warranty of any kind. Please see the [LICENSE](LICENSE) file for more details.

# Citation

If you use any of the code or data from this repository, please cite our [paper](#TODO) and, if applicable, any software dependencies (e.g. [LAVA](https://github.com/josefin-werme/LAVA)).
