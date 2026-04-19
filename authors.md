# Authors and Citation

## Authors

- **[Xingchi Li](https://www.xingchi.li)**. Author, maintainer,
  copyright holder. [](https://orcid.org/0009-0006-2493-0853)

- **[Xianyang Zhang](https://zhangxiany-tamu.github.io)**. Author,
  copyright holder.

## Citation

Source:
[`inst/CITATION`](https://github.com/doccstat/fastcpd/blob/main/inst/CITATION)

Li X, Zhang X (2024). “fastcpd: Fast Change Point Detection in R.”
[doi:10.48550/arXiv.2404.05933](https://doi.org/10.48550/arXiv.2404.05933).

    @Misc{,
      title = {fastcpd: Fast Change Point Detection in R},
      author = {Xingchi Li and Xianyang Zhang},
      year = {2024},
      doi = {10.48550/arXiv.2404.05933},
      publisher = {arXiv},
    }

Zhang X, Dawn T (2023). “Sequential Gradient Descent and Quasi-Newton's
Method for Change-Point Analysis.” In Ruiz, Francisco, Dy, Jennifer, van
de Meent, Jan-Willem (eds.), *Proceedings of The 26th International
Conference on Artificial Intelligence and Statistics*, volume 206 series
Proceedings of Machine Learning Research, 1129–1143.
<https://proceedings.mlr.press/v206/zhang23b.html>.

    @InProceedings{,
      title = {Sequential Gradient Descent and Quasi-Newton's Method for Change-Point Analysis},
      author = {Xianyang Zhang and Trisha Dawn},
      year = {2023},
      booktitle = {Proceedings of The 26th International Conference on Artificial Intelligence and Statistics},
      volume = {206},
      pages = {1129--1143},
      editor = {{Ruiz} and {Francisco} and {Dy} and {Jennifer} and {van de Meent} and {Jan-Willem}},
      series = {Proceedings of Machine Learning Research},
      month = {25--27 Apr},
      publisher = {PMLR},
      pdf = {https://proceedings.mlr.press/v206/zhang23b/zhang23b.pdf},
      url = {https://proceedings.mlr.press/v206/zhang23b.html},
      abstract = {One common approach to detecting change-points is minimizing a cost function over possible numbers and locations of change-points. The framework includes several well-established procedures, such as the penalized likelihood and minimum description length. Such an approach requires finding the cost value repeatedly over different segments of the data set, which can be time-consuming when (i) the data sequence is long and (ii) obtaining the cost value involves solving a non-trivial optimization problem. This paper introduces a new sequential updating method (SE) to find the cost value effectively. The core idea is to update the cost value using the information from previous steps without re-optimizing the objective function. The new method is applied to change-point detection in generalized linear models and penalized regression. Numerical studies show that the new approach can be orders of magnitude faster than the Pruned Exact Linear Time (PELT) method without sacrificing estimation accuracy.},
    }
