---
name: Bug report
about: Submit a bug report to help us improve fastcpd
---

### Tips for a helpful bug report:

* Please include a **minimal reproducible example** to demonstrate the bug.
  A reprex can be helpful. See
  [make a reprex](https://www.tidyverse.org/help/#reprex).
  Do not include session info unless it is explicitly asked for.

* If you can, use a small toy dataset that exposes the bug. If for some reason,
  the bug only occurs on your original data, try to limit the number of rows that
  are necessary to expose the bug. Share such data by copying `dput()` output
  rather than an external file.

* Unless the bug is about the theme, labels, scales or other plot decoration:
  please omit these from the code.

* Please check whether somebody has reported the same problem in the
  [issues](https://github.com/doccstat/fastcpd/issues).

Delete these instructions once you have read them.

---

I found a problem with ...

I expected ...

Here is the code to reproduce the bug:

```r
# copy your code to the clipboard and run:
reprex::reprex()
```

<!-- Copy and modified based on https://github.com/tidyverse/ggplot2/tree/main/.github/ISSUE_TEMPLATE -->
