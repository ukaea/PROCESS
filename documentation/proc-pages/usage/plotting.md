# Utilities

<h2>Plotting an MFILE</h2>

<h3>Summary document</h3>
`plot_proc` is used for plotting an overview of the results from an MFILE. It can be run using its own CLI:

```bash
python process/io/plot_proc.py -f path/to/MFILE.DAT
```

or through Process's main CLI (working, but still in development):

```bash
process -i path/to/IN.DAT --plot --mfile path/to/MFILE.DAT
``` 
<h3>Radial build</h3>

`plot_radial_build` is to plot the radial build of the machine in terms of bar segments. It can be run as follows:

```bash
python process/io/plot_radial_build.py -f path/to/MFILE.DAT
```
<figure markdown>
![radial_build_plot](../../images/radial_build_plot.png){ width="100%"}
<figcaption>Figure 1: Simple radial build plot </figcaption>
</figure>

