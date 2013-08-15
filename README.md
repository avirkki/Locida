<img align="right" src="Documentation/images/vttplain.png" />

Locida
======

An **R package** for NMR Local Intensity Difference Analysis

What is Locida Good For
--------------------

[Nuclear magnetic resonance (NMR) spectroscopy](http://en.wikipedia.org/wiki/Nuclear_magnetic_resonance_spectroscopy) is a measurement technique that utilises the magnetic properties of the atomic nuclei. NMR technology has became famous for its robustness and repeatability for assessing the amount of chemical compounds in various types of biomedical samples.

**Locida** is an *R software package* and a *Web tool* for quick and illustrative *group-wise comparisons of 1-D NMR spectral measurements*. The mathematical method, *local intensity difference analysis*, is also called Locida, since the differences in signals are found by inspecting the original (binned) spectral measurements. *The method seeks diffences between groups of samples*, and illustrates them on a graphical display showing also the original data for each individual sample. 

A running prototype of Locida can be tested at [bicbox.vtt.fi](http://bicbox.vtt.fi:8080/Locida).

![Locida_overview.png](Documentation/images/Locida_overview.png?raw=true)
**Figure.** An example Locida output comparing 1-D NMR spectrum measurements from treated Nicotiana Tabacum plant samples agains two different (but indistinguishable) control groups. (Click the image to enlarge it.)

*The method* was originally developed by [A Virkki](http://fi.linkedin.com/in/arhovirkki), [H Maaheimo](http://www.researchgate.net/profile/Hannu_Maaheimo/) and colleagues at [VTT Technical Research Centre of Finland](http://www.vtt.fi/?lang=en) to study the differences between different plant cultures in EU 7th Framework Research Programme [SmartCell](http://www.smart-cell.org/). The description of the method will be published in due course of that project, whereas interested users can already download the [Locida R Package](Rpkg/Locida_1.0.tar.gz?raw=true) for R command-line analyses, or set up their own Web-based installation by downloading the Locida UI Java Package by cloning the Locida project, and setting up the needed [RVaadin](https://github.com/avirkki/RVaadin) run-time environment to launch it. Locida Web interface is implemented in Java by using the [Vaadin Web framework](http://vaadin.com), and the inteface communicates with the [R language](http://www.r-project.org) through the [RVaadin](https://github.com/avirkki/RVaadin) library. Both the R package and the Web interface are free software and published under GPL v2.



