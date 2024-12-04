This work is based on the referenced paper, where I attempt to reproduce most of the results. Unlike the paper, I use only 100 channel realizations. For data detections for SER analysis, I use only one channel realization. I have included a PDF of my walkthrough of the derivations in the paper.

The results are compiled in result.pdf.

This project is for educational purposes, aiming to enhance my ability to understand academic papers and improve my simulation skills.

The modules included are as follows:

maincode.m: Contains the main script that calls other functions to produce the data presented in Figures 2, 4, and 5 of the referenced paper.
channel.m: A function that returns a discrete physical channel based on UPAs at the transmitter (Tx) and receiver (Rx).
combiner.m: A function that returns the combiner V to be used at the receiver for estimation.




I. Atzeni, A. Tölli, D. H. N. Nguyen and A. L. Swindlehurst, "Doubly 1-Bit Quantized Massive MIMO," 2023 57th Asilomar Conference on Signals, Systems, and Computers, Pacific Grove, CA, USA, 2023, pp. 465-469, doi: 10.1109/IEEECONF59524.2023.10476782. keywords: {Radio frequency;Power demand;Costs;Error analysis;Symbols;Massive MIMO;Frequency conversion;1-bit ADCs;1-bit DACs;doubly massive MIMO;(sub-)THz communications},

