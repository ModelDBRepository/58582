# Model of a Thalamic Network, According to the Paper

Golomb D, Ahissar E, Kleinfeld D.
Coding of stimulus frequency by latency in thalamic networks through
the interplay of GABAB-mediated feedback and stimulus shape.
*J Neurophysiol.* 2005 Nov 2; [Epub ahead of print]

## Program Composition

The program is composed of several C programs and files:

- **tlp.c**: main simulation program.
- **tplp.c**: part of the tlp.c program. It separates the columns of tlp.col.?? and calculates properties of the "PSTHs".
- **plp.c**: independent program. It separates the columns of tlp.col.??.
- **stm.c**: used for plotting the stimuli.
- **intnor.c**: Integrate the activity trace - modeling results (to compute the normalized firing rates).
- **nr.c**: routines from the library "Numerical Recipes" (only a few of them are needed for this program).

The file **tlp.n.d1** is the parameter file for tlp.c.
The file **tlp.nps** is the parameter file for tplp.c.
The file **plp.n** is the parameter file for plp.c.

## Experimental Results (from Sosnik et al., *J. Neurophysiol.*, 86:339-353, 2001)

- **intg.res.a3** - integrated PSTH for VPm.
- **intg.res.a4** - integrated PSTH for POm.
- **delay_exp** - delays.

The file **ptrint.com** is a Unix/Linux command procedure.

## Commands (for a Unix/Linux Operating System)

To compile the programs:

```
make tlp
make plp
make stm
make intnor
```

To run the programs and generate figure 7A-H in the paper:

```
ptrint.com
```

You should have the graphic program **xmgrace**.

If you have any questions or problems, E-mail me at golomb@bgu.ac.il.

David Golomb

---

2025-05-27 – Standardized to Markdown.