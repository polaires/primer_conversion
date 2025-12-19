# ggtools

## Descrption

This is a PERL script to idenity an optimal set of overhangs of the required size given experimental ligation assembly data.

## Usage examples

* Create a 20-overhang set using BsaI ligation assembly data

```
ggtools_optimize.pl --size 20 --overhang_size 4 BsaI.csv
```

* Create a 20-overhang set using BsaI ligation assembly data AND require that one of the overhangs in the resulting set is TGGC

```
ggtools_optimize.pl --size 20 --overhang_size 4 --olist TGGC BsaI.csv
```

* Create a 20-overhang set using BsaI ligation assembly data AND require that one of the overhangs in the resulting set contains TGGC or GGAT or TTGG

```
ggtools_optimize.pl --size 20 --overhang_size 4 --olist TGGC,GGAT,TTGG BsaI.csv
```

* Create a 20-overhang set using BsaI ligation assembly data AND require that one of the overhangs in the resulting set is either TGGC or GGAT or TTGG AND another overhang in the resulting set is AGTC or CCGC

```
ggtools_optimize.pl --size 20 --overhang_size 4 --olist TGGC,GGAT,TTGG --olist AGTC,CCGC BsaI.csv
```

## Citations

* Pryor J.M., Potapov V., Kucera R.B., Bilotti K., Cantor E.J. and Lohman G.J.S. (2020) Enabling one-pot Golden Gate assemblies of unprecedented complexity using data-optimized assembly design. PLOS ONE, 15(9):e0238592. doi: 10.1371/journal.pone.0238592.
