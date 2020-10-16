[node distance=.5cm, start chain=going below,]

(init\_global) initial guess; (sampling) sampling; (pot\_update)
calculate potential update; = [diamond, draw = black, very thick, text
width=8em, on chain, text badly centered, node distance=4cm, inner
sep=0pt]

(continue) converged?; (end) done;

(positive) ;

(continue.west) \|- (positive.east); (positive.north) \|-
(sampling.west);

(a) at ([xshift=12pt, yshift=-4pt] positive.east) yes; (b) at
([xshift=-4pt] continue.south) no;
