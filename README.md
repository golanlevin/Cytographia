# Cytographia
*An algorithmic neoincunabulum of xenocytology*

---
### About *Cytographia*

*Cytographia* is an elegy for species we will never know, or will never know again, expressed through generative illustrations from an imaginary book about imaginary microorganisms. The artwork depicts speculative cell structures whose appearances and movements arise emergently and in response to real-time user interactions. *Cytographia* was created using [p5.js](https://p5js.org/) (an open-source programming toolkit for the arts) and may run best in the Chrome browser on OSX. 

An algorithmic "neoincunabulum of xenocytology", *Cytographia* presents an interactive diagram of a one-celled organism, styled to evoke a hand-drawn engraving. Every aspect of this illustration is generated through custom code, including the simulated behavior and anatomy of the depicted creature, the calligraphic ductus of its lines, the asemic letterforms of its labels, and the virtual "paper" on which it is rendered. Each token in the *Cytographia* collection comprises a signature of sixteen such illustration plates; in addition to the primary "cover" page, the other fifteen pages of a token's folios can be browsed by clicking in the corners of the display. *Cytographia*'s drawings may be exported as high-resolution PNG images, or in an SVG vector format suitable for pen-plotting on A4 paper. 

The *Cytographia* project draws inspiration from several historically significant books that visually and methodically documented encounters with the unknown. These include Robert Hooke's *[Micrographia](https://en.wikipedia.org/wiki/Micrographia)* (1665), a landmark of scientific observation in which living cells were described for the first time; Edmund Fry's *[Pantographia](https://en.wikipedia.org/wiki/Pantographia)* (1799), an attempt to compile exemplars of all the world's writing systems; Ernst Haeckel's [*Kunstformen der Natur*](https://en.wikipedia.org/wiki/Kunstformen_der_Natur) (1899), a rich exploration of symmetry and structural hierarchy in natural forms; and Luigi Serafini's hallucinatory *[Codex Seraphinianus](https://en.wikipedia.org/wiki/Codex_Seraphinianus)* (1981), a visual encyclopedia of the artist's imagined world. 

*Cytographia* also brings together several long-time threads in Golan Levin's thirty-year *oeuvre* of software art, including research into interactive blobs (e.g. [*Polygona Nervosa*](https://objkt.com/asset/hicetnunc/56312), 1997); artificial life (e.g. [*Obzok*](https://www.youtube.com/watch?v=oVOCKzE2fZk), 2001), and the computational generation of asemic writing systems (e.g. [*Alphabet Synthesis Machine*](https://en.wikipedia.org/wiki/Alphabet_Synthesis_Machine), 2002). The generative letterforms in *Cytographia* are loosely based on 16th-century italic typefaces by [Ludovico degli Arrighi](https://en.wikipedia.org/wiki/Ludovico_Vicentino_degli_Arrighi). 

---
### Display Notes

*Cytographia* requires a modern browser with WebGL and hardware acceleration enabled. The recommended configuration for *Cytographia* is the Chrome browser on MacOS, in which the artwork is known to produce consistent and replicable results.

An index of *Cytographia*'s key commands can be displayed by pressing **h** (for "help"). These key commands provide access to functionality including file export; page navigation and zooming; an "auto-play" mode that periodically cycles through the token's pages; toggles that enable or disable various graphical options, potentially improving performance on some systems; and access to a playful "sandbox" mode, in which a collector or visitor can assemble the organelles of their own imaginary lifeform.

*Cytographia* was developed using JavaScript version Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) with AppleWebKit/537.36 (KHTML, like Gecko), and was tested in Google Chrome version 118.0.5993.70 (arm64) on a MacBook Pro (14-inch, 2021) with an Apple M1 Pro CPU, running macOS 12.6 (21G115). Please note that [the precision of the JavaScript Math library is implementation-dependent](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math). This means that browsers and/or operating systems other than the recommended configuration may produce different results for the same series of mathematical operations. The Cytographia project is particularly susceptible to this, since it constantly simulates millions of interactions between thousands of particles. Over time, these environment-based differences may accumulate to produce visibly different results.

---
### Acknowledgements

*Cytographia* incorporates or adapts the following code under the specified licenses: [*p5.js v.1.0.0*](https://p5js.org/) by The Processing Foundation (GPL); [*Fortune's Voronoi*](https://github.com/d3/d3-delaunay) by Mike Bostock (Observable/Mapbox) (ISC); [*GLSL Blend Modes*](https://github.com/jamieowen/glsl-blend) by Jamie Owen (MIT); [*GLSL Value Noise*](https://www.shadertoy.com/view/lsf3WH) and [*GLSL XOR Noise*](https://www.shadertoy.com/view/XtXXD8) by Inigo Quilez (MIT); [*2D Perlin Noise*](https://github.com/stegu/webgl-noise/blob/master/src/classicnoise2D.glsl) by Stefan Gustavson (MIT); [*ImageJ (Distance Transform)*](https://github.com/imagej/ImageJ/blob/master/ij/process/BinaryInterpolator.java) by NIH (Public Domain); [*ofPolyline*](https://github.com/openframeworks/openFrameworks/tree/master/libs/openFrameworks/graphics) in openFrameworks (MIT); and [*webgl-lines*](https://mattdesl.github.io/webgl-lines/expanded/gl-line-2d.js) by Matt DesLauriers (MIT). *Cytographia* additionally adapts the following code and/or media, whose authors appear not to have specified a license: [*GLSL p5jsShaderExamples*](https://github.com/aferriss/p5jsShaderExamples/) by Adam Ferriss; [*GLSL LA Font*](https://github.com/hi-ogawa/python-shader-app/tree/master/misc/la_font) by Hi Ogawa; [*Delaunay Triangulation*](https://editor.p5js.org/allison.parrish/sketches/BkhEmKKjW) by Allison Parrish; [*Closest Point of a Polygon*](https://codesandbox.io/s/elated-liskov-3v65c?file=/src/getClosestPointInsidePolygon.ts) by Daniel Neveux; [*Point, Line, Plane*](http://paulbourke.net/geometry/pointlineplane/) by Paul Bourke; [*Flocking/Boids*](https://p5js.org/examples/simulate-flocking.html) by Craig Reynolds, adapted by Daniel Shiffman; and [*Arrighi*](https://en.wikipedia.org/wiki/Ludovico_Vicentino_degli_Arrighi) (1523) by Ludovico Vicentino degli Arrighi. The following writings were also essential references for this work: [*Introduction to p5.js shaders*](https://itp-xstory.github.io/p5js-shaders/#/) (2019) by Casey Conchinha and Louise Lessél; [*Drawing Lines is Hard*](https://mattdesl.svbtle.com/drawing-lines-is-hard) (2015) by Matt DesLauriers; and [*Entering the Blobosphere: A Musing on Blobs*](http://www.lauraonsale.com/blob.html) (2019) by Laura Hyunjhee Kim. Special thanks to Zach Lieberman, Casey Reas, James Paterson, and Lingdong Huang.