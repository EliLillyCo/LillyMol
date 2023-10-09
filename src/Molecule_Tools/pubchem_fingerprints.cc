// clang-format off
/*
  Compute Pubchem Fingerprints per

PubChem Substructure Fingerprint
V1.3							   http://pubchem.ncbi.nlm.nih.gov

The PubChem System generates a binary substructure fingerprint for 
chemical structures.  These fingerprints are used by PubChem for 
similarity neighboring and similarity searching.

A substructure is a fragment of a chemical structure.  A 
fingerprint is an ordered list of binary (1/0) bits.  Each bit 
represents a Boolean determination of, or test for, the presence 
of, for example, an element count, a type of ring system, atom 
pairing, atom environment (nearest neighbors), etc., in a chemical 
structure.

The native format of the PubChem Substructure Fingerprint property 
is binary data with a four byte integer prefix, where this integer 
prefix indicates the length of the bit list.  For the ASN.1 and XML 
formatted data, this property is stored in a PC-InfoData container, 
as described by the PCSubstance ASN.1 definition or XML schema:
ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/

PC-InfoData is able to handle various types of data.  Each PC-
InfoData has a PC-Urn object (urn = universal resource name). Each 
property has a unique trio of "label", "name", and "datatype" 
definition (e.g., for PubChem Substructure Fingerprint, this is 
"Fingerprint", "SubStructure Keys", and "fingerprint", 
respectively).  The fingerprint binary data is hex-encoded, when 
provided in the XML or textual ASN.1 formats.

When exporting fingerprint information in the SD file format, the 
SD tag for the PubChem Substructure Fingerprint property is 
"PUBCHEM_CACTVS_SUBGRAPHKEYS".  The PubChem Substructure 
Fingerprint is Base64 encoded to provide a textual representation 
of the binary data.  For a description of the Base64 encoding and 
decoding algorithm specification, go to:
http://www.faqs.org/rfcs/rfc3548.html

Below is the description of each bit represented in the PubChem 
Substructure Fingerprint.  Some fingerprint bit descriptions are 
written in SMILES or SMARTS notation.  For additional information 
on SMARTS and SMILES, please go to:
http://en.wikipedia.org/wiki/Simplified_molecular_input_line_entry_
specification




PubChem Substructure Fingerprint Description


Section 1: Hierarchic Element Counts - These bits test for the 
presence or count of individual chemical atoms 
represented by their atomic symbol.

Bit Position	Bit Substructure
0	>= 4 H
1 	>= 8 H
2	>= 16 H
3	>= 32 H
4	>= 1 Li
5	>= 2 Li
6	>= 1 B
7	>= 2 B
8	>= 4 B
9	>= 2 C
10	>= 4 C
11	>= 8 C
12	>= 16 C
13	>= 32 C
14	>= 1 N
15	>= 2 N
16	>= 4 N
17	>= 8 N
18	>= 1 O
19	>= 2 O
20	>= 4 O
21	>= 8 O
22	>= 16 O
23	>= 1 F
24	>= 2 F
25	>= 4 F
26	>= 1 Na
27	>= 2 Na
28	>= 1 Si
29	>= 2 Si
30	>= 1 P
31	>= 2 P
32	>= 4 P
33	>= 1 S
34	>= 2 S
35	>= 4 S
36	>= 8 S
37	>= 1 Cl
38	>= 2 Cl
39	>= 4 Cl
40	>= 8 Cl
41	>= 1 K
42	>= 2 K
43	>= 1 Br
44	>= 2 Br
45	>= 4 Br
46	>= 1 I
47	>= 2 I
48	>= 4 I
49 	>= 1 Be
50	>= 1 Mg
51	>= 1 Al
52	>= 1 Ca
53	>= 1 Sc
54	>= 1 Ti
55	>= 1 V
56	>= 1 Cr
57	>= 1 Mn
58	>= 1 Fe
59	>= 1 Co
60	>= 1 Ni
61	>= 1 Cu
62	>= 1 Zn
63	>= 1 Ga
64	>= 1 Ge
65	>= 1 As
66	>= 1 Se
67	>= 1 Kr
68	>= 1 Rb
69	>= 1 Sr
70	>= 1 Y
71	>= 1 Zr
72	>= 1 Nb
73	>= 1 Mo
74	>= 1 Ru
75	>= 1 Rh
76	>= 1 Pd
77	>= 1 Ag
78	>= 1 Cd
79	>= 1 In
80	>= 1 Sn
81	>= 1 Sb
82	>= 1 Te
83	>= 1 Xe
84	>= 1 Cs
85	>= 1 Ba
86	>= 1 Lu
87	>= 1 Hf
88	>= 1 Ta
89	>= 1 W
90	>= 1 Re
91	>= 1 Os
92	>= 1 Ir
93	>= 1 Pt
94	>= 1 Au
95	>= 1 Hg
96	>= 1 Tl
97	>= 1 Pb
98	>= 1 Bi
99	>= 1 La
100	>= 1 Ce
101	>= 1 Pr
102	>= 1 Nd
103	>= 1 Pm
104	>= 1 Sm
105	>= 1 Eu
106	>= 1 Gd
107	>= 1 Tb
108	>= 1 Dy
109	>= 1 Ho
110	>= 1 Er
111	>= 1 Tm
112	>= 1 Yb
113	>= 1 Tc
114	>= 1 U



Section 2: Rings in a canonic Extended Smallest Set of Smallest 
Rings (ESSSR) ring set - These bits test for the presence 
or count of the described chemical ring system.  An ESSSR 
ring is any ring which does not share three consecutive 
atoms with any other ring in the chemical structure.  For 
example, naphthalene has three ESSSR rings (two phenyl 
fragments and the 10-membered envelope), while biphenyl 
will yield a count of only two ESSSR rings.

Bit Position	Bit Substructure
115	>= 1 any ring size 3
116	>= 1 saturated or aromatic carbon-only ring size 3
117	>= 1 saturated or aromatic nitrogen-containing ring size 3
118	>= 1 saturated or aromatic heteroatom-containing ring size 3
119	>= 1 unsaturated non-aromatic carbon-only ring size 3
120	>= 1 unsaturated non-aromatic nitrogen-containing ring size 3
121	>= 1 unsaturated non-aromatic heteroatom-containing ring size 3
122	>= 2 any ring size 3
123	>= 2 saturated or aromatic carbon-only ring size 3
124	>= 2 saturated or aromatic nitrogen-containing ring size 3
125	>= 2 saturated or aromatic heteroatom-containing ring size 3
126	>= 2 unsaturated non-aromatic carbon-only ring size 3
127	>= 2 unsaturated non-aromatic nitrogen-containing ring size 3
128	>= 2 unsaturated non-aromatic heteroatom-containing ring size 3

129	>= 1 any ring size 4
130	>= 1 saturated or aromatic carbon-only ring size 4
131	>= 1 saturated or aromatic nitrogen-containing ring size 4
132	>= 1 saturated or aromatic heteroatom-containing ring size 4
133	>= 1 unsaturated non-aromatic carbon-only ring size 4
134	>= 1 unsaturated non-aromatic nitrogen-containing ring size 4
135	>= 1 unsaturated non-aromatic heteroatom-containing ring size 4
136	>= 2 any ring size 4
137	>= 2 saturated or aromatic carbon-only ring size 4
138	>= 2 saturated or aromatic nitrogen-containing ring size 4
139	>= 2 saturated or aromatic heteroatom-containing ring size 4
140	>= 2 unsaturated non-aromatic carbon-only ring size 4
141	>= 2 unsaturated non-aromatic nitrogen-containing ring size 4
142	>= 2 unsaturated non-aromatic heteroatom-containing ring size 4

143	>= 1 any ring size 5
144	>= 1 saturated or aromatic carbon-only ring size 5
145	>= 1 saturated or aromatic nitrogen-containing ring size 5
146	>= 1 saturated or aromatic heteroatom-containing ring size 5
147	>= 1 unsaturated non-aromatic carbon-only ring size 5
148	>= 1 unsaturated non-aromatic nitrogen-containing ring size 5
149	>= 1 unsaturated non-aromatic heteroatom-containing ring size 5
150	>= 2 any ring size 5
151	>= 2 saturated or aromatic carbon-only ring size 5
152	>= 2 saturated or aromatic nitrogen-containing ring size 5
153	>= 2 saturated or aromatic heteroatom-containing ring size 5
154	>= 2 unsaturated non-aromatic carbon-only ring size 5
155	>= 2 unsaturated non-aromatic nitrogen-containing ring size 5
156	>= 2 unsaturated non-aromatic heteroatom-containing ring size 5
157	>= 3 any ring size 5
158	>= 3 saturated or aromatic carbon-only ring size 5
159	>= 3 saturated or aromatic nitrogen-containing ring size 5
160	>= 3 saturated or aromatic heteroatom-containing ring size 5
161	>= 3 unsaturated non-aromatic carbon-only ring size 5
162	>= 3 unsaturated non-aromatic nitrogen-containing ring size 5
163	>= 3 unsaturated non-aromatic heteroatom-containing ring size 5
164	>= 4 any ring size 5
165	>= 4 saturated or aromatic carbon-only ring size 5
166	>= 4 saturated or aromatic nitrogen-containing ring size 5
167	>= 4 saturated or aromatic heteroatom-containing ring size 5
168	>= 4 unsaturated non-aromatic carbon-only ring size 5
169	>= 4 unsaturated non-aromatic nitrogen-containing ring size 5
170	>= 4 unsaturated non-aromatic heteroatom-containing ring size 5
171	>= 5 any ring size 5
172	>= 5 saturated or aromatic carbon-only ring size 5
173	>= 5 saturated or aromatic nitrogen-containing ring size 5
174	>= 5 saturated or aromatic heteroatom-containing ring size 5
175	>= 5 unsaturated non-aromatic carbon-only ring size 5
176	>= 5 unsaturated non-aromatic nitrogen-containing ring size 5
177	>= 5 unsaturated non-aromatic heteroatom-containing ring size 5

178	>= 1 any ring size 6
179	>= 1 saturated or aromatic carbon-only ring size 6
180	>= 1 saturated or aromatic nitrogen-containing ring size 6
181	>= 1 saturated or aromatic heteroatom-containing ring size 6
182	>= 1 unsaturated non-aromatic carbon-only ring size 6
183	>= 1 unsaturated non-aromatic nitrogen-containing ring size 6
184	>= 1 unsaturated non-aromatic heteroatom-containing ring size 6
185	>= 2 any ring size 6
186	>= 2 saturated or aromatic carbon-only ring size 6
187	>= 2 saturated or aromatic nitrogen-containing ring size 6
188	>= 2 saturated or aromatic heteroatom-containing ring size 6
189	>= 2 unsaturated non-aromatic carbon-only ring size 6
190	>= 2 unsaturated non-aromatic nitrogen-containing ring size 6
191	>= 2 unsaturated non-aromatic heteroatom-containing ring size 6
192	>= 3 any ring size 6
193	>= 3 saturated or aromatic carbon-only ring size 6
194	>= 3 saturated or aromatic nitrogen-containing ring size 6
195	>= 3 saturated or aromatic heteroatom-containing ring size 6
196	>= 3 unsaturated non-aromatic carbon-only ring size 6
197	>= 3 unsaturated non-aromatic nitrogen-containing ring size 6
198	>= 3 unsaturated non-aromatic heteroatom-containing ring size 6
199	>= 4 any ring size 6
200	>= 4 saturated or aromatic carbon-only ring size 6
201	>= 4 saturated or aromatic nitrogen-containing ring size 6
202	>= 4 saturated or aromatic heteroatom-containing ring size 6
203	>= 4 unsaturated non-aromatic carbon-only ring size 6
204	>= 4 unsaturated non-aromatic nitrogen-containing ring size 6
205	>= 4 unsaturated non-aromatic heteroatom-containing ring size 6
206	>= 5 any ring size 6
207	>= 5 saturated or aromatic carbon-only ring size 6
208	>= 5 saturated or aromatic nitrogen-containing ring size 6
209	>= 5 saturated or aromatic heteroatom-containing ring size 6
210	>= 5 unsaturated non-aromatic carbon-only ring size 6
211	>= 5 unsaturated non-aromatic nitrogen-containing ring size 6
212	>= 5 unsaturated non-aromatic heteroatom-containing ring size 6

213	>= 1 any ring size 7
214	>= 1 saturated or aromatic carbon-only ring size 7
215	>= 1 saturated or aromatic nitrogen-containing ring size 7
216	>= 1 saturated or aromatic heteroatom-containing ring size 7
217	>= 1 unsaturated non-aromatic carbon-only ring size 7
218	>= 1 unsaturated non-aromatic nitrogen-containing ring size 7
219	>= 1 unsaturated non-aromatic heteroatom-containing ring size 7
220	>= 2 any ring size 7
221	>= 2 saturated or aromatic carbon-only ring size 7
222	>= 2 saturated or aromatic nitrogen-containing ring size 7
223	>= 2 saturated or aromatic heteroatom-containing ring size 7
224	>= 2 unsaturated non-aromatic carbon-only ring size 7
225	>= 2 unsaturated non-aromatic nitrogen-containing ring size 7
226	>= 2 unsaturated non-aromatic heteroatom-containing ring size 7

227	>= 1 any ring size 8
228	>= 1 saturated or aromatic carbon-only ring size 8
229	>= 1 saturated or aromatic nitrogen-containing ring size 8
230	>= 1 saturated or aromatic heteroatom-containing ring size 8
231	>= 1 unsaturated non-aromatic carbon-only ring size 8
232	>= 1 unsaturated non-aromatic nitrogen-containing ring size 8
233	>= 1 unsaturated non-aromatic heteroatom-containing ring size 8
234	>= 2 any ring size 8
235	>= 2 saturated or aromatic carbon-only ring size 8
236	>= 2 saturated or aromatic nitrogen-containing ring size 8
237	>= 2 saturated or aromatic heteroatom-containing ring size 8
238	>= 2 unsaturated non-aromatic carbon-only ring size 8
239	>= 2 unsaturated non-aromatic nitrogen-containing ring size 8
240	>= 2 unsaturated non-aromatic heteroatom-containing ring size 8

241	>= 1 any ring size 9
242	>= 1 saturated or aromatic carbon-only ring size 9
243	>= 1 saturated or aromatic nitrogen-containing ring size 9
244	>= 1 saturated or aromatic heteroatom-containing ring size 9
245	>= 1 unsaturated non-aromatic carbon-only ring size 9
246	>= 1 unsaturated non-aromatic nitrogen-containing ring size 9
247	>= 1 unsaturated non-aromatic heteroatom-containing ring size 9

248	>= 1 any ring size 10
249	>= 1 saturated or aromatic carbon-only ring size 10
250	>= 1 saturated or aromatic nitrogen-containing ring size 10
251	>= 1 saturated or aromatic heteroatom-containing ring size 10
252	>= 1 unsaturated non-aromatic carbon-only ring size 10
253	>= 1 unsaturated non-aromatic nitrogen-containing ring size 10
254	>= 1 unsaturated non-aromatic heteroatom-containing ring size 10

255	>= 1 aromatic ring
256	>= 1 hetero-aromatic ring
257	>= 2 aromatic rings
258	>= 2 hetero-aromatic rings
259	>= 3 aromatic rings
260	>= 3 hetero-aromatic rings
261	>= 4 aromatic rings
262	>= 4 hetero-aromatic rings



Section 3: Simple atom pairs - These bits test for the presence of 
patterns of bonded atom pairs, regardless of bond order 
or count.

Bit Position	Bit Substructure
263	Li-H
264	Li-Li
265	Li-B
266	Li-C
267	Li-O
268	Li-F
269	Li-P
270	Li-S
271	Li-Cl
272	B-H
273	B-B
274	B-C
275	B-N
276	B-O
277	B-F
278	B-Si
279	B-P
280	B-S
281	B-Cl
282	B-Br
283	C-H
284	C-C
285	C-N
286	C-O
287	C-F
288	C-Na
289	C-Mg
290	C-Al
291	C-Si
292	C-P
293	C-S
294	C-Cl
295	C-As
296	C-Se
297	C-Br
298	C-I
299	N-H
300	N-N
301	N-O
302	N-F
303	N-Si
304	N-P
305	N-S
306	N-Cl
307	N-Br
308	O-H
309	O-O
310	O-Mg
311	O-Na
312	O-Al
313	O-Si
314	O-P
315	O-K
316	F-P
317	F-S
318	Al-H
319	Al-Cl
320	Si-H
321	Si-Si
322	Si-Cl
323	P-H
324	P-P
325	As-H
326	As-As



Section 4: Simple atom nearest neighbors - These bits test for the 
presence of atom nearest neighbor patterns, regardless of 
bond order (denoted by "~") or count, but where bond 
aromaticity (denoted by ":") is significant.

Bit Position	Bit Substructure
327	C(~Br)(~C)
328	C(~Br)(~C)(~C)
329	C(~Br)(~H)
330	C(~Br)(:C)
331	C(~Br)(:N)
332	C(~C)(~C)
333	C(~C)(~C)(~C)
334	C(~C)(~C)(~C)(~C)
335	C(~C)(~C)(~C)(~H)
336	C(~C)(~C)(~C)(~N)
337	C(~C)(~C)(~C)(~O)
338	C(~C)(~C)(~H)(~N)
339	C(~C)(~C)(~H)(~O)
340	C(~C)(~C)(~N)
341	C(~C)(~C)(~O)
342	C(~C)(~Cl)
343	C(~C)(~Cl)(~H)
344	C(~C)(~H)
345	C(~C)(~H)(~N)
346	C(~C)(~H)(~O)
347	C(~C)(~H)(~O)(~O)
348	C(~C)(~H)(~P)
349	C(~C)(~H)(~S)
350	C(~C)(~I)
351	C(~C)(~N)
352	C(~C)(~O)
353	C(~C)(~S)
354	C(~C)(~Si)
355	C(~C)(:C)
356	C(~C)(:C)(:C)
357	C(~C)(:C)(:N)
358	C(~C)(:N)
359	C(~C)(:N)(:N)
360	C(~Cl)(~Cl)
361	C(~Cl)(~H)
362	C(~Cl)(:C)
363	C(~F)(~F)
364	C(~F)(:C)
365	C(~H)(~N)
366	C(~H)(~O)
367	C(~H)(~O)(~O)
368	C(~H)(~S)
369	C(~H)(~Si)
370	C(~H)(:C)
371	C(~H)(:C)(:C)
372	C(~H)(:C)(:N)
373	C(~H)(:N)
374	C(~H)(~H)(~H)
375	C(~N)(~N)
376	C(~N)(:C)
377	C(~N)(:C)(:C)
378	C(~N)(:C)(:N)
379	C(~N)(:N)
380	C(~O)(~O)
381	C(~O)(:C)
382	C(~O)(:C)(:C)
383	C(~S)(:C)
384	C(:C)(:C)
385	C(:C)(:C)(:C)
386	C(:C)(:C)(:N)
387	C(:C)(:N)
388	C(:C)(:N)(:N)
389	C(:N)(:N)
390	N(~C)(~C)
391	N(~C)(~C)(~C)
392	N(~C)(~C)(~H)
393	N(~C)(~H)
394	N(~C)(~H)(~N)
395	N(~C)(~O)
396	N(~C)(:C)
397	N(~C)(:C)(:C)
398	N(~H)(~N)
399	N(~H)(:C)
400	N(~H)(:C)(:C)
401	N(~O)(~O)
402	N(~O)(:O)
403	N(:C)(:C)
404	N(:C)(:C)(:C)
405	O(~C)(~C)
406	O(~C)(~H)
407	O(~C)(~P)
408	O(~H)(~S)
409	O(:C)(:C)
410	P(~C)(~C)
411	P(~O)(~O)
412	S(~C)(~C)
413	S(~C)(~H)
414	S(~C)(~O)
415	Si(~C)(~C)



Section 5: Detailed atom neighborhoods - These bits test for the 
presence of detailed atom neighborhood patterns, 
regardless of count, but where bond orders are specific, 
bond aromaticity matches both single and double bonds, 
and where "-", "=", and "#" matches a single bond, double 
bond, and triple bond order, respectively.

Bit Position	Bit Substructure
416	C=C
417	C#C
418	C=N
419	C#N
420	C=O
421	C=S
422	N=N
423	N=O
424	N=P
425	P=O
426	P=P
427	C(#C)(-C)
428	C(#C)(-H)
429	C(#N)(-C)
430	C(-C)(-C)(=C)
431	C(-C)(-C)(=N)
432	C(-C)(-C)(=O)
433	C(-C)(-Cl)(=O)
434	C(-C)(-H)(=C)
435	C(-C)(-H)(=N)
436	C(-C)(-H)(=O)
437	C(-C)(-N)(=C)
438	C(-C)(-N)(=N)
439	C(-C)(-N)(=O)
440	C(-C)(-O)(=O)
441	C(-C)(=C)
442	C(-C)(=N)
443	C(-C)(=O)
444	C(-Cl)(=O)
445	C(-H)(-N)(=C)
446	C(-H)(=C)
447	C(-H)(=N)
448	C(-H)(=O)
449	C(-N)(=C)
450	C(-N)(=N)
451	C(-N)(=O)
452	C(-O)(=O)
453	N(-C)(=C)
454	N(-C)(=O)
455	N(-O)(=O)
456	P(-O)(=O)
457	S(-C)(=O)
458	S(-O)(=O)
459	S(=O)(=O)



Section 6: Simple SMARTS patterns - These bits test for the 
presence of simple SMARTS patterns, regardless of count, 
but where bond orders are specific and bond aromaticity 
matches both single and double bonds.

Bit Position	Bit Substructure
460	C-C-C#C
461	O-C-C=N
462	O-C-C=O
463	N:C-S-[#1]
464	N-C-C=C
465	O=S-C-C
466	N#C-C=C
467	C=N-N-C
468	O=S-C-N
469	S-S-C:C
470	C:C-C=C
471	S:C:C:C
472	C:N:C-C
473	S-C:N:C
474	S:C:C:N
475	S-C=N-C
476	C-O-C=C
477	N-N-C:C
478	S-C=N-[#1]
479	S-C-S-C
480	C:S:C-C
481	O-S-C:C
482	C:N-C:C
483	N-S-C:C
484	N-C:N:C
485	N:C:C:N
486	N-C:N:N
487	N-C=N-C
488	N-C=N-[#1]
489	N-C-S-C
490	C-C-C=C
491	C-N:C-[#1]
492	N-C:O:C
493	O=C-C:C
494	O=C-C:N
495	C-N-C:C
496	N:N-C-[#1]
497	O-C:C:N
498	O-C=C-C
499	N-C:C:N
500	C-S-C:C
501	Cl-C:C-C
502	N-C=C-[#1]
503	Cl-C:C-[#1]
504	N:C:N-C
505	Cl-C:C-O
506	C-C:N:C
507	C-C-S-C
508	S=C-N-C
509	Br-C:C-C
510	[#1]-N-N-[#1]
511	S=C-N-[#1]
512	C-[As]-O-[#1]
513	S:C:C-[#1]
514	O-N-C-C
515	N-N-C-C
516	[#1]-C=C-[#1]
517	N-N-C-N
518	O=C-N-N
519	N=C-N-C
520	C=C-C:C
521	C:N-C-[#1]
522	C-N-N-[#1]
523	N:C:C-C
524	C-C=C-C
525	[As]-C:C-[#1]
526	Cl-C:C-Cl
527	C:C:N-[#1]
528	[#1]-N-C-[#1]
529	Cl-C-C-Cl
530	N:C-C:C
531	S-C:C-C
532	S-C:C-[#1]
533	S-C:C-N
534	S-C:C-O
535	O=C-C-C
536	O=C-C-N
537	O=C-C-O
538	N=C-C-C
539	N=C-C-[#1]
540	C-N-C-[#1]
541	O-C:C-C
542	O-C:C-[#1]
543	O-C:C-N
544	O-C:C-O
545	N-C:C-C
546	N-C:C-[#1]
547	N-C:C-N
548	O-C-C:C
549	N-C-C:C
550	Cl-C-C-C
551	Cl-C-C-O
552	C:C-C:C
553	O=C-C=C
554	Br-C-C-C
555	N=C-C=C
556	C=C-C-C
557	N:C-O-[#1]
558	O=N-C:C
559	O-C-N-[#1]
560	N-C-N-C
561	Cl-C-C=O
562	Br-C-C=O
563	O-C-O-C
564	C=C-C=C
565	C:C-O-C
566	O-C-C-N
567	O-C-C-O
568	N#C-C-C
569	N-C-C-N
570	C:C-C-C
571	[#1]-C-O-[#1]
572	N:C:N:C
573	O-C-C=C
574	O-C-C:C-C
575	O-C-C:C-O
576	N=C-C:C-[#1]
577	C:C-N-C:C
578	C-C:C-C:C
579	O=C-C-C-C
580	O=C-C-C-N
581	O=C-C-C-O
582	C-C-C-C-C
583	Cl-C:C-O-C
584	C:C-C=C-C
585	C-C:C-N-C
586	C-S-C-C-C
587	N-C:C-O-[#1]
588	O=C-C-C=O
589	C-C:C-O-C
590	C-C:C-O-[#1]
591	Cl-C-C-C-C
592	N-C-C-C-C
593	N-C-C-C-N
594	C-O-C-C=C
595	C:C-C-C-C
596	N=C-N-C-C
597	O=C-C-C:C
598	Cl-C:C:C-C
599	[#1]-C-C=C-[#1]
600	N-C:C:C-C
601	N-C:C:C-N
602	O=C-C-N-C
603	C-C:C:C-C
604	C-O-C-C:C
605	O=C-C-O-C
606	O-C:C-C-C
607	N-C-C-C:C
608	C-C-C-C:C
609	Cl-C-C-N-C
610	C-O-C-O-C
611	N-C-C-N-C
612	N-C-O-C-C
613	C-N-C-C-C
614	C-C-O-C-C
615	N-C-C-O-C
616	C:C:N:N:C
617	C-C-C-O-[#1]
618	C:C-C-C:C
619	O-C-C=C-C
620	C:C-O-C-C
621	N-C:C:C:N
622	O=C-O-C:C
623	O=C-C:C-C
624	O=C-C:C-N
625	O=C-C:C-O
626	C-O-C:C-C
627	O=[As]-C:C:C
628	C-N-C-C:C
629	S-C:C:C-N
630	O-C:C-O-C
631	O-C:C-O-[#1]
632	C-C-O-C:C
633	N-C-C:C-C
634	C-C-C:C-C
635	N-N-C-N-[#1]
636	C-N-C-N-C
637	O-C-C-C-C
638	O-C-C-C-N
639	O-C-C-C-O
640	C=C-C-C-C
641	O-C-C-C=C
642	O-C-C-C=O
643	[#1]-C-C-N-[#1]
644	C-C=N-N-C
645	O=C-N-C-C
646	O=C-N-C-[#1]
647	O=C-N-C-N
648	O=N-C:C-N
649	O=N-C:C-O
650	O=C-N-C=O
651	O-C:C:C-C
652	O-C:C:C-N
653	O-C:C:C-O
654	N-C-N-C-C
655	O-C-C-C:C
656	C-C-N-C-C
657	C-N-C:C-C
658	C-C-S-C-C
659	O-C-C-N-C
660	C-C=C-C-C
661	O-C-O-C-C
662	O-C-C-O-C
663	O-C-C-O-[#1]
664	C-C=C-C=C
665	N-C:C-C-C
666	C=C-C-O-C
667	C=C-C-O-[#1]
668	C-C:C-C-C
669	Cl-C:C-C=O
670	Br-C:C:C-C
671	O=C-C=C-C
672	O=C-C=C-[#1]
673	O=C-C=C-N
674	N-C-N-C:C
675	Br-C-C-C:C
676	N#C-C-C-C
677	C-C=C-C:C
678	C-C-C=C-C
679	C-C-C-C-C-C
680	O-C-C-C-C-C
681	O-C-C-C-C-O
682	O-C-C-C-C-N
683	N-C-C-C-C-C
684	O=C-C-C-C-C
685	O=C-C-C-C-N
686	O=C-C-C-C-O
687	O=C-C-C-C=O
688	C-C-C-C-C-C-C
689	O-C-C-C-C-C-C
690	O-C-C-C-C-C-O
691	O-C-C-C-C-C-N
692	O=C-C-C-C-C-C
693	O=C-C-C-C-C-O
694	O=C-C-C-C-C=O
695	O=C-C-C-C-C-N
696	C-C-C-C-C-C-C-C
697	C-C-C-C-C-C(C)-C
698	O-C-C-C-C-C-C-C
699	O-C-C-C-C-C(C)-C
700	O-C-C-C-C-C-O-C
701	O-C-C-C-C-C(O)-C
702	O-C-C-C-C-C-N-C
703	O-C-C-C-C-C(N)-C
704	O=C-C-C-C-C-C-C
705	O=C-C-C-C-C(O)-C
706	O=C-C-C-C-C(=O)-C
707	O=C-C-C-C-C(N)-C
708	C-C(C)-C-C
709	C-C(C)-C-C-C
710	C-C-C(C)-C-C
711	C-C(C)(C)-C-C
712	C-C(C)-C(C)-C



Section 7: Complex SMARTS patterns - These bits test for the 
presence of complex SMARTS patterns, regardless of count, 
but where bond orders and bond aromaticity are specific.

Bit Position	Bit Substructure
713	Cc1ccc(C)cc1
714	Cc1ccc(O)cc1
715	Cc1ccc(S)cc1
716	Cc1ccc(N)cc1
717	Cc1ccc(Cl)cc1
718	Cc1ccc(Br)cc1
719	Oc1ccc(O)cc1
720	Oc1ccc(S)cc1
721	Oc1ccc(N)cc1
722	Oc1ccc(Cl)cc1
723	Oc1ccc(Br)cc1
724	Sc1ccc(S)cc1
725	Sc1ccc(N)cc1
726	Sc1ccc(Cl)cc1
727	Sc1ccc(Br)cc1
728	Nc1ccc(N)cc1
729	Nc1ccc(Cl)cc1
730	Nc1ccc(Br)cc1
731	Clc1ccc(Cl)cc1
732	Clc1ccc(Br)cc1
733	Brc1ccc(Br)cc1
734	Cc1cc(C)ccc1
735	Cc1cc(O)ccc1
736	Cc1cc(S)ccc1
737	Cc1cc(N)ccc1
738	Cc1cc(Cl)ccc1
739	Cc1cc(Br)ccc1
740	Oc1cc(O)ccc1
741	Oc1cc(S)ccc1
742	Oc1cc(N)ccc1
743	Oc1cc(Cl)ccc1
744	Oc1cc(Br)ccc1
745	Sc1cc(S)ccc1
746	Sc1cc(N)ccc1
747	Sc1cc(Cl)ccc1
748	Sc1cc(Br)ccc1
749	Nc1cc(N)ccc1
750	Nc1cc(Cl)ccc1
751	Nc1cc(Br)ccc1
752	Clc1cc(Cl)ccc1
753	Clc1cc(Br)ccc1
754	Brc1cc(Br)ccc1
755	Cc1c(C)cccc1
756	Cc1c(O)cccc1
757	Cc1c(S)cccc1
758	Cc1c(N)cccc1
759	Cc1c(Cl)cccc1
760	Cc1c(Br)cccc1
761	Oc1c(O)cccc1
762	Oc1c(S)cccc1
763	Oc1c(N)cccc1
764	Oc1c(Cl)cccc1
765	Oc1c(Br)cccc1
766	Sc1c(S)cccc1
767	Sc1c(N)cccc1
768	Sc1c(Cl)cccc1
769	Sc1c(Br)cccc1
770	Nc1c(N)cccc1
771	Nc1c(Cl)cccc1
772	Nc1c(Br)cccc1
773	Clc1c(Cl)cccc1
774	Clc1c(Br)cccc1
775	Brc1c(Br)cccc1
776	CC1CCC(C)CC1
777	CC1CCC(O)CC1
778	CC1CCC(S)CC1
779	CC1CCC(N)CC1
780	CC1CCC(Cl)CC1
781	CC1CCC(Br)CC1
782	OC1CCC(O)CC1
783	OC1CCC(S)CC1
784	OC1CCC(N)CC1
785	OC1CCC(Cl)CC1
786	OC1CCC(Br)CC1
787	SC1CCC(S)CC1
788	SC1CCC(N)CC1
789	SC1CCC(Cl)CC1
790	SC1CCC(Br)CC1
791	NC1CCC(N)CC1
792	NC1CCC(Cl)CC1
793	NC1CCC(Br)CC1
794	ClC1CCC(Cl)CC1
795	ClC1CCC(Br)CC1
796	BrC1CCC(Br)CC1
797	CC1CC(C)CCC1
798	CC1CC(O)CCC1
799	CC1CC(S)CCC1
800	CC1CC(N)CCC1
801	CC1CC(Cl)CCC1
802	CC1CC(Br)CCC1
803	OC1CC(O)CCC1
804	OC1CC(S)CCC1
805	OC1CC(N)CCC1
806	OC1CC(Cl)CCC1
807	OC1CC(Br)CCC1
808	SC1CC(S)CCC1
809	SC1CC(N)CCC1
810	SC1CC(Cl)CCC1
811	SC1CC(Br)CCC1
812	NC1CC(N)CCC1
813	NC1CC(Cl)CCC1
814	NC1CC(Br)CCC1
815	ClC1CC(Cl)CCC1
816	ClC1CC(Br)CCC1
817	BrC1CC(Br)CCC1
818	CC1C(C)CCCC1
819	CC1C(O)CCCC1
820	CC1C(S)CCCC1
821	CC1C(N)CCCC1
822	CC1C(Cl)CCCC1
823	CC1C(Br)CCCC1
824	OC1C(O)CCCC1
825	OC1C(S)CCCC1
826	OC1C(N)CCCC1
827	OC1C(Cl)CCCC1
828	OC1C(Br)CCCC1
829	SC1C(S)CCCC1
830	SC1C(N)CCCC1
831	SC1C(Cl)CCCC1
832	SC1C(Br)CCCC1
833	NC1C(N)CCCC1
834	NC1C(Cl)CCCC1
835	NC1C(Br)CCCC1
836	ClC1C(Cl)CCCC1
837	ClC1C(Br)CCCC1
838	BrC1C(Br)CCCC1
839	CC1CC(C)CC1
840	CC1CC(O)CC1
841	CC1CC(S)CC1
842	CC1CC(N)CC1
843	CC1CC(Cl)CC1
844	CC1CC(Br)CC1
845	OC1CC(O)CC1
846	OC1CC(S)CC1
847	OC1CC(N)CC1
848	OC1CC(Cl)CC1
849	OC1CC(Br)CC1
850	SC1CC(S)CC1
851	SC1CC(N)CC1
852	SC1CC(Cl)CC1
853	SC1CC(Br)CC1
854	NC1CC(N)CC1
855	NC1CC(Cl)CC1
856	NC1CC(Br)CC1
857	ClC1CC(Cl)CC1
858	ClC1CC(Br)CC1
859	BrC1CC(Br)CC1
860	CC1C(C)CCC1
861	CC1C(O)CCC1
862	CC1C(S)CCC1
863	CC1C(N)CCC1
864	CC1C(Cl)CCC1
865	CC1C(Br)CCC1
866	OC1C(O)CCC1
867	OC1C(S)CCC1
868	OC1C(N)CCC1
869	OC1C(Cl)CCC1
870	OC1C(Br)CCC1
871	SC1C(S)CCC1
872	SC1C(N)CCC1
873	SC1C(Cl)CCC1
874	SC1C(Br)CCC1
875	NC1C(N)CCC1
876	NC1C(Cl)CC1
877	NC1C(Br)CCC1
878	ClC1C(Cl)CCC1
879	ClC1C(Br)CCC1
880	BrC1C(Br)CCC1



Decoding PubChem Fingerprints

PubChem fingerprints are currently 881 bits in length.  Binary data is stored in one 
byte increments.  The fingerprint is, therefore, 111 bytes in length (888 bits), 
which includes padding of seven bits at the end to complete the last byte.  A four-
byte prefix, containing the bit length of the fingerprint (881 bits), increases the 
stored PubChem fingerprint size to 115 bytes (920 bits).

When PubChem fingerprints are encoded in base64 format, the base64-encoded 
fingerprints are 156 bytes in length.  The last two bytes are padding so that the 
base64 length is divisible by four (156 bytes - 2 bytes = 154 bytes).  Each base64 
byte encodes six binary bits (154 bytes * 6 bits/byte = 924 bits).  The last four 
bits are padding to complete the last base64 byte (924 bits - 4 bits = 920 bits).  
The resulting 920 binary bits (115 bytes) are described in the previous paragraph.



Document Version History

V1.3 - 2009May01 - Updated introduction to describe how to identify 
the PubChem Substructure Fingerprint property in 
a PubChem Compound record.
V1.2 - 2007Aug30 - Added section on decoding PubChem fingerprints.
V1.1 - 2007Aug06 - Corrected and expanded documentation of bits 
with SMARTS patterns used.
V1.0 - 2005Dec02 - Initial release.

*/

// clang-format on

#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

using std::cerr;
using std::endl;

const char *prog_name = nullptr;

static int verbose = 0;

static int molecules_read = 0;
static int tdts_read = 0;

static int function_as_filter = 0;

static int write_as_array = 0;
static int write_zero_bits = 0;

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString fingerprint_tag("FPPB<");

static Accumulator_Int<int> bits_hit;

static IWString stream_for_bit_stem("PBFPB");

static IWString_and_File_Descriptor *write_bit_file = nullptr;

/*
  The query groups are differentiated as in the documentation
  because some need different treatment of hydrogens.
  Right now, only section 6 needs hydrogens
*/

static resizable_array_p<Substructure_Hit_Statistics> section4;
static resizable_array_p<Substructure_Hit_Statistics> section5;
static resizable_array_p<Substructure_Hit_Statistics> section6;
static resizable_array_p<Substructure_Hit_Statistics> section7;

static int external_queries_present = 0;

static int flush_after_every_molecule = 0;

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
// clang-format off
  cerr << "Computes Pubchem Fingerprints\n";
  cerr << "  -a            write as an array\n";
  cerr << "  -J <tag>      tag for fingerprints\n";
  cerr << "  -q ...        external query specifications\n";
  cerr << "  -X ...        miscellaneous options, enter '-X help'\n";
  cerr << "  -f            function as a TDT filter\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -X ...        other options, enter '-X help' for info\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -E ...        standard element specifications\n";
  cerr << "  -A ...        standard aromaticity specifications\n";
  cerr << "  -v            verbose output\n";
// clang-format on

  exit(rc);
}

static void
preprocess(Molecule &m) {
  if (reduce_to_largest_fragment)
    m.reduce_to_largest_fragment();

  if (chemical_standardisation.active())
    chemical_standardisation.process(m);

  m.remove_all_chiral_centres();
  m.revert_all_directional_bonds_to_non_directional();

  return;
}

static void
do_queries(Molecule_to_Match &target,
           resizable_array_p<Substructure_Hit_Statistics> &q,
           int *fp) {
  const int nq = q.number_elements();

  for (int i = 0; i < nq; i++) {
    Substructure_Hit_Statistics *qi = q[i];

    if (0 == qi->comment().length())
      continue;

    if (q[i]->substructure_search(target))
      fp[i] = 1;
  }

#ifdef COUNT_QUERIES_MATCHING
  int rc = 0;
  for (int i = 0; i < nq; ++i) {
    if (fp[i])
      rc++;
  }
  cerr << " hit " << rc << " of " << nq << " queries hit\n";
#endif

  return;
}

static int
do_external_queries(Molecule &m, int *fp) {
  Molecule_to_Match target1(&m);

  do_queries(target1, section4, fp + 327);
  do_queries(target1, section5, fp + 416);
  do_queries(target1, section7, fp + 713);

  m.make_implicit_hydrogens_explicit();
  m.compute_aromaticity_if_needed();

  Molecule_to_Match target2(&m);

  do_queries(target2, section6, fp + 460);

  return 1;
}

/*
  We need an index into where the bit for various pairs is stored.
*/

#define SIZE_OF_PAIR_XREF                                                      \
  ((HIGHEST_ATOMIC_NUMBER + 1) * (HIGHEST_ATOMIC_NUMBER + 1))

static int pair_xref[SIZE_OF_PAIR_XREF];

static void
initialise_pair_xref() {
  set_vector(pair_xref, SIZE_OF_PAIR_XREF, -1);

  pair_xref[3 * HIGHEST_ATOMIC_NUMBER + 3] = 264;
  pair_xref[3 * HIGHEST_ATOMIC_NUMBER + 3] = 264;
  pair_xref[3 * HIGHEST_ATOMIC_NUMBER + 5] = 265;
  pair_xref[5 * HIGHEST_ATOMIC_NUMBER + 3] = 265;
  pair_xref[3 * HIGHEST_ATOMIC_NUMBER + 6] = 266;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 3] = 266;
  pair_xref[3 * HIGHEST_ATOMIC_NUMBER + 9] = 267;
  pair_xref[8 * HIGHEST_ATOMIC_NUMBER + 3] = 267;
  pair_xref[3 * HIGHEST_ATOMIC_NUMBER + 8] = 268;
  pair_xref[9 * HIGHEST_ATOMIC_NUMBER + 3] = 268;
  pair_xref[3 * HIGHEST_ATOMIC_NUMBER + 15] = 269;
  pair_xref[15 * HIGHEST_ATOMIC_NUMBER + 3] = 269;
  pair_xref[3 * HIGHEST_ATOMIC_NUMBER + 16] = 270;
  pair_xref[16 * HIGHEST_ATOMIC_NUMBER + 3] = 270;
  pair_xref[3 * HIGHEST_ATOMIC_NUMBER + 17] = 271;
  pair_xref[17 * HIGHEST_ATOMIC_NUMBER + 3] = 271;
  pair_xref[5 * HIGHEST_ATOMIC_NUMBER + 5] = 273;
  pair_xref[5 * HIGHEST_ATOMIC_NUMBER + 5] = 273;
  pair_xref[5 * HIGHEST_ATOMIC_NUMBER + 6] = 274;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 5] = 274;
  pair_xref[5 * HIGHEST_ATOMIC_NUMBER + 7] = 275;
  pair_xref[7 * HIGHEST_ATOMIC_NUMBER + 5] = 275;
  pair_xref[5 * HIGHEST_ATOMIC_NUMBER + 8] = 276;
  pair_xref[8 * HIGHEST_ATOMIC_NUMBER + 5] = 276;
  pair_xref[5 * HIGHEST_ATOMIC_NUMBER + 9] = 277;
  pair_xref[9 * HIGHEST_ATOMIC_NUMBER + 5] = 277;
  pair_xref[5 * HIGHEST_ATOMIC_NUMBER + 14] = 278;
  pair_xref[14 * HIGHEST_ATOMIC_NUMBER + 5] = 278;
  pair_xref[5 * HIGHEST_ATOMIC_NUMBER + 15] = 279;
  pair_xref[15 * HIGHEST_ATOMIC_NUMBER + 5] = 279;
  pair_xref[5 * HIGHEST_ATOMIC_NUMBER + 16] = 280;
  pair_xref[16 * HIGHEST_ATOMIC_NUMBER + 5] = 280;
  pair_xref[5 * HIGHEST_ATOMIC_NUMBER + 17] = 281;
  pair_xref[17 * HIGHEST_ATOMIC_NUMBER + 5] = 281;
  pair_xref[5 * HIGHEST_ATOMIC_NUMBER + 35] = 282;
  pair_xref[35 * HIGHEST_ATOMIC_NUMBER + 5] = 282;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 6] = 284;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 6] = 284;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 7] = 285;
  pair_xref[7 * HIGHEST_ATOMIC_NUMBER + 6] = 285;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 8] = 286;
  pair_xref[8 * HIGHEST_ATOMIC_NUMBER + 6] = 286;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 9] = 287;
  pair_xref[9 * HIGHEST_ATOMIC_NUMBER + 6] = 287;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 11] = 288;
  pair_xref[11 * HIGHEST_ATOMIC_NUMBER + 6] = 288;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 12] = 289;
  pair_xref[12 * HIGHEST_ATOMIC_NUMBER + 6] = 289;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 13] = 290;
  pair_xref[13 * HIGHEST_ATOMIC_NUMBER + 6] = 290;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 14] = 291;
  pair_xref[14 * HIGHEST_ATOMIC_NUMBER + 6] = 291;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 15] = 292;
  pair_xref[15 * HIGHEST_ATOMIC_NUMBER + 6] = 292;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 16] = 293;
  pair_xref[16 * HIGHEST_ATOMIC_NUMBER + 6] = 293;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 17] = 294;
  pair_xref[17 * HIGHEST_ATOMIC_NUMBER + 6] = 294;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 33] = 295;
  pair_xref[33 * HIGHEST_ATOMIC_NUMBER + 6] = 295;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 34] = 296;
  pair_xref[34 * HIGHEST_ATOMIC_NUMBER + 6] = 296;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 35] = 297;
  pair_xref[35 * HIGHEST_ATOMIC_NUMBER + 6] = 297;
  pair_xref[6 * HIGHEST_ATOMIC_NUMBER + 53] = 298;
  pair_xref[53 * HIGHEST_ATOMIC_NUMBER + 6] = 298;
  pair_xref[7 * HIGHEST_ATOMIC_NUMBER + 7] = 300;
  pair_xref[7 * HIGHEST_ATOMIC_NUMBER + 7] = 300;
  pair_xref[7 * HIGHEST_ATOMIC_NUMBER + 8] = 301;
  pair_xref[8 * HIGHEST_ATOMIC_NUMBER + 7] = 301;
  pair_xref[7 * HIGHEST_ATOMIC_NUMBER + 9] = 302;
  pair_xref[9 * HIGHEST_ATOMIC_NUMBER + 7] = 302;
  pair_xref[7 * HIGHEST_ATOMIC_NUMBER + 14] = 303;
  pair_xref[14 * HIGHEST_ATOMIC_NUMBER + 7] = 303;
  pair_xref[7 * HIGHEST_ATOMIC_NUMBER + 15] = 304;
  pair_xref[15 * HIGHEST_ATOMIC_NUMBER + 7] = 304;
  pair_xref[7 * HIGHEST_ATOMIC_NUMBER + 16] = 305;
  pair_xref[16 * HIGHEST_ATOMIC_NUMBER + 7] = 305;
  pair_xref[7 * HIGHEST_ATOMIC_NUMBER + 17] = 306;
  pair_xref[17 * HIGHEST_ATOMIC_NUMBER + 7] = 306;
  pair_xref[7 * HIGHEST_ATOMIC_NUMBER + 35] = 307;
  pair_xref[35 * HIGHEST_ATOMIC_NUMBER + 7] = 307;
  pair_xref[8 * HIGHEST_ATOMIC_NUMBER + 8] = 309;
  pair_xref[8 * HIGHEST_ATOMIC_NUMBER + 8] = 309;
  pair_xref[8 * HIGHEST_ATOMIC_NUMBER + 12] = 310;
  pair_xref[12 * HIGHEST_ATOMIC_NUMBER + 8] = 310;
  pair_xref[8 * HIGHEST_ATOMIC_NUMBER + 11] = 311;
  pair_xref[11 * HIGHEST_ATOMIC_NUMBER + 8] = 311;
  pair_xref[8 * HIGHEST_ATOMIC_NUMBER + 13] = 312;
  pair_xref[13 * HIGHEST_ATOMIC_NUMBER + 8] = 312;
  pair_xref[8 * HIGHEST_ATOMIC_NUMBER + 14] = 313;
  pair_xref[14 * HIGHEST_ATOMIC_NUMBER + 8] = 313;
  pair_xref[8 * HIGHEST_ATOMIC_NUMBER + 15] = 314;
  pair_xref[15 * HIGHEST_ATOMIC_NUMBER + 8] = 314;
  pair_xref[9 * HIGHEST_ATOMIC_NUMBER + 15] = 316;
  pair_xref[15 * HIGHEST_ATOMIC_NUMBER + 9] = 316;
  pair_xref[9 * HIGHEST_ATOMIC_NUMBER + 16] = 317;
  pair_xref[16 * HIGHEST_ATOMIC_NUMBER + 9] = 317;
  pair_xref[13 * HIGHEST_ATOMIC_NUMBER + 17] = 319;
  pair_xref[17 * HIGHEST_ATOMIC_NUMBER + 13] = 319;
  pair_xref[14 * HIGHEST_ATOMIC_NUMBER + 14] = 321;
  pair_xref[14 * HIGHEST_ATOMIC_NUMBER + 14] = 321;
  pair_xref[14 * HIGHEST_ATOMIC_NUMBER + 17] = 322;
  pair_xref[17 * HIGHEST_ATOMIC_NUMBER + 14] = 322;
  pair_xref[15 * HIGHEST_ATOMIC_NUMBER + 15] = 324;
  pair_xref[15 * HIGHEST_ATOMIC_NUMBER + 15] = 324;
  pair_xref[33 * HIGHEST_ATOMIC_NUMBER + 33] = 326;
  pair_xref[33 * HIGHEST_ATOMIC_NUMBER + 33] = 326;

  return;
}

#define LI_PAIR_H 263
#define B_PAIR_H 272
#define C_PAIR_H 283
#define N_PAIR_H 299
#define O_PAIR_H 308
#define AL_PAIR_H 318
#define SI_PAIR_H 320
#define P_PAIR_H 323

// CACTVS defines Si-H but not S-H. I'll do S-H instead

#define S_PAIR_H 320
#define P_PAIR_H 323

static void
atom_pairs(Molecule &m, const atomic_number_t *z, int *fp) {
  int nb = m.nedges();

  for (int i = 0; i < nb; i++) {
    const Bond *b = m.bondi(i);

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    const atomic_number_t z1 = z[a1];
    const atomic_number_t z2 = z[a2];

    int ndx = HIGHEST_ATOMIC_NUMBER * z1 + z2;

    if (pair_xref[ndx] >= 0)
      fp[pair_xref[ndx]]++;
  }

  // Since we do not have explicit hydrogens, count them here

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    const int h = m.hcount(i);
    if (0 == h)
      continue;

    const atomic_number_t zi = z[i];

    if (6 == zi)
      fp[C_PAIR_H] += h;
    else if (7 == zi)
      fp[N_PAIR_H] += h;
    else if (8 == zi)
      fp[O_PAIR_H] += h;
    else if (16 == zi)
      fp[S_PAIR_H] += h;
    else if (15 == zi)
      fp[P_PAIR_H] += h;
    else if (3 == zi)
      fp[LI_PAIR_H] += h;
    else if (5 == zi)
      fp[B_PAIR_H] += h;
    else if (14 == zi)
      fp[SI_PAIR_H] += h;
    else if (15 == zi)
      fp[P_PAIR_H] += h;
  }
}

static int
contains_heteroatom(const Molecule &m, const atomic_number_t *z,
                               const Ring &r) {
  int n = r.number_elements();

  for (int i = 0; i < n; i++) {
    atom_number_t j = r[i];

    if (6 != z[j])
      return 1;
  }

  return 0;
}

static void
update_count_bits(const int *ecount, int *fp) {
#ifdef DEBUG_UPDATE_COUNT_BITS
  cerr << "Ecount C " << ecount[6] << endl;
  cerr << "Ecount N " << ecount[7] << endl;
  cerr << "Ecount O " << ecount[8] << endl;
  cerr << "Ecount H " << ecount[1] << endl;
#endif

  if (ecount[1] >= 4)
    fp[0] = 1;
  if (ecount[1] >= 8)
    fp[1] = 1;
  if (ecount[1] >= 16)
    fp[2] = 1;
  if (ecount[1] >= 32)
    fp[3] = 1;

  if (ecount[5] >= 1)
    fp[6] = 1;
  if (ecount[5] >= 2)
    fp[7] = 1;
  if (ecount[5] >= 4)
    fp[8] = 1;

  if (ecount[6] >= 2)
    fp[9] = 1;
  if (ecount[6] >= 4)
    fp[10] = 1;
  if (ecount[6] >= 8)
    fp[11] = 1;
  if (ecount[6] >= 16)
    fp[12] = 1;
  if (ecount[6] >= 32)
    fp[13] = 1;

  if (ecount[7] >= 1)
    fp[14] = 1;
  if (ecount[7] >= 2)
    fp[15] = 1;
  if (ecount[7] >= 4)
    fp[16] = 1;
  if (ecount[7] >= 8)
    fp[17] = 1;

  if (ecount[8] >= 1)
    fp[18] = 1;
  if (ecount[8] >= 2)
    fp[19] = 1;
  if (ecount[8] >= 4)
    fp[20] = 1;
  if (ecount[8] >= 8)
    fp[21] = 1;
  if (ecount[8] >= 16)
    fp[22] = 1;

  if (ecount[9] >= 1)
    fp[23] = 1;
  if (ecount[9] >= 2)
    fp[24] = 1;
  if (ecount[9] >= 4)
    fp[25] = 1;

  if (ecount[11] >= 1) // Na
    fp[26] = 1;
  if (ecount[11] >= 2)
    fp[27] = 1;

  if (ecount[14] >= 1)
    fp[28] = 1;
  if (ecount[14] >= 2)
    fp[29] = 1;

  if (ecount[15] >= 1)
    fp[30] = 1;
  if (ecount[15] >= 2)
    fp[31] = 1;
  if (ecount[15] >= 4)
    fp[32] = 1;

  if (ecount[16] >= 1)
    fp[33] = 1;
  if (ecount[16] >= 2)
    fp[34] = 1;
  if (ecount[16] >= 4)
    fp[35] = 1;
  if (ecount[16] >= 8)
    fp[36] = 1;

  if (ecount[17] >= 1)
    fp[37] = 1;
  if (ecount[17] >= 2)
    fp[38] = 1;
  if (ecount[17] >= 4)
    fp[39] = 1;
  if (ecount[17] >= 8)
    fp[40] = 1;

  if (ecount[20] >= 1) // Ca
    fp[52] = 1;

  if (ecount[33] >= 1)
    fp[65] = 1;

  if (ecount[34] >= 1)
    fp[66] = 1;

  if (ecount[35] >= 1)
    fp[43] = 1;
  if (ecount[35] >= 2)
    fp[44] = 1;
  if (ecount[35] >= 4)
    fp[45] = 1;

  if (ecount[53] >= 1)
    fp[46] = 1;
  if (ecount[53] >= 2)
    fp[47] = 1;
  if (ecount[53] >= 4)
    fp[48] = 1;

  return;
}

class Ring_Statistics {
  private:
    int _saturated_or_aromatic_carbon_only;
    int _saturated_or_aromatic_nitrogen_containing;
    int _saturated_or_aromatic_heteroatom_containing;
    int _unsaturated_non_aromatic_carbon_only;
    int _unsaturated_non_aromatic_nitrogen_containing;
    int _unsaturated_non_aromatic_heteroatom_containing;

    int _rings_examined;

    int _last_ring_was_aromatic;

  // private functions

    void _default_values();

  public:
    Ring_Statistics();

    void reset() { _default_values(); }

    void extra(Molecule &m, const atomic_number_t *z, const Ring &r);

    void insert_results(int threshold, int *fp) const;

    int last_ring_was_aromatic() const { return _last_ring_was_aromatic; }
};

Ring_Statistics::Ring_Statistics() {
  _default_values();

  return;
}

void
Ring_Statistics::_default_values() {
  _saturated_or_aromatic_carbon_only = 0;
  _saturated_or_aromatic_nitrogen_containing = 0;
  _saturated_or_aromatic_heteroatom_containing = 0;
  _unsaturated_non_aromatic_carbon_only = 0;
  _unsaturated_non_aromatic_nitrogen_containing = 0;
  _unsaturated_non_aromatic_heteroatom_containing = 0;

  _rings_examined = 0;
  _last_ring_was_aromatic = 0;

  return;
}

void
Ring_Statistics::extra(Molecule &m, const atomic_number_t *z,
                       const Ring &r) {
  _rings_examined++;
  _last_ring_was_aromatic = 0;

  const int n = r.number_elements();

  int ncarbon = 0;
  int nnitrogen = 0;
  int fully_saturated = 1;

  static constexpr int kMaybeAromatic = -2;

  int aromatic;
  if (r.is_aromatic())
    aromatic = 1;
  else
    aromatic = kMaybeAromatic;

  for (int i = 0; i < n; i++) {
    atom_number_t j = r[i];

    if (6 == z[j])
      ncarbon++;
    else if (7 == z[j])
      nnitrogen++;

    if (kMaybeAromatic != aromatic) // already set
      ;
    else if (!m.is_aromatic(j))
      aromatic = 0;

    //  First check the atoms for being fully saturated. If they are all
    //  fully saturated, there will be no need to check the bonds

    if (0 == fully_saturated) // already turned off, no need to check
      ;
    else if (m.ncon(j) <
             m.nbonds(j)) // atom unsaturation found, will later check bonds
      fully_saturated = 0;
  }

  if (kMaybeAromatic == aromatic)
    aromatic = 1;

  // cerr << "Checking " << r << " arom " << aromatic << " saturated " <<
  // fully_saturated << endl;

  _last_ring_was_aromatic = aromatic;

  if (fully_saturated || aromatic) // no need to check bonds within the ring
    ;
  else // must check the within ring bonds
  {
    fully_saturated = 1; // assume OK
    for (Ring_Bond_Iterator i(r); i != r.zend(); i++) {
      atom_number_t ra1 = i.a1();
      //    cerr << "Iterator processing atoms " << ra1 << " and " << i.a2() <<
      //    endl;

      if (m.ncon(ra1) == m.nbonds(ra1))
        continue;

      atom_number_t ra2 = i.a2();

      if (m.ncon(ra1) == m.nbonds(ra1))
        continue;

      const Atom *a1 = m.atomi(ra1);

      const Bond *b = a1->bond_to_atom(ra1, ra2);

      if (nullptr == b) // dealing with a NON-SSSR ring
        continue;

      if (b->is_single_bond())
        continue;

      fully_saturated = 0; // found multiple bond, ring is unsaturated
      break;
    }
  }

  // cerr << "Ring of size " << n << " saturated " << fully_saturated << "
  // aromatic " << aromatic << " ncarbon " << ncarbon << endl;

  if (fully_saturated || aromatic) {
    if (ncarbon == n)
      _saturated_or_aromatic_carbon_only++;
    else {
      if (nnitrogen)
        _saturated_or_aromatic_nitrogen_containing++;
      _saturated_or_aromatic_heteroatom_containing++;
    }
  } else if (0 == fully_saturated && !aromatic) {
    if (n == ncarbon)
      _unsaturated_non_aromatic_carbon_only++;
    else {
      if (nnitrogen)
        _unsaturated_non_aromatic_nitrogen_containing++;
      _unsaturated_non_aromatic_heteroatom_containing++;
    }
  }

  return;
}

void
Ring_Statistics::insert_results(int threshold, int *fp) const {
  if (_rings_examined < threshold)
    return;

  fp[0]++;

  if (_saturated_or_aromatic_carbon_only >= threshold)
    fp[1]++;
  if (_saturated_or_aromatic_nitrogen_containing >= threshold)
    fp[2]++;
  if (_saturated_or_aromatic_heteroatom_containing >= threshold)
    fp[3]++;
  if (_unsaturated_non_aromatic_carbon_only >= threshold)
    fp[4]++;
  if (_unsaturated_non_aromatic_nitrogen_containing >= threshold)
    fp[5]++;
  if (_unsaturated_non_aromatic_heteroatom_containing >= threshold)
    fp[6]++;

  return;
}

#define LARGEST_RING_SIZE_PROCESSED 10

static Ring_Statistics ring_statistics[LARGEST_RING_SIZE_PROCESSED + 1];

static void
construct_walk_around_outside_strongly_fused(const Molecule &m,
                                             atom_number_t zatom,
                                             int *in_system,
                                             Set_of_Atoms &s) {
  int initial_flag = in_system[zatom]; // will be 1 or 2 (nrings value)

  in_system[zatom] = 0;

  const Atom *a = m.atomi(zatom);

  s.add(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++) {
    atom_number_t j = a->other(zatom, i);

    if (0 == in_system[j])
      continue;

    if (1 == initial_flag) // can go anywhere
      construct_walk_around_outside_strongly_fused(m, j, in_system, s);
    else if (2 == initial_flag && 1 == in_system[j])
      construct_walk_around_outside_strongly_fused(m, j, in_system, s);
  }

  return;
}

static void
construct_walk_around_outside(const Molecule &m,
                                          atom_number_t zatom, int *in_system,
                                          Set_of_Atoms &s) {
  in_system[zatom] = 0;

  const Atom *a = m.atomi(zatom);

  s.add(zatom);

  int acon = a->ncon();

  for (int i = 0; i < acon; i++) {
    atom_number_t j = a->other(zatom, i);

    if (in_system[j])
      construct_walk_around_outside(m, j, in_system, s);
  }

  return;
}

static void
process_fused_system(Molecule &m, const atomic_number_t *z,
                                 const Ring &r, int *tmp, int &aromatic_ring,
                                 int &hetero_aromatic_ring) {
#ifdef DEBUG_PROCESS_FUSED_SYSTEM
  cerr << "Examining fused ring " << r << endl;
#endif

  int matoms = m.natoms();

  int n = r.fused_ring_neighbours();

  // We deliberately restrict ourselves to pairs of rings

  for (int i = 0; i < n; i++) {
    const Ring *rj = r.fused_neighbour(i);

    //  cerr << "Ring j is " << (*rj) << endl;
    if (rj->ring_number() < r.ring_number())
      continue;

    int strongly_fused = 0;
    if (1 ==
        r.largest_number_of_bonds_shared_with_another_ring()) // definitely no
                                                              // strong fusions
                                                              // here
      ;
    else if (r.compute_bonds_shared_with(*rj) > 1)
      strongly_fused = 1;

    set_vector(tmp, matoms, 0);
    r.set_vector(tmp, 1);
    rj->increment_vector(tmp, 1);

    int atoms_in_both_rings = count_non_zero_occurrences_in_array(tmp, matoms);

#ifdef DEBUG_PROCESS_FUSED_SYSTEM
    cerr << "There are " << atoms_in_both_rings << " atoms in both: strong? "
         << strongly_fused << endl;
#endif

    if (atoms_in_both_rings > LARGEST_RING_SIZE_PROCESSED)
      continue;

    Ring s;
    if (strongly_fused)
      construct_walk_around_outside_strongly_fused(m, r[0], tmp, s);
    else
      construct_walk_around_outside(m, r[0], tmp, s);

    int ring_size = s.number_elements();

    ring_statistics[ring_size].extra(m, z, s);

    if (ring_statistics[ring_size].last_ring_was_aromatic()) {
      aromatic_ring++;
      if (contains_heteroatom(m, z, s))
        hetero_aromatic_ring++;
    }
  }

  return;
}

#define PUBCHEM_NBITS 888

static int
pubchem_fingerprints(Molecule &m, int *fp) {
  set_vector(fp, PUBCHEM_NBITS, 0);

  int *ecount = new_int(HIGHEST_ATOMIC_NUMBER + 1);
  std::unique_ptr<int[]> free_ecount(ecount);

  int matoms = m.natoms();

  atomic_number_t *z = new atomic_number_t[matoms];
  std::unique_ptr<atomic_number_t[]> free_z(z);

  m.atomic_numbers(z);

  for (int i = 0; i < matoms; i++) {
    const atomic_number_t zi = z[i];

    if (zi >= 0) {
      ecount[zi]++;
      ecount[1] += m.implicit_hydrogens(i);
    }
  }

  m.compute_aromaticity_if_needed();

  update_count_bits(ecount, fp);
  atom_pairs(m, z, fp);

  for (int i = 0; i <= LARGEST_RING_SIZE_PROCESSED; i++) {
    ring_statistics[i].reset();
  }

  m.compute_aromaticity_if_needed();

  int nr = m.nrings();

  int aromatic_ring = 0;
  int hetero_aromatic_ring = 0;

  for (int i = 0; i < nr; i++) {
    const Ring *ri = m.ringi(i);

    int ring_size = ri->number_elements();

    if (ring_size >= LARGEST_RING_SIZE_PROCESSED)
      continue;

    //  cerr << "Ezaminging ring " << (*ri) << endl;

    ring_statistics[ring_size].extra(m, z, *ri);

    if (ring_statistics[ring_size].last_ring_was_aromatic()) {
      aromatic_ring++;
      if (contains_heteroatom(m, z, *ri))
        hetero_aromatic_ring++;
    }
  }

  int *tmp = new int[matoms];
  std::unique_ptr<int[]> free_tmp(tmp);

  (void)m.ring_membership();

  for (int i = 0; i < nr; i++) {
    const Ring *ri = m.ringi(i);

    if (!ri->is_fused())
      continue;

    process_fused_system(m, z, *ri, tmp, aromatic_ring, hetero_aromatic_ring);
  }

  if (aromatic_ring >= 1)
    fp[255] = 1;
  if (hetero_aromatic_ring >= 1)
    fp[256] = 1;

  if (aromatic_ring >= 2)
    fp[257] = 1;
  if (hetero_aromatic_ring >= 2)
    fp[258] = 1;

  if (aromatic_ring >= 3)
    fp[259] = 1;
  if (hetero_aromatic_ring >= 3)
    fp[260] = 1;

  if (aromatic_ring >= 4)
    fp[261] = 1;
  if (hetero_aromatic_ring >= 4)
    fp[262] = 1;

  ring_statistics[3].insert_results(1, fp + 115);
  ring_statistics[3].insert_results(2, fp + 122);

  ring_statistics[4].insert_results(1, fp + 129);
  ring_statistics[4].insert_results(2, fp + 136);

  ring_statistics[5].insert_results(1, fp + 143);
  ring_statistics[5].insert_results(2, fp + 150);
  ring_statistics[5].insert_results(3, fp + 157);
  ring_statistics[5].insert_results(4, fp + 164);
  ring_statistics[5].insert_results(5, fp + 171);

  ring_statistics[6].insert_results(1, fp + 178);
  ring_statistics[6].insert_results(2, fp + 185);
  ring_statistics[6].insert_results(3, fp + 192);
  ring_statistics[6].insert_results(4, fp + 199);
  ring_statistics[6].insert_results(5, fp + 206);

  ring_statistics[7].insert_results(1, fp + 213);
  ring_statistics[7].insert_results(2, fp + 220);

  ring_statistics[8].insert_results(1, fp + 227);
  ring_statistics[8].insert_results(2, fp + 234);

  ring_statistics[9].insert_results(1, fp + 241);

  ring_statistics[10].insert_results(1, fp + 248);

  if (external_queries_present)
    do_external_queries(m, fp);

  return 1;
}

static int
write_any_subsets(Molecule &m, const int *fp,
                             IWString_and_File_Descriptor *write_bit_file) {
  for (int i = 0; i < PUBCHEM_NBITS; i++) {
    if (0 == fp[i])
      continue;

    if (!write_bit_file[i].is_open())
      continue;

    write_bit_file[i] << m.smiles() << ' ' << m.name() << '\n';
    write_bit_file[i].write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
do_write_array(Molecule &m, const int *fp,
                          IWString_and_File_Descriptor &output) {
  output << m.smiles() << ' ' << m.name() << '\n';

  for (int i = 0; i < PUBCHEM_NBITS; i++) {
    if (fp[i])
      output << i << " 1\n";
    else if (write_zero_bits)
      output << i << " 0\n";
  }

  output << "|\n";

  return 1;
}

// After each molecule is processed, we might need to flush
// `output`, or write if it gets too full.
void
MaybeFlush(IWString_and_File_Descriptor& output) {
  if (flush_after_every_molecule) {
    output.flush();
  } else {
    output.write_if_buffer_holds_more_than(32768);
  }
}

static int
pubchem_fingerprints(Molecule &m,
                                IWString_and_File_Descriptor &output) {
  int fp[PUBCHEM_NBITS + 1];

  if (!pubchem_fingerprints(m, fp)) {
    cerr << "Fatal error processing '" << m.name() << "'\n";
    return 0;
  }

  if (verbose) {
    int nset = count_non_zero_occurrences_in_array(fp, PUBCHEM_NBITS);
    bits_hit.extra(nset);
  }

  if (nullptr != write_bit_file)
    write_any_subsets(m, fp, write_bit_file);

  if (write_as_array)
    return do_write_array(m, fp, output);

  if (!function_as_filter) {
    m.remove_explicit_hydrogens();
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  IW_Bits_Base b;

  b.construct_from_array_of_ints(fp, PUBCHEM_NBITS);
  IWString tmp;
  b.daylight_ascii_representation_including_nset_info(tmp);

  output << fingerprint_tag << tmp << ">\n";

  if (!function_as_filter)
    output << "|\n";

  MaybeFlush(output);

  return 1;
}

static int
pubchem_fingerprints(data_source_and_type<Molecule> &input,
                                IWString_and_File_Descriptor &output) {
  Molecule *m;
  while (nullptr != (m = input.next_molecule())) {
    molecules_read++;

    std::unique_ptr<Molecule> free_m(m);

    preprocess(*m);

    if (!pubchem_fingerprints(*m, output))
      return 0;

    output.write_if_buffer_holds_more_than(32768);
  }

  return 1;
}

static int
pubchem_fingerprints_filter_record(const const_IWSubstring &smiles,
                                   IWString_and_File_Descriptor &output) {
  Molecule m;

  if (!m.build_from_smiles(smiles)) {
    cerr << "Invalid smiles '" << smiles << "'\n";
    return 0;
  }

  return pubchem_fingerprints(m, output);
}

static int
pubchem_fingerprints_filter(iwstring_data_source &input,
                                       IWString_and_File_Descriptor &output) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    output << buffer << '\n';

    output.write_if_buffer_holds_more_than(4096);

    if ("|" == buffer) {
      tdts_read++;
      continue;
    }

    if (!buffer.starts_with(smiles_tag))
      continue;

    buffer.remove_leading_chars(smiles_tag.length());
    buffer.chop();

    if (!pubchem_fingerprints_filter_record(buffer, output)) {
      cerr << "Fatal error processing smiles '" << buffer << "', line "
           << input.lines_read() << endl;
      return 0;
    }
  }

  return 1;
}

static int
pubchem_fingerprints_filter(const char *fname,
                                       IWString_and_File_Descriptor &output) {
  iwstring_data_source input(fname);

  if (!input.ok()) {
    cerr << "Cannot open filter source '" << fname << "'\n";
    return 0;
  }

  return pubchem_fingerprints_filter(input, output);
}

static int
pubchem_fingerprints(const char *fname, FileType input_type,
                                IWString_and_File_Descriptor &output) {
  assert(nullptr != fname);

  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 1)
    input.set_verbose(1);

  return pubchem_fingerprints(input, output);
}

static int
remove_non_matches(resizable_array_p<Substructure_Hit_Statistics> &q,
                              IW_STL_Hash_Set &only_process) {
  int rc = 0;

  if (0 == only_process.size())
    return 0;

  for (int i = q.number_elements() - 1; i >= 0; i--) {
    const IWString c = q[i]->comment();

    if (only_process.contains(c))
      continue;

    q[i]->set_comment("");
    rc++;
  }

  return rc;
}

static int
read_queries(const_IWSubstring &fname,
             resizable_array_p<Substructure_Hit_Statistics> &queries) {
  assert(fname.starts_with("section"));

  fname.remove_leading_chars(9);

  if (!smarts_from_file(fname, queries, verbose)) {
    cerr << "Cannot read queries from '" << fname << "'\n";
    return 0;
  }

  int nq = queries.number_elements();

  assert(nq > 0);

  for (int i = 0; i < nq; i++) {
    queries[i]->set_max_matches_to_find(1);
    queries[i]->set_save_matched_atoms(0);
  }

  return nq;
}

static int
initialise_write_bit_files(const Command_Line &cl, char flag) {
  write_bit_file = new IWString_and_File_Descriptor[PUBCHEM_NBITS];

  int i = 0;
  const_IWSubstring s;
  while (cl.value(flag, s, i++)) {
    int b;

    if (!s.numeric_value(b) || b < 0 || b > PUBCHEM_NBITS) {
      cerr << "Invalid bit number '" << s << "'\n";
      return 0;
    }

    IWString fname;
    fname << stream_for_bit_stem << b << ".smi";

    if (!write_bit_file[b].open(fname.null_terminated_chars())) {
      cerr << "Cannot open '" << fname << "'\n";
      return 0;
    }
  }

  return i;
}

static void
DisplayDashXOptions(std::ostream& output) {
  output << " -X flush            flush output after each molecule\n";

  ::exit(0);
}

static int
pubchem_fingerprints(int argc, char **argv) {
  Command_Line cl(argc, argv, "vA:E:i:g:lfazq:b:O:J:X:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, verbose, 'A')) {
      cerr << "Cannot initialise aromaticity specifications\n";
      usage(5);
    }
  } else
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('E')) {
    if (!process_elements(cl, verbose, 'E')) {
      cerr << "Cannot initialise elements\n";
      return 6;
    }
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1,
                                                              'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  initialise_pair_xref();

  if (cl.option_present('a')) {
    write_as_array = 1;

    if (verbose)
      cerr << "Will write as an array\n";

    if (cl.option_present('z'))
      write_zero_bits = 1;
  }

  if (cl.option_present('q')) {
    int i = 0;
    const_IWSubstring q;
    while (cl.value('q', q, i++)) {
      int rc;
      if (q.starts_with("section4"))
        rc = read_queries(q, section4);
      else if (q.starts_with("section5"))
        rc = read_queries(q, section5);
      else if (q.starts_with("section6"))
        rc = read_queries(q, section6);
      else if (q.starts_with("section7"))
        rc = read_queries(q, section7);
      else {
        cerr << "Unrecognised query specification '" << q << "'\n";
        return 3;
      }

      if (!rc) {
        cerr << "Cannot read queries '" << q << "'\n";
        return i + 1;
      }

      if (verbose)
        cerr << "Read " << rc << " sectional queries\n";
    }
    external_queries_present = 1;
  }

  if (cl.option_present('b')) {
    if (!initialise_write_bit_files(cl, 'b')) {
      cerr << "Cannot initialise bit subset file(s), -b option\n";
      usage(3);
    }
  }

  // Not sure this is useful.
  if (cl.option_present('O')) {
    IW_STL_Hash_Set only_process;
    int i = 0;
    IWString o;
    while (cl.value('O', o, i++)) {
      only_process.insert(o);
    }

    remove_non_matches(section4, only_process);
    remove_non_matches(section5, only_process);
    remove_non_matches(section6, only_process);
    remove_non_matches(section7, only_process);

    int queries_remaining = 0;

    queries_remaining = section4.number_elements() +
                        section5.number_elements() +
                        section6.number_elements() + section7.number_elements();

    if (0 == queries_remaining)
      cerr << "Warning, all external queries suppressed\n";
    else if (verbose)
      cerr << "Subset selection leaves " << queries_remaining << " queries\n";
  }

  if (cl.option_present('J')) {
    if (write_as_array) {
      cerr << "The -a and -J options are mutually inconsistent\n";
      usage(3);
    }

    cl.value('J', fingerprint_tag);

    if (verbose)
      cerr << "Fingerprints written with tag '" << fingerprint_tag << "'\n";

    if (!fingerprint_tag.ends_with('<'))
      fingerprint_tag << '<';
  }

  if (cl.option_present('X')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('X', x, i); ++i) {
      if (x == "flush") {
        flush_after_every_molecule = 1;
        if (verbose) {
          cerr << "Will flush output after every molecule\n";
        }
      } else if (x == "help") {
        DisplayDashXOptions(cerr);
      } else {
        cerr << "Unrecognised -X qualifier '" << x << "'\n";
        DisplayDashXOptions(cerr);
      }
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (function_as_filter)
    ;
  else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-"))
    input_type = FILE_TYPE_SMI;
  else if (!all_files_recognised_by_suffix(cl))
    return 4;

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  if (cl.option_present('f')) {
    function_as_filter = 1;

    if (!pubchem_fingerprints_filter(cl[0], output))
      rc = 2;
  } else {
    for (int i = 0; i < cl.number_elements(); i++) {
      if (!pubchem_fingerprints(cl[i], input_type, output)) {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules\n";

    cerr << "Fingerprints had between " << bits_hit.minval() << " and "
         << bits_hit.maxval() << " bits set, ave " << bits_hit.average()
         << endl;
  }

  if (nullptr != write_bit_file)
    delete[] write_bit_file;

  return rc;
}

int
main(int argc, char **argv) {
  prog_name = argv[0];

  int rc = pubchem_fingerprints(argc, argv);

  return rc;
}
