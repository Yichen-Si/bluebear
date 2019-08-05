ibs0 & phasing
# bluebear

// Aug. 5. Notes:
— Existing problem in the unconditional approach:

1 Simple global criteria to execute flips may be difficult to decide.
1.1 Future flips are not aware of past flips but treat the result so far as the given haplotypes.

2 Many flips can occur at one position.
2.1 At each position, there is a weighted directional graph representing the switch & IBD relation.
Currently, the length of extension is not considered when deciding which individual to flip so it greedily flip according to net degree. It does not seem reasonable to decide by length since switches are not independent and there can by a single wrongly phased variants.
It can be slow when #sample increases.

3 Do we want to treat variants differently according to frequency?

4 Suffix pbwt does not go forward, so two dense int matrices has to be kept in memory when processing a chunk.


