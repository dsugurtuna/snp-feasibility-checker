# Why it's built this way

## The problem

Planning a recall-by-genotype study starts with two questions: is the SNP typed for enough people,
and how many people should be in each genotype group? In a cohort genotyped on several arrays, both
answers depend on which arrays carry the SNP and how many participants each array covers.

## Design choices

**Why track participants per array, not just SNPs per array?** Because "the SNP is on an array" is
not the same as "the SNP is typed for the cohort". A SNP on the smaller of two arrays is typed for a
fraction of people, and that fraction is the right N for any estimate.

**Why Hardy-Weinberg expectations at all?** Because at the planning stage you often have an allele
frequency from a reference population but no genotype counts yet. HWE turns that into an expected
count quickly. It is a planning estimate, and the README says so.

**Why round instead of truncate?** Because truncation after floating-point arithmetic loses whole
people at random. It was a real bug here: 395 expected carriers where the exact answer is 396.

**Why raise on a wrong manifest column?** Because the failure mode otherwise is an empty array and a
report that says nothing is available. A loud error costs a minute; a quiet wrong answer can cancel
a study.

**Why no dependencies?** Because this is set arithmetic and one formula. Anyone should be able to
read it and check it by hand.

## Questions worth asking

**"HWE assumes a lot. When would these numbers be badly wrong?"**
For rare variants (small expected counts are noisy), variants under selection, and cohorts with
mixed ancestry, where population structure produces fewer heterozygotes than HWE predicts. When the
cohort is already genotyped, observed counts from PLINK `--freq` or `--hardy` are better, which is
the first roadmap item.

**"Your expected carriers are not people I can recall. What is missing?"**
Consent status, withdrawal, eligibility criteria, whether contact details are current, and response
rate. Each is a separate multiplier with its own uncertainty. I keep them out of the genetic
expectation so that each can be seen and challenged on its own.

**"What happens with a participant typed on two arrays?"**
They are counted twice in `genotyped_samples`. Fixing that properly needs participant-level data,
not array-level counts, which is what
[ld-linkage-mapper](https://github.com/dsugurtuna/ld-linkage-mapper)'s participant mapper works
from.

## What's next

- Read observed genotype counts from PLINK output.
- Add per-SNP call-rate thresholds.
- Show attrition factors next to, not inside, the genetic expectation.
