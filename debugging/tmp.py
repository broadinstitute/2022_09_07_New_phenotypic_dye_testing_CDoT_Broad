#!/usr/bin/env jupyter
"""
Script to understand why the Saguaro profiles are so different to the other ones.
Phase seems to be very important, as it is not present in the batch 1 cppipe.
"""

with open("batch1.txt", "r") as f:
    batch1 = [line.strip("\n") for line in f.readlines()]

with open("batch5.txt", "r") as f:
    batch5 = [line.strip("\n") for line in f.readlines()]

orig_channels = ["Brightfield", "DNA", "AGP", "Actin", "ER", "Mito", "RNA"]
sag_channels = ["Brightfield", "DNA", "_647", "CL488R", "CL488Y", "CL561YE"]

mapper = {k: v for k, v in zip(sag_channels, orig_channels)}

batch5_translated = []
for x in batch5:
    result = x
    for k, v in mapper.items():
        result = result.replace(k, v)
    batch5_translated.append(result)

set_batch5 = set(batch5_translated)
set_batch1 = set(batch1)

n_intersection = len(set_batch1.intersection(set_batch5))
# n_intersection
# Out: 6348
n_difference = len(set_batch1.difference(set_batch5))
# n_difference
# Out: 1333

batch1_noextra = [
    x
    for x in batch1
    if list(set(orig_channels).difference(mapper.values()))[0] not in x
]
# len(batch1_noextra)
# Out: 6545

# In [274]: len(batch1) - len(batch1_noextra)
# Out: 1136

in_saguaro_only = set(batch5_translated).difference(batch1_noextra)
# len( in_saguaro_only )
# Out: 1204

nocorr = sorted(
    ["_".join(x.split("_")[1:]) for x in in_saguaro_only if "Correlation" not in x]
)
# Check how many are not derived from correlation
# Correlation was changed between batches
# len(nocorr)
# Out: 958

# nophase = [x for x in nocorr if "Phase" not in x]
# In [426]: len(nophase)
# Out: 27
# Most of the difference in features comes from the phase channel, present in one but not the best
