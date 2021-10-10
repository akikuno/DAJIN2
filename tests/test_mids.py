import os
from src.DAJIN2.preprocess import mids

inputdir = os.path.join("tests", "data", "preprocess",
                        "mids", "samfile_to_mids", "input")
input_samfiles = [os.path.join(inputdir, _) for _ in os.listdir(inputdir)]

samfile = input_samfiles[-3]

print(mids.samfile_to_mids(samfile))
