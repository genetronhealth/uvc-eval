
This repository contains the scripts and code used to evaluate the performance of UVC and other variant callers.
Filenames of the main evaluation scripts are matched with the regular expression ( */eval-*.sh ).
The raw database (mainly in the datafiles/ and tools/ directories), SRA, FASTQ, and/or BAM files must be manually downloaded by the user into the corresponding directories that are supposed to contain the input to the main evaluation scripts.
The reason is that the optimal method to download the raw data may depend on the geographical location, ISP, network bandwidth, network latency, etc. of the user,
  so the method to download data is left to the user.

In total, there are about 5 terabytes of compressed raw data. 
Moreover, even more data are derived from such compressed data for evaluation. 
In addition, several variant callers are evaluated.
Therefore, the whole evaluation can take a very long time to finish.

