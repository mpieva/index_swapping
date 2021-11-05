# Index swapping

## Synopsis

This repository provides a perl script for identifying cross-contamination 
events that occurred between double-index libraries due to index swapping 
during library amplification or sequencing:

```sh
cross_contamination.pl
```

The script includes a help menu with further instructions and explanations 
(accessible by executing the script without specifying an input file).

## Software Dependencies

The script has been successfully used with perl 5 (version 30) and has no 
further dependencies on perl modules or external software. 

## Example

In what follows we provide an example for how this script was used to determine 
cross-contamination using a demultiplexing report file formatted as described 
in the help menu of the script. The exemplary input file 
(`demultiplexing_report_B52992.txt`) is provided as part of the documentation. 

```sh
./cross_contamination.pl demultiplexing_report_B52992.txt >report_file.txt
```

This script generates text output on the shell (identified cross-contamination 
events: (i) by event, and (ii) by library), which can be redirected to a file 
as indicated above.
