This pipeline takes params.greeting, which defaults to the string Hello world!, and splits it into individual words in the SPLITLETTERS process. 
Each word is written to a separate file, named chunk_aa, chunk_ab, chunk_acand so on. These files are picked up as the process output.

The second process CONVERTTOUPPER takes the output channel from the first process as its input. The use of the operator .flatten() here is to
split the SPLITLETTERS output channel element that contains two files into two separate elements to be put through the CONVERTTOUPPERprocess,
else they would be treated as a single element. The CONVERTTOUPPER process thus launches two tasks, one for each element. 
The bash script uses cat to print the file contents and tr to convert to upper-case. It takes the resulting standard-out as the process output channel.
