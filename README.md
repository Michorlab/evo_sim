# evo_sim

This software package is a framework for simulating stochastic evolutionary processes. Currently it implements both branching process and Moran process simulations. However, this software can be modified to implement arbitrary birth-death processes, while retaining the same input/output and multithreading structures. Similarly, additional output data writers, mutational behaviors, reproduction strategies, inheritance models, and clone-associated data (such as barcodes) can be added to the software while retaining the functionality of the rest of the software. The strength of this package is its modularity and ease of customization.

## Command-line interface and file types
The command line call format is: evo_sim -i [input file path] -o [output file folder path] -m [simulation type] -n [number of threads]

All of the above command line inputs are required. The simulation type is currently either "branching" or "moran". If there is an error in the command line inputs, the program will print to the console and exit. If there is an error with the input file format, a message detailing the error will print to a file in the output directory with extension ".eevo".

Input text files have a format detailed below and are of file extension ".ievo". Output text files have formats that depend on what data they are recording, and have file extension ".oevo".

## Input file formatting
Individual lines in the input file are read as separate commands. These commands can be in any order, but one mistake in the format of any of the commands will result in an error. To introduce a comment line, begin the line with the pound sign ("#"). 

There are currently 5 general types of valid commands in an input file. These are listed below and identified by the initial string that begins that type of line. For specific instructions, look at the code and example input files.

1. sim_params commands. These are simulation parameters and include the number of trials and information on how mutations are handled. The num_simulations, mut_handler_type, and mut_handler_params parameters are required.
2. pop_params commands. These are cell population parameters. Currently, the death rate and maximum cell types parameters are required.
3. writer commands. These are optional and determine what data from the simulation will be written to output files.
4. listener commands. These are optional and determine what stopping conditions each simulation trial will have. Simulation trials will always stop when there are no cells left in the population.
5. clone commands. These determine what clones are present initially. At least one clone command is required.
